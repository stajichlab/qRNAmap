use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::{Path, PathBuf};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const ZSTD_MAGIC: [u8; 4] = [0x28, 0xb5, 0x2f, 0xfd];

#[derive(Debug)]
pub enum InputError {
    /// An I/O failure. `path` is empty when built through `From<io::Error>`.
    Io {
        path: PathBuf,
        source: io::Error,
    },
    Truncated {
        path: PathBuf,
    },
    Corrupt {
        path: PathBuf,
        detail: String,
    },
    /// A 0-byte file whose name ends in `.gz` or `.zst`: not valid compressed data.
    Empty {
        path: PathBuf,
    },
}

impl fmt::Display for InputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InputError::Io { path, source } if path.as_os_str().is_empty() => {
                write!(f, "input error: {source}")
            }
            InputError::Io { path, source } => {
                write!(f, "input error: {}: {source}", path.display())
            }
            InputError::Truncated { path } => {
                write!(f, "compressed input is truncated: {}", path.display())
            }
            InputError::Corrupt { path, detail } => {
                write!(f, "corrupt compressed input: {}: {detail}", path.display())
            }
            InputError::Empty { path } => {
                write!(f, "compressed input is empty (0 bytes): {}", path.display())
            }
        }
    }
}

impl std::error::Error for InputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            InputError::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl From<io::Error> for InputError {
    fn from(e: io::Error) -> Self {
        InputError::Io {
            path: PathBuf::new(),
            source: e,
        }
    }
}

/// Compression of a file, from its magic bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Plain,
    Gzip,
    Zstd,
}

/// Detect the format from the magic bytes. The magic bytes win over the file
/// extension, so a `.txt` file that holds gzip data is decoded, and a `.gz`
/// file that holds plain text (or nothing) is read as plain.
fn detect(head: &[u8]) -> Format {
    if head.starts_with(&GZIP_MAGIC) {
        Format::Gzip
    } else if head.starts_with(&ZSTD_MAGIC) {
        Format::Zstd
    } else {
        Format::Plain
    }
}

fn truncated(path: &Path) -> io::Error {
    io::Error::new(
        io::ErrorKind::UnexpectedEof,
        InputError::Truncated {
            path: path.to_path_buf(),
        },
    )
}

fn corrupt(path: &Path, detail: String) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        InputError::Corrupt {
            path: path.to_path_buf(),
            detail,
        },
    )
}

/// Map a decoder error: EOF is truncation; data errors are corruption;
/// other kinds (real device errors) pass through.
fn map_decoder_error(path: &Path, e: io::Error) -> io::Error {
    match e.kind() {
        io::ErrorKind::UnexpectedEof => truncated(path),
        io::ErrorKind::InvalidInput | io::ErrorKind::InvalidData | io::ErrorKind::Other => {
            corrupt(path, e.to_string())
        }
        _ => e,
    }
}

/// zstd decoder with error mapping.
struct ZstdMap<R: BufRead> {
    inner: zstd::Decoder<'static, R>,
    path: PathBuf,
}

impl<R: BufRead> Read for ZstdMap<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner
            .read(buf)
            .map_err(|e| map_decoder_error(&self.path, e))
    }
}

enum GzState<R: BufRead> {
    Active(flate2::bufread::GzDecoder<R>),
    Between(R),
    Done,
    /// A read failed. Every later read returns the same error again (R29-8), so a caller that
    /// reads on after an error never sees a clean end of file.
    Failed(Failure),
}

/// What a failed read reported, kept so it can be returned again.
#[derive(Clone)]
enum Failure {
    Truncated,
    Corrupt(String),
    Other(io::ErrorKind, String),
}

impl Failure {
    fn of(e: &io::Error) -> Failure {
        match e.get_ref().and_then(|i| i.downcast_ref::<InputError>()) {
            Some(InputError::Truncated { .. }) => Failure::Truncated,
            Some(InputError::Corrupt { detail, .. }) => Failure::Corrupt(detail.clone()),
            _ => Failure::Other(e.kind(), e.to_string()),
        }
    }

    fn error(&self, path: &Path) -> io::Error {
        match self {
            Failure::Truncated => truncated(path),
            Failure::Corrupt(d) => corrupt(path, d.clone()),
            Failure::Other(k, m) => io::Error::new(*k, m.clone()),
        }
    }
}

/// Multi-member gzip reader. It handles members one by one so that an EOF in
/// the header of a later member (trailing bytes) is told apart from an EOF in
/// a payload (truncation).
struct GzMap<R: BufRead> {
    state: GzState<R>,
    members: usize,
    path: PathBuf,
}

impl<R: BufRead> GzMap<R> {
    fn new(r: R, path: PathBuf) -> Self {
        GzMap {
            state: GzState::Between(r),
            members: 0,
            path,
        }
    }
}

impl<R: BufRead> Read for GzMap<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if let GzState::Failed(f) = &self.state {
            return Err(f.error(&self.path));
        }
        if buf.is_empty() {
            return Ok(0);
        }
        let r = self.read_inner(buf);
        if let Err(e) = &r {
            self.state = GzState::Failed(Failure::of(e));
        }
        r
    }
}

impl<R: BufRead> GzMap<R> {
    fn read_inner(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        loop {
            match std::mem::replace(&mut self.state, GzState::Done) {
                GzState::Done | GzState::Failed(_) => return Ok(0),
                GzState::Between(mut r) => {
                    if r.fill_buf()?.is_empty() {
                        return Ok(0);
                    }
                    self.members += 1;
                    self.state = GzState::Active(flate2::bufread::GzDecoder::new(r));
                }
                GzState::Active(mut d) => match d.read(buf) {
                    Ok(0) => self.state = GzState::Between(d.into_inner()),
                    Ok(n) => {
                        self.state = GzState::Active(d);
                        return Ok(n);
                    }
                    Err(e) => {
                        if d.header().is_none() && self.members > 1 {
                            return Err(corrupt(
                                &self.path,
                                "trailing bytes after last member".to_string(),
                            ));
                        }
                        return Err(map_decoder_error(&self.path, e));
                    }
                },
            }
        }
    }
}

/// Read up to 4 bytes, looping over short reads, until 4 bytes or EOF.
fn read_head<R: Read>(r: &mut R) -> io::Result<Vec<u8>> {
    let mut head = Vec::with_capacity(4);
    let mut tmp = [0u8; 4];
    while head.len() < 4 {
        match r.read(&mut tmp[..4 - head.len()]) {
            Ok(0) => break,
            Ok(n) => head.extend_from_slice(&tmp[..n]),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(head)
}

/// The compression of the file at `path`, from its first bytes (not its name).
pub fn file_format(path: &Path) -> Result<Format, InputError> {
    let mut f = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let head = read_head(&mut f).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(detect(&head))
}

/// Open a plain, gzip (multi-member) or zstd file for buffered reading.
/// Decoding is lazy: a truncated stream is reported while reading, as an
/// `io::Error` of kind `UnexpectedEof` that wraps `InputError::Truncated`.
/// A corrupt stream is an `InvalidData` error that wraps `InputError::Corrupt`.
pub fn open_input(path: &Path) -> Result<Box<dyn BufRead>, InputError> {
    let file = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    open_reader(file, path)
}

/// Same as `open_input` for any reader. `path` is used in error messages.
/// The first bytes are read in a loop, so a reader that returns short reads
/// is still detected correctly.
pub fn open_reader<R: Read + 'static>(
    mut r: R,
    path: &Path,
) -> Result<Box<dyn BufRead>, InputError> {
    let io_err = |source| InputError::Io {
        path: path.to_path_buf(),
        source,
    };
    let head = read_head(&mut r).map_err(io_err)?;
    let name = path.to_string_lossy();
    if head.is_empty() && (name.ends_with(".gz") || name.ends_with(".zst")) {
        return Err(InputError::Empty {
            path: path.to_path_buf(),
        });
    }
    let format = detect(&head);
    let stream = BufReader::new(Cursor::new(head).chain(r));
    let path_buf = path.to_path_buf();
    Ok(match format {
        Format::Plain => Box::new(stream),
        Format::Gzip => Box::new(BufReader::new(GzMap::new(stream, path_buf))),
        Format::Zstd => Box::new(BufReader::new(ZstdMap {
            inner: zstd::Decoder::with_buffer(stream).map_err(io_err)?,
            path: path_buf,
        })),
    })
}
