use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const ZSTD_MAGIC: [u8; 4] = [0x28, 0xb5, 0x2f, 0xfd];

#[derive(Debug)]
pub enum InputError {
    Io(io::Error),
    Truncated { path: PathBuf },
}

impl fmt::Display for InputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InputError::Io(e) => write!(f, "input error: {e}"),
            InputError::Truncated { path } => {
                write!(f, "compressed input is truncated: {}", path.display())
            }
        }
    }
}

impl std::error::Error for InputError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            InputError::Io(e) => Some(e),
            InputError::Truncated { .. } => None,
        }
    }
}

impl From<io::Error> for InputError {
    fn from(e: io::Error) -> Self {
        InputError::Io(e)
    }
}

#[derive(Debug, PartialEq, Eq)]
enum Format {
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

/// Wraps a decoder and maps its end-of-file errors to `InputError::Truncated`.
struct TruncationMap<R> {
    inner: R,
    path: PathBuf,
}

impl<R: Read> Read for TruncationMap<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf).map_err(|e| {
            if e.kind() == io::ErrorKind::UnexpectedEof {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    InputError::Truncated {
                        path: self.path.clone(),
                    },
                )
            } else {
                e
            }
        })
    }
}

/// Open a plain, gzip (multi-member) or zstd file for buffered reading.
/// Decoding is lazy: a truncated stream is reported while reading, as an
/// `io::Error` of kind `UnexpectedEof` that wraps `InputError::Truncated`.
pub fn open_input(path: &Path) -> Result<Box<dyn BufRead>, InputError> {
    let mut file = BufReader::new(File::open(path)?);
    let format = detect(file.fill_buf()?);
    let path_buf = path.to_path_buf();
    Ok(match format {
        Format::Plain => Box::new(file),
        Format::Gzip => Box::new(BufReader::new(TruncationMap {
            inner: flate2::bufread::MultiGzDecoder::new(file),
            path: path_buf,
        })),
        Format::Zstd => Box::new(BufReader::new(TruncationMap {
            inner: zstd::Decoder::with_buffer(file)?,
            path: path_buf,
        })),
    })
}
