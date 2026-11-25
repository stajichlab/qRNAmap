//! Build script for `qrnamap-cli`.
//!
//! `serve.rs` embeds the web assets of `dist-site/` with `rust-embed`. The derive is a macro, so
//! cargo does not see that folder as an input of the crate: after the site is rebuilt, cargo
//! would keep the binary of the earlier build, with the assets of that build. This script watches
//! the folder and stamps the crate with a hash of its contents, so a changed site gives a new
//! stamp, a new fingerprint and a rebuild that embeds the new assets.

use std::collections::hash_map::DefaultHasher;
use std::fs;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

/// The folder with the built web assets, the same one `serve.rs` embeds.
const ASSETS: &str = "../../dist-site";

fn main() {
    let assets = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap()).join(ASSETS);

    let mut files = Vec::new();
    collect(&assets, &assets, &mut files);
    files.sort();

    // Cargo watches the folder itself and every file in it. The folder alone is not enough: cargo
    // compares the modification time of the folder, which does not change when a file in it is
    // only rewritten. A folder that does not exist keeps the build script fresh on every build,
    // but the stamp below does not change, so the crate is not rebuilt for nothing.
    println!("cargo:rerun-if-changed={}", assets.display());
    for (name, _, _) in &files {
        println!("cargo:rerun-if-changed={}", assets.join(name).display());
    }

    println!("cargo:rustc-env=QRNAMAP_ASSETS_STAMP={}", stamp(&files));
}

/// A hash of the names, sizes and modification times of the files under `dir`. A missing folder
/// gives the hash of no file.
fn stamp(files: &[(String, u64, u64)]) -> String {
    let mut hasher = DefaultHasher::new();
    for (name, size, mtime) in files {
        name.hash(&mut hasher);
        size.hash(&mut hasher);
        mtime.hash(&mut hasher);
    }
    format!("{:016x}", hasher.finish())
}

/// Collect the name, size and modification time of every file under `dir`. The name is relative
/// to `root` and uses `/` as the separator, so the stamp does not depend on the platform.
fn collect(root: &Path, dir: &Path, out: &mut Vec<(String, u64, u64)>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let Ok(meta) = entry.metadata() else {
            continue;
        };
        if meta.is_dir() {
            collect(root, &path, out);
        } else if meta.is_file() {
            let name = path
                .strip_prefix(root)
                .unwrap_or(&path)
                .to_string_lossy()
                .replace('\\', "/");
            let mtime = meta
                .modified()
                .ok()
                .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
                .map(|d| d.as_nanos() as u64)
                .unwrap_or(0);
            out.push((name, meta.len(), mtime));
        }
    }
}
