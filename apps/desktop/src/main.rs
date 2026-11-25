// Prevents additional console window on Windows in release
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use tauri::{WebviewUrl, WebviewWindowBuilder};

#[tauri::command]
fn open_url(url: String) -> Result<(), String> {
    if url.starts_with("http://") || url.starts_with("https://") {
        open::that(&url).map_err(|e| e.to_string())
    } else {
        Err("disallowed scheme".to_string())
    }
}

fn main() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![open_url])
        .setup(|app| {
            #[allow(unused_variables)]
            let window = WebviewWindowBuilder::new(app, "main", WebviewUrl::default())
                .title("qRNAmap")
                .inner_size(1100.0, 820.0)
                .min_inner_size(800.0, 600.0)
                .resizable(true)
                .on_navigation(|url| {
                    if url.scheme() == "http" || url.scheme() == "https" {
                        let _ = open::that(url.as_str());
                        return false;
                    }
                    true
                })
                .on_new_window(|url, _features| {
                    if url.scheme() == "http" || url.scheme() == "https" {
                        let _ = open::that(url.as_str());
                    }
                    tauri::webview::NewWindowResponse::Deny
                })
                .build()?;

            #[cfg(debug_assertions)]
            {
                let _ = window.open_devtools();
            }
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running qRNAmap desktop application");
}
