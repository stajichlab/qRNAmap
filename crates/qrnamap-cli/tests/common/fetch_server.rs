//! Local scripted HTTP server for tests. Copied as is from `crates/genome-fetch/tests/common/mod.rs` (ruling R6: no shared dev crate).
//!
//! `Server::start(handler)` binds 127.0.0.1 on a free port and serves each request with `handler(request, own_base_url)`.
//! The server records every request (path and headers, header names in lower case) for later checks.
#![allow(dead_code)]

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::{Arc, Mutex};
use std::thread;

#[derive(Debug, Clone)]
pub struct Request {
    pub path: String,
    /// (lower-case name, value)
    pub headers: Vec<(String, String)>,
}

impl Request {
    pub fn header(&self, name: &str) -> Option<&str> {
        let name = name.to_ascii_lowercase();
        self.headers
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
    }
}

#[derive(Debug, Clone)]
pub struct Response {
    pub status: u16,
    pub headers: Vec<(String, String)>,
    pub body: Vec<u8>,
}

impl Response {
    pub fn ok(body: impl Into<Vec<u8>>) -> Response {
        Response {
            status: 200,
            headers: vec![],
            body: body.into(),
        }
    }
    pub fn status(status: u16) -> Response {
        Response {
            status,
            headers: vec![],
            body: vec![],
        }
    }
    /// Close the connection without sending a response (a transport failure for the client).
    pub fn close() -> Response {
        Response {
            status: 0,
            headers: vec![],
            body: vec![],
        }
    }
    pub fn redirect(location: &str) -> Response {
        Response {
            status: 302,
            headers: vec![("Location".into(), location.into())],
            body: vec![],
        }
    }
}

pub struct Server {
    /// `http://127.0.0.1:<port>`
    pub base: String,
    pub port: u16,
    requests: Arc<Mutex<Vec<Request>>>,
}

impl Server {
    pub fn start<F>(handler: F) -> Server
    where
        F: Fn(&Request, &str) -> Response + Send + 'static,
    {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
        let port = listener.local_addr().unwrap().port();
        let base = format!("http://127.0.0.1:{port}");
        let requests = Arc::new(Mutex::new(Vec::new()));
        let (rec, own) = (requests.clone(), base.clone());
        thread::spawn(move || {
            for stream in listener.incoming() {
                let Ok(mut stream) = stream else { continue };
                let Some(req) = read_request(&mut stream) else {
                    continue;
                };
                rec.lock().unwrap().push(req.clone());
                let resp = handler(&req, &own);
                if resp.status == 0 {
                    continue; // status 0 means: close the connection without any response
                }
                let _ = write_response(&mut stream, &resp);
            }
        });
        Server {
            base,
            port,
            requests,
        }
    }

    /// All requests seen so far.
    pub fn requests(&self) -> Vec<Request> {
        self.requests.lock().unwrap().clone()
    }
}

fn read_request(stream: &mut TcpStream) -> Option<Request> {
    let mut buf = Vec::new();
    let mut byte = [0u8; 1];
    while !buf.ends_with(b"\r\n\r\n") {
        if stream.read(&mut byte).ok()? == 0 {
            return None;
        }
        buf.push(byte[0]);
    }
    let text = String::from_utf8_lossy(&buf).to_string();
    let mut lines = text.split("\r\n");
    let path = lines.next()?.split_whitespace().nth(1)?.to_string();
    let headers = lines
        .filter_map(|l| l.split_once(':'))
        .map(|(n, v)| (n.trim().to_ascii_lowercase(), v.trim().to_string()))
        .collect();
    Some(Request { path, headers })
}

fn write_response(stream: &mut TcpStream, r: &Response) -> std::io::Result<()> {
    let mut head = format!(
        "HTTP/1.1 {} X\r\nContent-Length: {}\r\nConnection: close\r\n",
        r.status,
        r.body.len()
    );
    for (n, v) in &r.headers {
        head.push_str(&format!("{n}: {v}\r\n"));
    }
    head.push_str("\r\n");
    stream.write_all(head.as_bytes())?;
    stream.write_all(&r.body)?;
    stream.flush()
}
