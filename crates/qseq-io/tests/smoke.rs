#[test]
fn crate_builds() {
    assert_eq!(env!("CARGO_PKG_NAME"), "qseq-io");
}
