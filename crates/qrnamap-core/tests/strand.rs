//! Strandedness report tests. Inputs are hand-written counts.

use qrnamap_core::{assigned_fraction, parse_featurecounts_summary, Strandedness, STRAND_THRESHOLD};

#[test]
fn forward_library_suggests_forward() {
    let r = assigned_fraction((950, 50, 10));
    assert_eq!(r.suggested, Strandedness::Forward);
    assert!((r.sense_fraction - 0.95).abs() < 1e-12);
    assert!((r.antisense_fraction - 0.05).abs() < 1e-12);
    assert!((r.ambiguous_fraction - 10.0 / 1010.0).abs() < 1e-12);
    assert!(!r.low_fraction_warning);
}

#[test]
fn reverse_library_suggests_reverse() {
    let r = assigned_fraction((30, 970, 0));
    assert_eq!(r.suggested, Strandedness::Reverse);
    assert!(!r.low_fraction_warning);
}

#[test]
fn unstranded_library_warns_low_fraction() {
    let r = assigned_fraction((510, 490, 5));
    assert_eq!(r.suggested, Strandedness::Unstranded);
    assert!(r.low_fraction_warning);
    assert!((r.sense_fraction - 0.51).abs() < 1e-12);
}

#[test]
fn threshold_boundary_is_not_stranded() {
    // Exactly 0.8 does not exceed the threshold.
    assert_eq!(STRAND_THRESHOLD, 0.8);
    let r = assigned_fraction((800, 200, 0));
    assert_eq!(r.suggested, Strandedness::Unstranded);
    assert!(!r.low_fraction_warning, "0.8 is not below 0.8");
    let r = assigned_fraction((801, 199, 0));
    assert_eq!(r.suggested, Strandedness::Forward);
}

#[test]
fn no_assigned_reads_is_unstranded_with_warning() {
    let r = assigned_fraction((0, 0, 7));
    assert_eq!(r.suggested, Strandedness::Unstranded);
    assert!(r.low_fraction_warning);
    assert_eq!(r.sense_fraction, 0.0);
    assert_eq!(r.antisense_fraction, 0.0);
    assert_eq!(r.ambiguous_fraction, 1.0);
    let r = assigned_fraction((0, 0, 0));
    assert_eq!(r.ambiguous_fraction, 0.0);
}

#[test]
fn featurecounts_summary_parses() {
    let text = "Status\t/tmp/x.sam\nAssigned\t9500\nUnassigned_Unmapped\t300\r\nUnassigned_Ambiguity\t12\n\n";
    let m = parse_featurecounts_summary(text).unwrap();
    assert_eq!(m["Assigned"], 9500);
    assert_eq!(m["Unassigned_Unmapped"], 300);
    assert_eq!(m["Unassigned_Ambiguity"], 12);
    assert_eq!(m.len(), 3);
}

#[test]
fn featurecounts_summary_rejects_bad_input() {
    assert!(parse_featurecounts_summary("").is_err());
    assert!(parse_featurecounts_summary("Assigned\t5\n").is_err());
    assert!(parse_featurecounts_summary("Status\tx\nAssigned\tfive\n").is_err());
}
