"""Tests of compare_parity.py with tiny hand-made pairs: equal passes, a difference fails.

Run: python3 -m pytest -q validation/browser-spike0/test_compare_parity.py
"""
import gzip
import importlib.util
import sys
from pathlib import Path

spec = importlib.util.spec_from_file_location("compare_parity", Path(__file__).with_name("compare_parity.py"))
cp = importlib.util.module_from_spec(spec)
sys.modules["compare_parity"] = cp
spec.loader.exec_module(cp)

ABUND = ("target_id\tlength\teff_length\test_counts\ttpm\n"
         "t1\t1000\t800\t10\t5000\n"
         "t2\t2000\t1800\t20\t5000\n")
SAM = ("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chrI\tLN:230218\n"
       "@PG\tID:hisat2\tPN:hisat2\tCL:\"hisat2-align-s -x a\"\n"
       "r1\t99\tchrI\t100\t60\t10M\t=\t200\t110\tACGTACGTAC\tIIIIIIIIII\n"
       "r1\t147\tchrI\t200\t60\t10M\t=\t100\t-110\tACGTACGTAC\tIIIIIIIIII\n")
FC = ("# Program:featureCounts v2.1.1; Command:\"featureCounts\" \"-o\" \"fc.tsv\"\n"
      "Geneid\tChr\tStart\tEnd\tStrand\tLength\tyeast.sam\n"
      "g1\tchrI\t1\t100\t+\t100\t7\n"
      "g2\tchrI\t200\t300\t-\t101\t0\n")
FCS = "Status\tyeast.sam\nAssigned\t7\nUnassigned_NoFeatures\t1\n"
DE = ("\"gene\",\"baseMean\",\"log2FoldChange\",\"lfcSE\",\"stat\",\"pvalue\",\"padj\"\n"
      "\"g1\",100.5,1.25,0.2,6.25,4.1e-10,8.2e-10\n"
      "\"g2\",0,NA,NA,NA,NA,NA\n"
      "\"g3\",50,-0.5,0.3,-1.6667,0.0956,0.0956\n")


def w(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def test_kallisto_equal_and_unequal(tmp_path):
    a = w(tmp_path / "a.tsv", ABUND)
    assert cp.compare_kallisto(a, w(tmp_path / "b.tsv", ABUND))[0]
    ok, detail = cp.compare_kallisto(a, w(tmp_path / "c.tsv", ABUND.replace("\t20\t5000", "\t20.0001\t5000")))
    assert not ok and "est_counts" in detail
    ok, detail = cp.compare_kallisto(a, w(tmp_path / "d.tsv", ABUND.replace("t2\t2000\t1800\t20\t5000", "t2\t2000\t1800\t20\t5001")))
    assert not ok and "tpm" in detail
    # eff_length is not one of the checked columns
    assert cp.compare_kallisto(a, w(tmp_path / "e.tsv", ABUND.replace("\t800\t", "\t801\t")))[0]


def test_sam_ignores_pg_only(tmp_path):
    a = w(tmp_path / "a.sam", SAM)
    assert cp.compare_sam(a, w(tmp_path / "b.sam.gz", SAM.replace("-x a", "-x other/path")))[0]
    ok, detail = cp.compare_sam(a, w(tmp_path / "c.sam", SAM.replace("\t200\t60\t10M", "\t201\t60\t10M")))
    assert not ok and "line 4" in detail  # line number after @PG removal
    ok, _ = cp.compare_sam(a, w(tmp_path / "d.sam", SAM + "r2\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"))
    assert not ok


def test_featurecounts_ignores_command_line(tmp_path):
    a, s = w(tmp_path / "a.tsv", FC), w(tmp_path / "a.summary", FCS)
    b = w(tmp_path / "b.tsv", FC.replace("\"fc.tsv\"", "\"/x/fc.tsv\""))
    assert cp.compare_featurecounts(a, b, s, w(tmp_path / "b.summary", FCS))[0]
    c = w(tmp_path / "c.tsv", FC.replace("\t101\t0", "\t101\t1"))
    ok, detail = cp.compare_featurecounts(a, c, s, s)
    assert not ok and "counts line" in detail
    ok, detail = cp.compare_featurecounts(a, b, s, w(tmp_path / "c.summary", FCS.replace("\t1\n", "\t2\n")))
    assert not ok and "summary" in detail


def test_deseq2_tolerance(tmp_path):
    a = w(tmp_path / "a.csv", DE)
    ok, detail = cp.compare_deseq2(a, w(tmp_path / "b.csv", DE))
    assert ok and "max_abs=0" in detail
    # within relative 1e-6
    assert cp.compare_deseq2(a, w(tmp_path / "c.csv", DE.replace("100.5,", "100.50000001,")))[0]
    # outside tolerance
    ok, detail = cp.compare_deseq2(a, w(tmp_path / "d.csv", DE.replace("1.25,", "1.2501,")))
    assert not ok and "log2FoldChange" in detail
    # NA pattern change
    ok, detail = cp.compare_deseq2(a, w(tmp_path / "e.csv", DE.replace("NA,NA,NA,NA,NA", "0,NA,NA,NA,NA")))
    assert not ok and "NA pattern" in detail


def test_evaluate_reports_failed_step(tmp_path):
    art = tmp_path / "parity" / "spike0-parity-firefox"
    w(art / "kallisto.steps.tsv",
      "engine\ttool\tstep\tstatus\texit_code\tseconds\targs\tdetail\toutputs\n"
      "firefox\tkallisto\tindex\terror\t\t1.0\tindex\tout of memory\t\n")
    row = cp.evaluate(tmp_path / "native", tmp_path / "parity", "firefox", "kallisto")
    assert row["result"] == "error" and "out of memory" in row["detail"]
    row = cp.evaluate(tmp_path / "native", tmp_path / "parity", "webkit", "hisat2")
    assert row["result"] == "error" and "not found" in row["detail"]
