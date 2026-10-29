use std::fs;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::Command;

#[test]
fn reads_from_stdin() {
    let input = "chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1\nchr1\t.\ttranscript\t1\t10\t.\t+\t.\tParent=gene1;ID=tx1\n";

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_gffsort-rs"));
    cmd.arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());

    let mut child = cmd.spawn().unwrap();
    child.stdin.as_mut().unwrap().write_all(input.as_bytes()).unwrap();
    drop(child.stdin.take());

    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1"));
    assert!(stdout.contains("chr1\t.\ttranscript\t1\t10\t.\t+\t.\tParent=gene1;ID=tx1"));
}

#[test]
fn sorts_naturally_and_extracts_fasta() {
    let dir = std::env::temp_dir().join(format!("gffsort-rs-{}", std::process::id()));
    let _ = fs::create_dir_all(&dir);
    let path = dir.join("sample.gff3");
    let content = "chr2\t.\tgene\t1\t10\t.\t+\t.\tID=gene2\nchr10\t.\tgene\t1\t10\t.\t+\t.\tID=gene10\nchr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1\n##FASTA\n>chr1\nACGT\n";
    fs::write(&path, content).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_gffsort-rs"))
        .arg("--chr_order")
        .arg("natural")
        .arg("--extract_FASTA")
        .arg(&path)
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.starts_with("chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1"));
    assert!(stdout.contains("chr2\t.\tgene\t1\t10\t.\t+\t.\tID=gene2"));

    let fasta_path = dir.join("sample.fasta");
    assert!(fasta_path.exists());
    let fasta = fs::read_to_string(fasta_path).unwrap();
    assert!(fasta.contains(">chr1"));
    assert!(fasta.contains("ACGT"));

    let _ = fs::remove_file(path);
    let _ = fs::remove_file(dir.join("sample.fasta"));
    let _ = fs::remove_dir(dir);
}

#[test]
fn precise_mode_orders_parents_before_children() {
    let input = "chr1\t.\ttranscript\t1\t10\t.\t+\t.\tID=tx1;Parent=gene1\nchr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1\nchr1\t.\texon\t1\t10\t.\t+\t.\tID=exon1;Parent=tx1\n";

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_gffsort-rs"));
    cmd.arg("--precise")
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());

    let mut child = cmd.spawn().unwrap();
    child.stdin.as_mut().unwrap().write_all(input.as_bytes()).unwrap();
    drop(child.stdin.take());

    let output = child.wait_with_output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.find("ID=gene1").unwrap() < stdout.find("ID=tx1").unwrap());
    assert!(stdout.find("ID=tx1").unwrap() < stdout.find("ID=exon1").unwrap());
}

