use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};

use clap::Parser;

#[derive(Debug, Parser)]
#[command(
    name = "gffsort-rs",
    version,
    about = "Sort GFF3 files for tabix indexing",
    long_about = "Sort GFF3 files for tabix indexing. Outputs to stdout by default."
)]
struct Cli {
    #[arg(value_name = "INPUT", help = "Input GFF3 file; use '-' or omit for stdin")]
    input: Option<PathBuf>,

    #[arg(long, default_value = "alphabet", help = "Chromosome order: alphabet, natural, or original")]
    chr_order: String,

    #[arg(long, action = clap::ArgAction::SetTrue, help = "Run in precise mode")]
    precise: bool,

    #[arg(long, action = clap::ArgAction::SetTrue, help = "Extract inline FASTA sequences to a .fasta file")]
    extract_fasta: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ChrOrder {
    Alphabet,
    Natural,
    Original,
}

#[derive(Debug, Default)]
struct ParsedGff {
    chromosomes: Vec<String>,
    records_by_chr: HashMap<String, HashMap<u64, Vec<String>>>,
    fasta: String,
    pragma_lines: Vec<String>,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let cli = Cli::parse();
    let chr_order = parse_chr_order(&cli.chr_order)?;

    let input = match cli.input.as_deref() {
        Some(path) if path == Path::new("-") => None,
        Some(path) => Some(path.to_path_buf()),
        None => None,
    };

    let input_data = read_input(input.as_deref())?;
    let mut parsed = parse_gff(input_data.as_deref(), &cli, input.as_deref())?;

    let mut chromosome_names = parsed.chromosomes;
    match chr_order {
        ChrOrder::Alphabet => chromosome_names.sort(),
        ChrOrder::Natural => chromosome_names.sort_by(|a, b| natural_cmp(a, b)),
        ChrOrder::Original => {}
    }

    let stdout = io::stdout();
    let mut out = stdout.lock();

    for line in &parsed.pragma_lines {
        writeln!(out, "{line}").map_err(|e| format!("Failed writing output: {e}"))?;
    }

    for chr in chromosome_names {
        let Some(positions) = parsed.records_by_chr.get(&chr) else {
            continue;
        };

        let mut sorted_positions: Vec<u64> = positions.keys().copied().collect();
        sorted_positions.sort_unstable();

        for pos in sorted_positions {
            let Some(lines) = positions.get(&pos) else {
                continue;
            };

            if lines.len() == 1 {
                writeln!(out, "{}", lines[0]).map_err(|e| format!("Failed writing output: {e}"))?;
                continue;
            }

            if cli.precise {
                let sorted_lines = topo_sort_same_position(lines);
                for line in sorted_lines {
                    writeln!(out, "{line}").map_err(|e| format!("Failed writing output: {e}"))?;
                }
            } else {
                for line in lines {
                    writeln!(out, "{line}").map_err(|e| format!("Failed writing output: {e}"))?;
                }
            }
        }
    }

    Ok(())
}

fn parse_chr_order(value: &str) -> Result<ChrOrder, String> {
    match value.to_ascii_lowercase().as_str() {
        "alphabet" => Ok(ChrOrder::Alphabet),
        "natural" => Ok(ChrOrder::Natural),
        "original" => Ok(ChrOrder::Original),
        other => Err(format!(
            "Unknown option: --chr_order={}. Only [alphabet] [natural] [original] are allowed",
            other
        )),
    }
}

fn read_input(input_path: Option<&Path>) -> Result<Option<String>, String> {
    match input_path {
        Some(path) if path == Path::new("-") => {
            let mut buffer = String::new();
            io::stdin()
                .lock()
                .read_to_string(&mut buffer)
                .map_err(|err| format!("Failed to read stdin: {err}"))?;
            Ok(Some(buffer))
        }
        Some(path) => {
            fs::read_to_string(path)
                .map(Some)
                .map_err(|err| format!("Failed to read input file '{}': {err}", path.display()))
        }
        None => {
            let mut buffer = String::new();
            io::stdin()
                .lock()
                .read_to_string(&mut buffer)
                .map_err(|err| format!("Failed to read stdin: {err}"))?;
            Ok(Some(buffer))
        }
    }
}

fn parse_gff(input: Option<&str>, cli: &Cli, input_path: Option<&Path>) -> Result<ParsedGff, String> {
    let mut parsed = ParsedGff::default();
    let mut encountered_fasta = false;
    let input = input.unwrap_or("");

    for raw_line in input.lines() {
        let line = raw_line.trim_end_matches('\r');

        if line == "##FASTA" {
            encountered_fasta = true;
            continue;
        }

        if encountered_fasta {
            if cli.extract_fasta && line != "##FASTA" {
                parsed.fasta.push_str(line);
                parsed.fasta.push('\n');
            }
            continue;
        }

        if line.starts_with('#') {
            if line.chars().all(|ch| ch == '#') {
                continue;
            }
            parsed.pragma_lines.push(line.to_string());
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let chr = fields[0].to_string();
        let pos = match fields[3].parse::<u64>() {
            Ok(v) => v,
            Err(_) => continue,
        };

        let chromosome_entry = parsed.records_by_chr.entry(chr.clone()).or_default();
        chromosome_entry.entry(pos).or_default().push(line.to_string());

        if !parsed.chromosomes.contains(&chr) {
            parsed.chromosomes.push(chr);
        }
    }

    if cli.extract_fasta && !parsed.fasta.trim().is_empty() {
        let fasta_path = fasta_output_path(input_path);
        fs::write(&fasta_path, &parsed.fasta)
            .map_err(|err| format!("Failed to write extracted FASTA to '{}': {err}", fasta_path.display()))?;
        eprintln!("The inline FASTA sequences were extracted to file {}", fasta_path.display());
    }

    Ok(parsed)
}

fn fasta_output_path(input_path: Option<&Path>) -> PathBuf {
    match input_path {
        Some(path) if path != Path::new("-") => {
            let mut out = path.to_path_buf();
            out.set_extension("fasta");
            out
        }
        _ => PathBuf::from("stdin.fasta"),
    }
}

fn topo_sort_same_position(lines: &[String]) -> Vec<String> {
    let mut id_to_line: HashMap<String, String> = HashMap::new();
    let mut parent_to_children: HashMap<String, Vec<String>> = HashMap::new();
    let mut node_order: Vec<String> = Vec::new();

    for line in lines {
        let attrs = line.rsplit_once('\t').map_or("", |(_, attrs)| attrs);
        let attributes = parse_attributes(attrs);

        let id = attributes.get("ID").and_then(|values| values.first().cloned());
        let parents = attributes.get("Parent").cloned().unwrap_or_default();

        let key = id.clone().unwrap_or_else(|| line.clone());
        if !id_to_line.contains_key(&key) {
            node_order.push(key.clone());
        }
        id_to_line.insert(key.clone(), line.clone());

        for parent in parents {
            parent_to_children.entry(parent.clone()).or_default().push(key.clone());
        }
    }

    let sorted_ids = topo_sort(&node_order, &parent_to_children);
    let mut result = Vec::with_capacity(sorted_ids.len());
    for id in sorted_ids {
        if let Some(line) = id_to_line.get(&id) {
            result.push(line.clone());
        }
    }
    result
}

fn topo_sort(node_order: &[String], parent_to_children: &HashMap<String, Vec<String>>) -> Vec<String> {
    let mut indegree: HashMap<String, usize> = HashMap::new();
    let mut children_map: HashMap<String, Vec<String>> = HashMap::new();
    let mut all_nodes: HashSet<String> = HashSet::new();

    for node in node_order {
        indegree.insert(node.clone(), 0);
        all_nodes.insert(node.clone());
    }

    for (parent, children) in parent_to_children {
        for child in children {
            all_nodes.insert(child.clone());
            if !indegree.contains_key(child) {
                indegree.insert(child.clone(), 0);
            }

            children_map.entry(parent.clone()).or_default().push(child.clone());
            *indegree.get_mut(child).unwrap() += 1;
        }
    }

    let mut queue: VecDeque<String> = VecDeque::new();
    for node in &all_nodes {
        if indegree.get(node).copied().unwrap_or(0) == 0 {
            queue.push_back(node.clone());
        }
    }

    let mut result = Vec::new();
    while let Some(node) = queue.pop_front() {
        result.push(node.clone());

        if let Some(children) = children_map.get(&node) {
            for child in children {
                let deg = indegree.get_mut(child).unwrap();
                *deg = deg.saturating_sub(1);
                if *deg == 0 {
                    queue.push_back(child.clone());
                }
            }
        }
    }

    if result.len() != all_nodes.len() {
        let mut remaining: Vec<String> = all_nodes.iter().cloned().collect();
        remaining.sort();
        for node in remaining {
            if !result.contains(&node) {
                result.push(node);
            }
        }
    }

    result
}

fn parse_attributes(attributes: &str) -> HashMap<String, Vec<String>> {
    let mut result: HashMap<String, Vec<String>> = HashMap::new();
    for part in attributes.split(';') {
        let trimmed = part.trim();
        if trimmed.is_empty() {
            continue;
        }

        let Some((key, value)) = trimmed.split_once('=') else {
            continue;
        };

        let key = key.trim();
        let value = value.trim();
        if key.is_empty() || value.is_empty() {
            continue;
        }

        let entries: Vec<String> = value
            .split(',')
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(str::to_string)
            .collect();

        if entries.is_empty() {
            continue;
        }

        result.insert(key.to_string(), entries);
    }
    result
}

fn natural_cmp(a: &str, b: &str) -> Ordering {
    let a_parts = split_natural(a);
    let b_parts = split_natural(b);

    for (left, right) in a_parts.iter().zip(b_parts.iter()) {
        match left.cmp(right) {
            Ordering::Equal => {}
            other => return other,
        }
    }

    a_parts.len().cmp(&b_parts.len())
}

fn split_natural(value: &str) -> Vec<NaturalChunk> {
    let mut parts = Vec::new();
    let mut current = String::new();
    let mut last_was_digit = None;

    for ch in value.chars() {
        let is_digit = ch.is_ascii_digit();
        if last_was_digit == Some(is_digit) {
            current.push(ch);
        } else {
            if !current.is_empty() {
                parts.push(parse_chunk(&current));
            }
            current.clear();
            current.push(ch);
            last_was_digit = Some(is_digit);
        }
    }

    if !current.is_empty() {
        parts.push(parse_chunk(&current));
    }

    parts
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum NaturalChunk {
    Text(String),
    Number(u64),
}

impl PartialOrd for NaturalChunk {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for NaturalChunk {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (NaturalChunk::Text(left), NaturalChunk::Text(right)) => {
                left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase())
            }
            (NaturalChunk::Number(left), NaturalChunk::Number(right)) => left.cmp(right),
            (NaturalChunk::Text(_), NaturalChunk::Number(_)) => Ordering::Greater,
            (NaturalChunk::Number(_), NaturalChunk::Text(_)) => Ordering::Less,
        }
    }
}

fn parse_chunk(value: &str) -> NaturalChunk {
    if value.chars().all(|ch| ch.is_ascii_digit()) {
        NaturalChunk::Number(value.parse::<u64>().unwrap_or(0))
    } else {
        NaturalChunk::Text(value.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_order_is_used_for_chromosomes() {
        let mut values = vec!["Chr7".to_string(), "Chr1".to_string(), "Chr10".to_string(), "Chr2".to_string()];
        values.sort_by(|a, b| natural_cmp(a, b));
        assert_eq!(values, vec!["Chr1", "Chr2", "Chr7", "Chr10"]);
    }

    #[test]
    fn precise_mode_keeps_parents_before_children() {
        let lines = vec![
            "chr1\t.\tgene\t1\t10\t.\t+\t.\tID=gene1".to_string(),
            "chr1\t.\ttranscript\t1\t10\t.\t+\t.\tParent=gene1;ID=tx1".to_string(),
            "chr1\t.\texon\t1\t10\t.\t+\t.\tParent=tx1;ID=exon1".to_string(),
        ];

        let sorted = topo_sort_same_position(&lines);
        assert_eq!(sorted[0], lines[0]);
        assert_eq!(sorted[1], lines[1]);
        assert_eq!(sorted[2], lines[2]);
    }

    #[test]
    fn parses_gff_attributes_with_spaces() {
        let attrs = "ID=gene1; Name=gene one; Parent=gene1,other";
        let parsed = parse_attributes(attrs);
        assert_eq!(parsed.get("Name").unwrap()[0], "gene one");
        assert_eq!(parsed.get("Parent").unwrap(), &vec!["gene1", "other"]);
    }
}
