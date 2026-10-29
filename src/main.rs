use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::env;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;

const USAGE: &str = r#"gffsort-rs [OPTIONS] input.file.gff3 > output.sort.gff3

Options:
  --help
  --precise           Run in precise mode, keeping parent features before children when they share a start position.
  --chr_order VALUE   Select chromosome ordering: alphabet, natural, original [default: alphabet]
  --extract_FASTA     Extract inline FASTA sequences to a sibling .fasta file.
"#;

#[derive(Clone, Copy, Debug)]
enum ChrOrder {
    Alphabet,
    Natural,
    Original,
}

#[derive(Debug)]
struct Options {
    help: bool,
    precise: bool,
    chr_order: ChrOrder,
    extract_fasta: bool,
}

#[derive(Default)]
struct ParsedGff {
    chromosomes: Vec<String>,
    records_by_chr: HashMap<String, HashMap<u64, Vec<String>>>,
    fasta: String,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let mut args = env::args().skip(1);
    let mut options = Options {
        help: false,
        precise: false,
        chr_order: ChrOrder::Alphabet,
        extract_fasta: false,
    };
    let mut input: Option<PathBuf> = None;

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--help" => options.help = true,
            "--precise" => options.precise = true,
            "--extract_FASTA" => options.extract_fasta = true,
            "--chr_order" => {
                let value = args.next().ok_or_else(|| "Missing value for --chr_order".to_string())?;
                let order = match value.to_ascii_lowercase().as_str() {
                    "alphabet" => ChrOrder::Alphabet,
                    "natural" => ChrOrder::Natural,
                    "original" => ChrOrder::Original,
                    other => {
                        return Err(format!(
                            "Unknown option: --chr_order={}. Only [alphabet] [natural] [original] are allowed",
                            other
                        ));
                    }
                };
                options.chr_order = order;
            }
            other if other.starts_with("--") => {
                return Err(format!("Unknown option: {other}"));
            }
            other => {
                if input.is_some() {
                    return Err("Expected a single input file path".to_string());
                }
                input = Some(PathBuf::from(other));
            }
        }
    }

    if options.help {
        print!("{USAGE}");
        return Ok(());
    }

    let input_path = input.ok_or_else(|| "Missing input file".to_string())?;
    let parsed = parse_gff(&input_path, &options)?;

    let mut chromosome_names = parsed.chromosomes;
    match options.chr_order {
        ChrOrder::Alphabet => chromosome_names.sort(),
        ChrOrder::Natural => chromosome_names.sort_by(|a, b| natural_cmp(a, b)),
        ChrOrder::Original => {}
    }

    let mut output = io::stdout().lock();
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
                writeln!(output, "{}", lines[0]).map_err(|e| format!("Failed writing output: {e}"))?;
                continue;
            }

            if options.precise {
                let sorted_lines = topo_sort_same_position(lines);
                for line in sorted_lines {
                    writeln!(output, "{line}").map_err(|e| format!("Failed writing output: {e}"))?;
                }
            } else {
                for line in lines {
                    writeln!(output, "{line}").map_err(|e| format!("Failed writing output: {e}"))?;
                }
            }
        }
    }

    Ok(())
}

fn parse_gff(input_path: &PathBuf, options: &Options) -> Result<ParsedGff, String> {
    let mut parsed = ParsedGff::default();
    let mut encountered_fasta = false;

    let contents = fs::read_to_string(input_path)
        .map_err(|e| format!("Failed to read input file '{}': {e}", input_path.display()))?;

    for raw_line in contents.lines() {
        let line = raw_line.trim_end_matches('\r');

        if line == "##FASTA" {
            encountered_fasta = true;
            continue;
        }

        if encountered_fasta {
            if options.extract_fasta && line != "##FASTA" {
                parsed.fasta.push_str(line);
                parsed.fasta.push('\n');
            }
            continue;
        }

        if line.starts_with('#') {
            if line.chars().all(|ch| ch == '#') {
                continue;
            }
            println!("{line}");
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
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

    if options.extract_fasta && !parsed.fasta.trim().is_empty() {
        let fasta_path = format!("{}{}.fasta", input_path.to_string_lossy(), "");
        fs::write(&fasta_path, &parsed.fasta)
            .map_err(|e| format!("Failed to write extracted FASTA to '{}': {e}", fasta_path))?;
        eprintln!("The inline FASTA sequences were extracted to file {}", fasta_path);
    }

    Ok(parsed)
}

fn topo_sort_same_position(lines: &[String]) -> Vec<String> {
    let mut id_to_line: HashMap<String, String> = HashMap::new();
    let mut parent_to_children: HashMap<String, Vec<String>> = HashMap::new();
    let mut node_order: Vec<String> = Vec::new();

    for line in lines {
        let attrs = line.rsplit_once('\t').map_or("", |(_, attrs)| attrs);
        let id = extract_attribute(attrs, "ID");
        let parents = extract_attribute_values(attrs, "Parent");

        let key = if let Some(ref id_value) = id {
            id_value.clone()
        } else {
            line.clone()
        };

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
    let mut child_to_parents: HashMap<String, Vec<String>> = HashMap::new();
    let mut indegree: HashMap<String, usize> = HashMap::new();
    let mut nodes: Vec<String> = Vec::new();

    for node in node_order {
        indegree.insert(node.clone(), 0);
        nodes.push(node.clone());
    }

    for (parent, children) in parent_to_children {
        for child in children {
            if !indegree.contains_key(child) {
                indegree.insert(child.clone(), 0);
                nodes.push(child.clone());
            }
            child_to_parents
                .entry(child.clone())
                .or_default()
                .push(parent.clone());
            *indegree.get_mut(child).unwrap() += 1;
        }
    }

    let mut queue = VecDeque::new();
    for node in &nodes {
        if indegree.get(node) == Some(&0) {
            queue.push_back(node.clone());
        }
    }

    let mut result = Vec::new();
    while let Some(node) = queue.pop_front() {
        result.push(node.clone());
        if let Some(children) = parent_to_children.get(&node) {
            for child in children {
                let deg = indegree.get_mut(child).unwrap();
                *deg = deg.saturating_sub(1);
                if *deg == 0 {
                    queue.push_back(child.clone());
                }
            }
        }
    }

    if result.len() != nodes.len() {
        for node in &nodes {
            if !result.contains(node) {
                result.push(node.clone());
            }
        }
    }

    result
}

fn extract_attribute(attrs: &str, key: &str) -> Option<String> {
    let prefix = format!("{key}=");
    for part in attrs.split(';') {
        if let Some(value) = part.strip_prefix(&prefix) {
            return Some(value.trim().to_string());
        }
    }
    None
}

fn extract_attribute_values(attrs: &str, key: &str) -> Vec<String> {
    let mut values = Vec::new();
    let prefix = format!("{key}=");
    for part in attrs.split(';') {
        if let Some(value) = part.strip_prefix(&prefix) {
            for candidate in value.split(',') {
                let cleaned = candidate.trim();
                if !cleaned.is_empty() {
                    values.push(cleaned.to_string());
                }
            }
        }
    }
    values
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
    let mut result = Vec::new();
    let mut current = String::new();
    let mut last_is_digit = None;

    for ch in value.chars() {
        let is_digit = ch.is_ascii_digit();
        if last_is_digit == Some(is_digit) {
            current.push(ch);
        } else {
            if !current.is_empty() {
                result.push(parse_chunk(&current));
            }
            current.clear();
            current.push(ch);
            last_is_digit = Some(is_digit);
        }
    }

    if !current.is_empty() {
        result.push(parse_chunk(&current));
    }

    result
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
            (NaturalChunk::Text(left), NaturalChunk::Text(right)) => left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase()),
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
    fn precise_mode_respects_parent_before_child() {
        let lines = vec![
            "chr1\t.	gene\t1\t10\t.\t+\t.\tID=gene1".to_string(),
            "chr1\t.\ttranscript\t1\t10\t.\t+\t.\tParent=gene1;ID=tx1".to_string(),
            "chr1\t.\texon\t1\t10\t.\t+\t.\tParent=tx1;ID=exon1".to_string(),
        ];

        let sorted = topo_sort_same_position(&lines);
        assert_eq!(sorted[0], lines[0]);
        assert_eq!(sorted[1], lines[1]);
        assert_eq!(sorted[2], lines[2]);
    }
}

