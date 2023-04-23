use std::cmp;
use std::env;
use std::path::PathBuf;
use std::collections::HashMap;
use bio::io::{fasta, fastq};
use itertools::Itertools;
use glob;
use std::fmt;

mod utils;


#[derive(Debug)]
struct ExtendBridgeResult {
    maskable_walk_seq: Vec<u8>,
    scores: Vec<usize>,
    total_score: usize,
}

#[derive(Debug, Clone)]
struct SearchResult {
    name: String,
    kmers: Vec<Vec<u8>>,
    weights: Vec<usize>,
    gaps: Vec<(usize, usize)>
}

impl fmt::Display for SearchResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // let kmers_as_string = self.kmers
        //     .iter()
        //     .map(|kmer| String::from_utf8(kmer.to_owned()).expect("ruh roh"))
        //     .collect_vec();
        // write!(
        //     f,
        //     "SearchResult {{ name: {}, kmers: {:?}, weights: {:?}, gaps: {:?} }}",
        //     self.name, kmers_as_string, self.weights, self.gaps
        // )
        write!(
            f,
            "SearchResult {{ name: {} }}",
            self.name
        )
    }
}


impl SearchResult {
    fn score(&self) -> usize {
        self.weights.iter().sum()
    }

    fn coverage(&self) -> f64 {
        let non_zero = self.weights.iter().filter(|x| **x != 0).collect_vec();
        non_zero.len() as f64 / self.weights.len() as f64
    }

    fn total_score(&self) -> f64 {
        self.score() as f64 * self.coverage()
    }

    fn sequence(&self) -> Vec<u8> {
        self.kmers.concat()
    }
}


#[derive(Debug)]
struct Graph {
    nodes: HashMap<Vec<u8>, usize>,
    kmer_size: usize,
}


impl Graph {
    fn from(fastq_files: &Vec<PathBuf>, kmer_size: usize, node_min_coverage: usize) -> Self {
        let reads = fastq_files
            .iter()
            .map(|file| fastq::Reader::from_file(file)
                .expect("fastq ruh roh")
                .records()
                .map(|record| record.expect("Error during fasta parsing").seq().to_vec())
                .collect_vec()
            )
            .concat();
        println!("num reads {}", reads.len());

        let kmers = reads
            .iter()
            .map(|read| read
                .windows(kmer_size)
                .map(|kmer| kmer.to_vec())
                .collect_vec()
            )
            .concat();
        println!("num kmers {}", kmers.len());

        // convert sequence into kmers, get unique counts, and filter minimmum coverage
        let mut nodes = kmers.into_iter().counts();
        println!("num nodes {}", nodes.len());

        // filter read coverage
        nodes = nodes
            .into_iter()
            .filter(|(_kmer, count)| count >= &node_min_coverage)
            .collect();
        println!("num nodes after filtering {}", nodes.len());

        Self {
            nodes,
            kmer_size
        }
    }

    fn get_weights(&self, kmers: &Vec<Vec<u8>>) -> Vec<usize> {
        let mut weights = Vec::new();
        for kmer in kmers {
            let weight = match self.nodes.get(kmer) {
                Some(x) => *x,
                _ => 0
            };
            weights.push(weight)
        }
        weights
    }

    fn get_search_result(&self, name: String, kmers: Vec<Vec<u8>>) -> SearchResult {
        // let kmers = self.get_kmers(&sequence);
        let weights = self.get_weights(&kmers);
        SearchResult {
            name,
            gaps: utils::get_gaps(&weights),
            weights,
            kmers,
        }
    }

    fn get_alt_scores(&self, kmer: &Vec<u8>, reverse: bool) -> Vec<&usize> {
        let nucleotides = [b"A", b"C", b"G", b"T"];
        let kmer = match reverse {
            false => kmer.to_owned(),
            true => kmer.to_owned().into_iter().rev().collect_vec()
        };
        nucleotides
            .iter()
            .map(|nt| {
                let mut node = kmer[..self.kmer_size-1].to_vec();
                node.push(nt[0]);
                node = match reverse {
                    false => node,
                    true => node.into_iter().rev().collect()
                };
                self.nodes.get(&node).unwrap_or(&0)
            })
            .collect_vec()
    }

    fn get_suboptimal_branches(&self, result: &SearchResult, reverse: bool) -> Vec<usize> {
        result.kmers
            .iter()
            .zip(result.weights.iter())
            .enumerate()
            .filter(|(_idx, (kmer, _weight))| self.nodes.contains_key(*kmer))
            .filter(|(_idx, (kmer, weight))| {
                let scores = self.get_alt_scores(&kmer, reverse);
                match scores.iter().max() {
                    Some(max_score) => max_score > weight,
                    _ => false
                }
            })
            .map(|(idx, _)| match reverse {
                true => idx + 1,
                false => idx,
            })
            .collect_vec()
    }

    fn extend_gaps(&self, mut search_result: SearchResult) -> SearchResult {
        let max_trim_size = self.kmer_size + 1;
        let forward_suboptimal_branches = self.get_suboptimal_branches(&search_result, false);
        let reverse_suboptimal_branches = self.get_suboptimal_branches(&search_result, true);
        let mut new_gaps = Vec::new();
        for (start, end) in search_result.gaps.iter() {
            let mut new_gap = (start.clone(), end.clone());
            for new_start in 0..cmp::max(0, start.checked_sub(max_trim_size).unwrap_or(0)) {
                if forward_suboptimal_branches.contains(&new_start) {
                    new_gap.0 = new_start;
                    break;
                }
            }

            for new_end in 0..cmp::min(end + max_trim_size, search_result.sequence().len()) {
                if reverse_suboptimal_branches.contains(&new_end) {
                    new_gap.1 = new_end;
                    break;
                }
            }
            new_gaps.push(new_gap)
        }
        search_result.gaps = utils::merge_overlapping_intervals(&mut new_gaps);
        search_result
    }

    fn extend_bridge(&self, start_kmer: &Vec<u8>, gap_size: usize, reverse: bool) -> ExtendBridgeResult {
        let nucleotides = [b"A", b"C", b"G", b"T"];
        let mut walk_seq = match reverse {
            false => start_kmer.to_owned(),
            true => start_kmer.to_owned().into_iter().rev().collect()
        };
        let mut scores: Vec<usize> = Vec::new();
        while walk_seq.len() < (gap_size + start_kmer.len()) {
            let walk_kmer = walk_seq[walk_seq.len()-(self.kmer_size-1)..].to_vec();
            let alt_scores = self.get_alt_scores(&walk_kmer, reverse);
            let max_score = alt_scores.iter().max().expect("ruh roh");
            if max_score == &&0 {
                break;
            } else {
                let idx = alt_scores
                    .iter()
                    .position(|x| x == max_score)
                    .expect("ruh roh");
                scores.push(**max_score);
                walk_seq.push(nucleotides[idx][0]);
            }
        }

        // trim and mask
        let mask = (0..gap_size.checked_sub(walk_seq.len()).unwrap_or(0))
            .map(|_| b"-"[0])
            .collect_vec();
        walk_seq = walk_seq[self.kmer_size..].to_vec();
        walk_seq.extend(mask);
        println!("{:?}", walk_seq);
        let maskable_walk_seq = walk_seq[..walk_seq.len()-self.kmer_size+1].to_vec();

        // reverse
        match reverse {
            true => {
                scores = scores.into_iter().rev().collect();
                ExtendBridgeResult {
                    maskable_walk_seq: maskable_walk_seq.into_iter().rev().collect(),
                    total_score: scores.iter().sum(),
                    scores,
                }
            },
            false => {
                ExtendBridgeResult {
                    maskable_walk_seq,
                    total_score: scores.iter().sum(),
                    scores
                }
            }
        }
    }

    fn spread_weights(&self, mut search_result: SearchResult) -> SearchResult {
        let mut new_weights = vec![0, 0, 0];
        new_weights.extend(&search_result.weights);
        new_weights.extend(vec![0, 0, 0]);
        search_result.weights = new_weights
            .windows(self.kmer_size)
            .map(|x| *x.iter().max().expect("ruh roh"))
            .collect();
        search_result
    }

    fn graph_walk(&self, mut search_result: SearchResult) -> SearchResult {
        let mut all_masks = Vec::new();
        search_result = self.extend_gaps(search_result);
        for (start, end) in search_result.gaps.iter() {
            let gap_size = end.checked_sub(*start).unwrap_or(0);

            // get bridges and choose maximum total score
            let mut bridges = Vec::new();
            if start > &0 {
                let start_kmer = &search_result.kmers[start - 1];
                let forward_bridge = self.extend_bridge(start_kmer, gap_size, false);
                bridges.push(forward_bridge);
            }
            if end < &search_result.kmers.len() {
                let start_kmer = &search_result.kmers[*end];
                let reverse_bridge = self.extend_bridge(start_kmer, gap_size, true);
                bridges.push(reverse_bridge);
            }
            let max_bridge = bridges
                .iter()
                .rev()
                .max_by_key(|x| x.total_score)
                .unwrap_or(&bridges[0]);

            // define substring (possibly empty) of gap subject to masking.
            let gap_seq = &search_result.sequence()[*start..*end];
            let maskable_ref_seq = &gap_seq[self.kmer_size-1..];
            for idx in 0..maskable_ref_seq.len() {
                if maskable_ref_seq[idx] != max_bridge.maskable_walk_seq[idx] {
                    let val = start + self.kmer_size - 1 + idx;
                    all_masks.push(val);
                }
            }

            // update weights for entire gap (not just the maskable portion).
            search_result.weights.splice(*start..*end, max_bridge.scores.clone());
        }
        search_result = self.spread_weights(search_result);
        for mask in all_masks.iter() {
            search_result.weights[*mask] = 0
        }
        search_result
    }

    fn select_refs(&self, reference_fasta_path: PathBuf, min_ref_coverage: f64, use_top_candidates: usize) -> Vec<SearchResult> {
        fasta::Reader::from_file(reference_fasta_path)
            .expect("Unable to read reference FASTA")
            .records()
            .map(|record| {
                let record = record.expect("Unable to access FASTA record");
                let name = record.id().to_string();
                let kmers = record
                    .seq()
                    .windows(self.kmer_size)
                    .map(|kmer| kmer.to_vec())
                    .collect_vec();
                self.get_search_result(name, kmers)
            })
            .filter(|x| x.coverage() > min_ref_coverage)
            .sorted_by(|a, b| a
                .coverage()
                .partial_cmp(&b.coverage())
                .expect("Unable to sort Results by coverage score")
            )
            .enumerate()
            .filter(|(idx, _x)| match use_top_candidates {
                0 => true,
                val => idx < &val
            })
            .map(|(_, x)| self.graph_walk(x))
            .sorted_by(|a, b|
                a.total_score().partial_cmp(&b.total_score()).expect("Unable to sort Search Results by total score"))
            .collect_vec()
    }
}


fn main() {
    let fastq_files = glob::glob("/Users/oliver.tucher/Code/Scratch/rust-practice/*.fastq")
        .expect("Unable to fetch FASTQ files")
        .map(|file| file.expect("Unable to open FASTQ file"))
        .collect_vec();
    println!("fastq files {:?}", fastq_files);

    let graph = Graph::from(
        &fastq_files,
        4,
        1
    );
    let reference_fasta_path = env::current_dir().unwrap().join("panviral_steam_reference_20221112.fasta");
    let results = graph.select_refs(reference_fasta_path, 0.0, 0);
    println!("{}", results[0]);
}
