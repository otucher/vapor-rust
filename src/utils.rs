use std::cmp;
use itertools::Itertools; // 0.9.0


pub fn get_gaps(weights: &Vec<usize>) -> Vec<(usize, usize)> {
    let mut gaps = Vec::new();

    // add last non-zero value to capture gap for weights that end in 0
    let mut weights_plus_one = weights.to_owned();
    weights_plus_one.push(1);

    for (i, weight) in weights_plus_one.iter().skip(1).enumerate() {
        // i is behind 1 due to skip
        let current_is_zero = *weight == 0;
        let previous_is_zero = weights[i] == 0;
        if current_is_zero != previous_is_zero {
            gaps.push(i + 1)
        }
    }
    gaps.into_iter().tuples().collect()
}


pub fn merge_overlapping_intervals(intervals: &mut Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if intervals.is_empty() {
        intervals.to_vec()
    } else {
        intervals.sort_by_key(|x| x.0);
        let mut results = vec![intervals[0]];
        for interval in intervals.iter().skip(1) {
            let j = results.len() - 1;
            if interval.0 >= results[j].0 && interval.0 <= results[j].1 {
                results[j].1 = cmp::max(interval.1, results[j].0);
            } else {
                results.push(*interval);
            }
        }
        results
    }
}
