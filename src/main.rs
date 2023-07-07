use std::collections::hash_set::HashSet;

mod hll;

fn main() {
    let s = std::fs::read_to_string("./shakespeare.txt")
        .expect("load file")
        .parse::<String>()
        .expect("parse file");
    let shakespeare = s.split_whitespace();

    let mut hll = hll::HyperLogLog::new().expect("hll construction");
    println!("bits: {}", hll.b);
    println!("m_counters: {}", hll.m_counters);

    // add all of our words
    for word in shakespeare.clone() {
        hll.add(word);
    }

    let estimated_count = hll.count() as f64;
    let actual_count = shakespeare
        .into_iter()
        .collect::<HashSet<_>>()
        .iter()
        .count() as f64;

    let estimated_error = hll.error();
    let error = (actual_count - estimated_count).abs() / (actual_count);

    println!("Estimated {estimated_count:.0} err({estimated_error:.5})");
    println!("Actual {actual_count:.0} err({error:.5})");
}
