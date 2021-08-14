use clap::*;
use indicatif::ProgressBar;
use rand::{seq::SliceRandom, thread_rng};
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::io::{self, BufRead};

#[derive(Debug)]
struct Snp {
    scaffold: String,
    pos: i64,

    name: String,
    a_ref: u8,
    a_alt: u8,
    ps: Vec<String>,
}

fn incompatible(snp1: &Snp, snp2: &Snp, threshold: i64) -> bool {
    if snp1.scaffold != snp2.scaffold {
        false
    } else {
        (snp1.pos - snp2.pos).abs() < threshold
    }
}

fn is_ko(snp: &Snp, current: &[Snp], threshold: i64) -> bool {
    current.par_iter().any(|old| incompatible(snp, old, threshold))
}

fn main() {
    let args = App::new("Sample")
        .arg(
            Arg::with_name("INPUT_BEAGLE")
                .help("Set the input beagle file to use")
                .required(true),
        )
        .arg(
            Arg::with_name("OUTPUT")
                .short("o")
                .long("out")
                .help("Sets the output file")
                .default_value("out.beagle"),
        )
        .arg(
            Arg::with_name("THRESHOLD")
                .short("t")
                .long("threshold")
                .help("The minimum distance (in bp) between two SNPs to sample")
                .default_value("2000"),
        )
        .get_matches();

    let filename = args.value_of("INPUT_BEAGLE").unwrap();
    let outfile = args.value_of("OUTPUT").unwrap();
    let threshold = value_t!(args, "THRESHOLD", i64).unwrap();

    println!("Reading {}", &filename);
    let mut snps: Vec<Snp> = io::BufReader::new(File::open(filename).unwrap())
        .lines()
        .map(|l| l.expect("Could not parse line"))
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|l| {
            let mut s = l.split('\t');
            let name = s.next().expect(&format!("No name found in {}", l));
            let mut ns = name.split('_');
            Snp {
                scaffold: ns.nth(2).unwrap().to_owned(),
                pos: ns.next().unwrap().parse().unwrap(),
                name: name.to_owned(),
                a_ref: s.next().unwrap().parse().unwrap(),
                a_alt: s.next().unwrap().parse().unwrap(),
                ps: s.map(str::to_owned).collect(),
            }
        })
        .collect();

    println!("Shuffling input");
    snps.shuffle(&mut thread_rng());

    let todo = (0.25 * snps.len() as f32 + 1.) as usize;
    let mut done = vec![snps.pop().unwrap()];
    let bar = ProgressBar::new(snps.len() as u64);
    bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:80} {pos:>8}/{len:8} ({eta} remaining)"),
    );

    println!("Sampling {} SNPs >{}bp apart from {}", todo, threshold, snps.len());
    while let Some(new) = snps.pop() {
        bar.inc(1);
        if !is_ko(&new, &done, threshold) {
            done.push(new);
        }

        if done.len() > todo {
            break;
        }
    }

    if done.len() < todo {
        eprintln!(
            "Not enough valid SNPs; stopping at {} instead of {}",
            done.len(),
            todo
        );
    }
    done.sort_by_key(|snp| (snp.scaffold.to_owned(), snp.pos));

    println!("Writing result to {}", &outfile);
    let mut out = File::create(&outfile).expect(&format!("Can't create `{}`", &outfile));
    for snp in done {
        out.write_all(
            format!(
                "{}\t{}\t{}\t{}\n",
                snp.name,
                snp.a_ref,
                snp.a_alt,
                snp.ps.join("\t"),
            )
            .as_bytes(),
        )
        .unwrap();
    }
}
