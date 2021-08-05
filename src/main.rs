use clap::{App, Arg, SubCommand};
use indicatif::ProgressBar;
use rand::{seq::SliceRandom, thread_rng};
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::io::{self, BufRead};

const THRESHOLD: i64 = 2000;
#[derive(Debug)]
struct Snp {
    scaffold: String,
    pos: i64,

    name: String,
    a_ref: u8,
    a_alt: u8,
    p11: f32,
    p12: f32,
    p13: f32,
    p21: f32,
    p22: f32,
    p23: f32,
}

fn incompatible(snp1: &Snp, snp2: &Snp) -> bool {
    if snp1.scaffold != snp2.scaffold {
        false
    } else {
        (snp1.pos - snp2.pos).abs() < THRESHOLD
    }
}

fn is_ko(snp: &Snp, current: &[Snp]) -> bool {
    current.par_iter().any(|old| incompatible(snp, old))
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
        .get_matches();

    let filename = args.value_of("INPUT_BEAGLE").unwrap();
    let outfile = args.value_of("OUTPUT").unwrap();

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
                p11: s.next().unwrap().parse().unwrap(),
                p12: s.next().unwrap().parse().unwrap(),
                p13: s.next().unwrap().parse().unwrap(),
                p21: s.next().unwrap().parse().unwrap(),
                p22: s.next().unwrap().parse().unwrap(),
                p23: s.next().unwrap().parse().unwrap(),
            }
        })
        .collect();

    println!("Shuffling input");
    snps.shuffle(&mut thread_rng());

    let todo = (0.25 * snps.len() as f32 + 1.) as usize;
    let mut done = vec![snps.pop().unwrap()];
    let bar = ProgressBar::new(snps.len() as u64);
    bar.set_style(indicatif::ProgressStyle::default_bar().template("[{elapsed_precise}] {bar:80} {pos:>8}/{len:8} ({eta} remaining)"));

    println!("Sampling {} SNPs from {}", todo, snps.len());
    while let Some(new) = snps.pop() {
        bar.inc(1);
        if !is_ko(&new, &done) {
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
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                snp.name,
                snp.a_ref,
                snp.a_alt,
                snp.p11,
                snp.p12,
                snp.p13,
                snp.p21,
                snp.p22,
                snp.p23,
            )
                .as_bytes(),
        )
            .unwrap();
    }
}
