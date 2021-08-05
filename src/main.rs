use rand::{thread_rng, seq::SliceRandom};
use std::{env, io::Write};
use std::fs::File;
use std::io::{self, BufRead};
use rayon::prelude::*;

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
    let mut args = env::args();
    let _ = args.next();
    let filename = args.next().expect("Please specify input file name");
    let outfile = args.next().expect("Please specify output file name");

    println!("Reading {}; writing result to {}", &filename, &outfile);
    let mut snps: Vec<Snp> = io::BufReader::new(File::open(filename).unwrap())
        .lines()
        .map(|l| {
            let l = l.unwrap();
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
    snps.shuffle(&mut thread_rng());

    let todo = (0.25 * snps.len() as f32 + 1.) as usize;
    let mut done = vec![snps.pop().unwrap()];

    while let Some(new) = snps.pop() {
        if !is_ko(&new, &done) {
            done.push(new);
        }

        if done.len() > todo {
            break;
        }
    }

    if done.len() < todo {
        eprintln!("Not enough valid SNPs; stopping at {} instead of {}", done.len(), todo);
    }
    done.sort_by_key(|snp| (snp.scaffold.to_owned(), snp.pos));

    let mut out = File::create(&outfile).expect(&format!("Can't create `{}`", &outfile));
    for snp in done {
        out.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                              snp.name, snp.a_ref, snp.a_alt,
                              snp.p11, snp.p12, snp.p13,
                              snp.p21, snp.p22, snp.p23,
        ).as_bytes()).unwrap();
    }
}
