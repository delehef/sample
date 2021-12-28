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

fn main() -> Result<()> {
    let args = App::new("Sample")
        .arg(
            Arg::with_name("INPUT_BEAGLE")
                .long("--in-beagle")
                .takes_value(true)
                .help("Set the input file to use in the beagle format")
        )
        .arg(
            Arg::with_name("INPUT_SIMPLE")
                .long("--in-simple")
                .takes_value(true)
                .help("Set the input file to use in the scaffold:position format")
        )
        .group(ArgGroup::with_name("INPUT")
               .args(&["INPUT_BEAGLE", "INPUT_SIMPLE"])
               .required(true))
        .arg(
            Arg::with_name("OUTPUT")
                .short("o")
                .long("out")
                .help("Sets the output file")
                .takes_value(true)
                .default_value("out.snps"),
        )
        .arg(
            Arg::with_name("OUTPUT_FORMAT")
                .long("format")
                .takes_value(true)
                .possible_values(&["beagle", "simple"])
                .default_value("simple")
        )
        .arg(
            Arg::with_name("THRESHOLD")
                .short("t")
                .long("threshold")
                .help("The minimum distance (in bp) between two SNPs to sample")
                .takes_value(true)
                .default_value("2000"),
        )
        .arg(
            Arg::with_name("QUANTITY")
                .short("Q")
                .long("quantity")
                .help("The maximum proportion of SNPs to sample. 1.0 means sample all, 0.5 means sample half.")
                .takes_value(true)
                .default_value("1")
        )
        .get_matches();

    let mut snps: Vec<Snp> = if args.is_present("INPUT_BEAGLE") {
        let filename = value_t!(args, "INPUT_BEAGLE", String)?;
        println!("Reading {}", &filename);
        io::BufReader::new(File::open(filename).unwrap())
            .lines()
            .map(|l| l.expect("Could not read line"))
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
            .collect()
    } else if args.is_present("INPUT_SIMPLE"){
        let filename = value_t!(args, "INPUT_SIMPLE", String)?;
        println!("Reading {}", &filename);
        io::BufReader::new(File::open(&filename).unwrap())
            .lines()
            .map(|l| l.expect("Could not read line"))
            .filter(|l| !l.is_empty())
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|l| {
                let mut s = l.split('\t');
                Snp {
                    scaffold: s.next().expect(&format!("No scaffold found in {}", l)).to_owned(),
                    pos: s.next().expect(&format!("No position found in {}", l)).parse().unwrap(),
                    name: String::new(),
                    a_ref: 0,
                    a_alt: 0,
                    ps: Vec::new(),
                }
            })
            .collect()
    } else {
        panic!("Please set an input file with --in-beagle or --in-simple")
    };
    let outfile = value_t!(args, "OUTPUT", String)?;
    let threshold = value_t!(args, "THRESHOLD", i64)?;


    println!("Shuffling input");
    snps.shuffle(&mut thread_rng());

    let bar = ProgressBar::new(snps.len() as u64);
    let mut done = vec![snps.pop().unwrap()];
    bar.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:80} {pos:>8}/{len:8} ({eta} remaining)"),
    );


    let quantity = value_t!(args, "QUANTITY", f32)?;
    let todo = (quantity * snps.len() as f32 + 1.) as usize;
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
    let format = value_t!(args, "OUTPUT_FORMAT", String)?;
    for snp in done {
        out.write_all(format_snp_out(&snp, &format).as_bytes())
        .unwrap();
    }
    Ok(())
}

fn format_snp_out(snp: &Snp, format: &str) -> String {
    match format {
        "simple" => format!("{}:{}\n", &snp.scaffold, snp.pos),
        "beagle" => format!("{}\t{}\t{}\t{}\n", snp.name, snp.a_ref, snp.a_alt, snp.ps.join("\t")),
        _ => unimplemented!(),
    }
}
