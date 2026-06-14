
use std::env;
use std::fs::File;
use std::error::Error;
use std::path::PathBuf;
use std::process::ExitCode;
use std::collections::HashSet;
use std::io::{self, BufRead, BufReader, Write, BufWriter};


struct Config {
    in_fpath: Option<PathBuf>,
    out_fpath: Option<PathBuf>,
}

impl Config {
    fn new() -> Config {
        let mut args = env::args();
        args.next(); // pass slement 0

        let in_fpath = match args.next() {
            Some(arg) => {
                PathBuf::from(arg)
            },
            None => {
                return Config{
                    in_fpath: None,
                    out_fpath: None,
                };
            },
        };

        let out_fpath = match args.next() {
            Some(arg) => {
                PathBuf::from(arg)
            },
            None => {
                return Config {
                    in_fpath: Some(in_fpath), // TODO: check existence, readability
                    out_fpath: None,
                };
            },
        };

        Config {
            in_fpath: Some(in_fpath),
            out_fpath: Some(out_fpath), // TODO: check dir existence
        }
    }

    fn print_self(&self) {
        eprintln!("Config:");
        eprint!("  in_fpath: ");
        match &self.in_fpath {
            Some(fpath) => {
                eprintln!("`{}`", fpath.display());
            },
            None => {
                eprintln!("None: reading from stdin");
            }
        }
        eprint!("  out_fpath: ");
        match &self.out_fpath {
            Some(fpath) => {
                eprintln!("`{}`", fpath.display());
            },
            None => {
                eprintln!("None: writing to stdout");
            }
        }
        eprintln!("");
    }
}


fn main() -> ExitCode {
    let config = Config::new();
    config.print_self();

    // TODO: make outdir

    let exit_code = match read_filter_and_write(&config) {
        Ok(..) => ExitCode::SUCCESS,
        Err(..) => ExitCode::FAILURE,
    };

    exit_code
}



fn read_filter_and_write(config: &Config) -> Result<(), ()> {
    // See input data structure: https://ftp.ncbi.nlm.nih.gov/refseq/release/README
    // Example input line:
    //    "7\tAzorhizobium caulinodans\tNZ_JBAFXE010000103.1\tbacteria|complete\tna\t18178"

    const ACC_COL_IDX: usize = 2;
    const DIR_COL_IDX: usize = 3;

    let reader: Box<dyn BufRead> = get_reader(&config)?;
    let mut writer: Box<dyn Write> = get_writer(&config)?;

    let unwanted_prefixes = make_unwanted_prefix_set();
    let wanted_dirs = make_wanted_dir_vec();

    for line in reader.lines() {
        if let Err(err) = line {
            eprintln!("Error reading input line: {err}");
            return Err(());
        }
        let line_str = line.unwrap();

        let col_values: Vec<&str> = line_str.split('\t').collect();

        let acc_prefix = &col_values[ACC_COL_IDX][0..3];
        if unwanted_prefixes.contains(acc_prefix) {
            continue;
        }

        let dir_str = &col_values[DIR_COL_IDX];
        if !check_wanted_dir(dir_str, &wanted_dirs) {
            continue;
        }

        if let Err(err) = output_line(&mut writer, &line_str) {
            eprintln!("Error writing line `{:?}`", line_str);
            eprintln!("{err}");
            return Err(());
        };
    }

    Ok(())
}

fn get_reader(config: &Config) -> Result<Box<dyn BufRead>, ()> {
    match &config.in_fpath {
        Some(fpath) => {
            match File::open(fpath) {
                Ok(handle) => {
                    return Ok(Box::new(BufReader::new(handle)));
                },
                Err(err) => {
                    eprintln!("Cannot open input file: {err}");
                    return Err(());
                }
            }
        },
        None => Ok(
            Box::new(io::stdin().lock())
        ),
    }
}

fn get_writer(config: &Config) -> Result<Box<dyn Write>, ()> {
    match &config.out_fpath {
        Some(fpath) => {
            match File::create(fpath) {
                Ok(handle) => {
                    return Ok(Box::new(BufWriter::new(handle)));
                },
                Err(err) => {
                    eprintln!("Cannot open output file: {err}");
                    return Err(());
                }
            }
        },
        None => Ok(
            Box::new(io::stdout().lock())
        )
    }
}

fn make_unwanted_prefix_set() -> HashSet<&'static str> {
    // Unwanted means non-genomic
    // https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
    HashSet::from([
        "NM_", "NR_", "XM_", "XR_", "AP_", "NP_", "YP_", "XP_", "WP_",
    ])
}

fn make_wanted_dir_vec() -> Vec<&'static str> {
    vec![
        "bacteria",
        "archaea",
    ]
}

fn check_wanted_dir(dir_str: &str,
                    wanted_dirs: &Vec<&str>) -> bool {
    for dir in dir_str.split('|') {
        if wanted_dirs.contains(&dir) {
            return true;
        }
    }
    false
}

fn output_line(writer: &mut Box<dyn Write>,
               out_string: &str) -> Result<(), Box<dyn Error>> {
    let out_strings = vec![
        out_string,
        "\n",
    ];
    for s in out_strings {
        writer.write(&s.as_bytes())?;
    }
    Ok(())
}
