
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

        // TODO: make it a function to test it
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

fn output_line(writer: &mut impl Write,
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


#[cfg(test)]
mod tests_get_reader {
    use super::*;
    use std::io::Read;
    use std::fs;

    #[test]
    fn test_get_reader_with_file() {
        let dir = tempfile::tempdir().unwrap();
        let fpath = dir.path().join("input.txt");
        fs::write(&fpath, "hello\nworld").unwrap();

        let config = Config {
            in_fpath: Some(fpath),
            out_fpath: None,
        };
        let mut reader = get_reader(&config).unwrap();
        let mut content = String::new();
        reader.read_to_string(&mut content).unwrap();
        assert_eq!(content, "hello\nworld");
    }

    #[test]
    fn test_get_reader_nonexistent_file() {
        let config = Config {
            in_fpath: Some(PathBuf::from("/tmp/__nonexistent_ribogrove_test__")),
            out_fpath: None,
        };
        assert!(get_reader(&config).is_err());
    }

    #[test]
    fn test_get_reader_stdin() {
        let config = Config {
            in_fpath: None,
            out_fpath: None,
        };
        assert!(get_reader(&config).is_ok());
    }
}

#[cfg(test)]
mod tests_get_writer {
    use super::*;
    use std::fs;

    #[test]
    fn test_get_writer_with_file() {
        let dir = tempfile::tempdir().unwrap();
        let fpath = dir.path().join("output.txt");

        let config = Config {
            in_fpath: None,
            out_fpath: Some(fpath.clone()),
        };
        {
            let mut writer = get_writer(&config).unwrap();
            writer.write_all(b"test data").unwrap();
        }

        let content = fs::read_to_string(&fpath).unwrap();
        assert_eq!(content, "test data");
    }

    #[test]
    fn test_get_writer_bad_path() {
        let config = Config {
            in_fpath: None,
            out_fpath: Some(PathBuf::from("/nonexistent_dir_12345/out.txt")),
        };
        assert!(get_writer(&config).is_err());
    }

    #[test]
    fn test_get_writer_stdout() {
        let config = Config {
            in_fpath: None,
            out_fpath: None,
        };
        assert!(get_writer(&config).is_ok());
    }
}

#[cfg(test)]
mod tests_check_wanted_dir {
    use super::*;

    #[test]
    fn test_check_wanted_dir_single_match() {
        let wanted = vec!["bacteria", "archaea"];
        assert!(check_wanted_dir("bacteria", &wanted));
    }

    #[test]
    fn test_check_wanted_dir_pipe_contains_match() {
        let wanted = vec!["bacteria", "archaea"];
        assert!(check_wanted_dir("bacteria|complete", &wanted));
    }

    #[test]
    fn test_check_wanted_dir_no_match() {
        let wanted = vec!["bacteria", "archaea"];
        assert!(!check_wanted_dir("fungi", &wanted));
    }

    #[test]
    fn test_check_wanted_dir_pipe_no_match() {
        let wanted = vec!["bacteria", "archaea"];
        assert!(!check_wanted_dir("fungi|complete", &wanted));
    }

    #[test]
    fn test_check_wanted_dir_second_part_match() {
        let wanted = vec!["bacteria", "archaea"];
        assert!(check_wanted_dir("complete|bacteria", &wanted));
    }

    #[test]
    fn test_check_wanted_dir_empty_string() {
        let wanted = vec!["bacteria", "archaea"];
        assert!(!check_wanted_dir("", &wanted));
    }

    #[test]
    fn test_check_wanted_dir_empty_wanted() {
        let wanted: Vec<&str> = vec![];
        assert!(!check_wanted_dir("bacteria", &wanted));
    }
}

#[cfg(test)]
mod tests_output_line {
    use super::*;
    use std::fs;

    #[test]
    fn test_output_line_writes_string_and_newline() {
        let dir = tempfile::tempdir().unwrap();
        let fpath = dir.path().join("out.txt");
        let file = fs::File::create(&fpath).unwrap();
        let mut writer: Box<dyn Write> = Box::new(std::io::BufWriter::new(file));

        output_line(&mut writer, "hello").unwrap();
        drop(writer);

        let content = fs::read_to_string(&fpath).unwrap();
        assert_eq!(content, "hello\n");
    }

    #[test]
    fn test_output_line_writes_newline_only_for_empty_string() {
        let dir = tempfile::tempdir().unwrap();
        let fpath = dir.path().join("out.txt");
        let file = fs::File::create(&fpath).unwrap();
        let mut writer: Box<dyn Write> = Box::new(std::io::BufWriter::new(file));

        output_line(&mut writer, "").unwrap();
        drop(writer);

        let content = fs::read_to_string(&fpath).unwrap();
        assert_eq!(content, "\n");
    }
}
