
import os
import sys
import glob
import gzip
import time
import subprocess as sp

from Bio import SeqIO
from pandas import Series

from src.file_navigation import get_asm_data_dir_path, \
                                get_asm_report_fpath, \
                                get_genome_seqannot_fpath


class GenomeDownloader:

    N_ATTEMPTS      = 3
    FAIL_SLEEP_TIME = 3

    def __init__(self, asm_sum_row: Series, outdir):
        self.assembly_accession = str(asm_sum_row['asm_acc'])
        self.asm_basedir_url    = str(asm_sum_row['ftp_path'])
        self.asm_name           = self._parse_asm_name(
            self.asm_basedir_url, self.assembly_accession
        )
        self.genome_dirpath = get_asm_data_dir_path(self.assembly_accession, outdir)
        self.asm_report_fpath = get_asm_report_fpath(
            self.assembly_accession,
            outdir
        )
        self.asm_report_check_time_fpath = self.asm_report_fpath + '.modify_timestamp'
        self.asm_seqannot_fpath = get_genome_seqannot_fpath(
            self.assembly_accession,
            outdir
        )
        self.asm_seqannot_check_time_fpath = self.asm_seqannot_fpath + '.modify_timestamp'
    # end def


    def _parse_asm_name(self, ftp_path, asm_acc):
        # Dir url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/012/635/GCF_900012635.1_Pyrococcus_chitonophagus_genome_sequence
        #  We need this:                                                                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # We get asm_name from ftp_path because asm_sum_row['asm_name'] is sometimes inconsistent with actual URLs.
        #  Some examples of inconsistency:
        #    GCF_902166805.1, GCF_000577895.1, GCF_900095285.1, GCF_905187425.1
        return ftp_path \
            .split('/')[-1] \
            .replace('{}_'.format(asm_acc), '')
    # end def


    def try_download(self):
        already_here = self._check_if_genome_is_alredy_here()
        if already_here:
            return DownloadStatus(DownloadStatus.ALREADY_HERE)
        # end if

        success, err_msg = False, ''
        for _ in range(self.N_ATTEMPTS):
            try:
                self._download(already_here)
                self._check_downloaded_files()
            except DownloadError as err:
                err_msg = str(err)
                time.sleep(self.FAIL_SLEEP_TIME)
            else:
                success = True
                break
            # end try
        # end while

        if success:
            return DownloadStatus(DownloadStatus.DOWNLOADED)
        else:
            return DownloadStatus(DownloadStatus.FAILED, err_msg)
        # end if
    # end def

    def _download(self, already_here):
        if not already_here and self._genome_dir_exists():
            self._empty_genome_dir()
        # end if
        self._download_genome_files()
    # end def

    def _check_if_genome_is_alredy_here(self):
        try:
            self._check_genome_files_exist()
            self._check_asm_report_content()
            self._check_seqannot_content()
        except DownloadError:
            return False
        else:
            return True
        # end try
    # end def

    def _check_downloaded_files(self):
        self._check_genome_files_exist()
        self._check_asm_report_content()
        self._check_seqannot_content()
    # end def

    def _check_genome_files_exist(self):
        if not self._genome_dir_exists():
            raise FilesMissingError('Genome dir does not exist')
        # end if

        asm_report_is_ok = self._check_asm_report_exists()
        seqannot_is_ok = self._check_seqannot_file_exists()
        if not (asm_report_is_ok and seqannot_is_ok):
            raise FilesMissingError('Cannot find genome files')
        # end if
    # end def

    def _genome_dir_exists(self):
        return os.path.isdir(self.genome_dirpath)
    # end def

    def _check_asm_report_exists(self):
        exists = os.path.isfile(self.asm_report_fpath)
        if not exists:
            return False
        # end if
        not_empty = os.path.getsize(self.asm_report_fpath) > 0
        return not_empty
    # end def

    def _check_seqannot_file_exists(self):
        exists = os.path.isfile(self.asm_seqannot_fpath)
        if not exists:
            return False
        # end if
        not_empty = os.path.getsize(self.asm_seqannot_fpath) > 0
        return not_empty
    # end def

    def _check_asm_report_content(self):
        if self._asm_report_modified():
            with open(self.asm_report_fpath, 'rt') as infile:
                if 'Error 404' in infile.read():
                    raise Error404('Assembly report file not found at NCBI server: 404.')
                # end if
            # end with
            self._rewrite_asm_report_timestamp()
        # end if
    # end def

    def _asm_report_modified(self):
        if not os.path.exists(self.asm_report_check_time_fpath):
            return True
        # end if
        return os.path.getmtime(self.asm_report_check_time_fpath) \
               < os.path.getmtime(self.asm_report_fpath)
    # end def

    def _rewrite_asm_report_timestamp(self):
        if os.path.exists(self.asm_report_check_time_fpath):
            os.unlink(self.asm_report_check_time_fpath)
        # end if
        with open(self.asm_report_check_time_fpath, 'wt') as file:
            pass
        # end with
    # end def

    def _check_seqannot_content(self):
        if self._asm_seqannot_modified():
            try:
                self._test_gb_gz_file(self.asm_seqannot_fpath)
            except gzip.BadGzipFile as err:
                if self._seqannot_file_is_404():
                    raise Error404('.gbff.gz file not found at NCBI server: 404.')
                else:
                    raise gzip.BadGzipFile('.gbff.gz file is not a gzip file')
                # end if
            except NoSeqsError as err:
                raise NoSeqsError('.gbff.gz file is not a GenBank file.')
            # end try
            self._rewrite_asm_seqannot_timestamp()
        # end if
    # end def

    def _asm_seqannot_modified(self):
        if not os.path.exists(self.asm_seqannot_check_time_fpath):
            return True
        # end if
        return os.path.getmtime(self.asm_seqannot_check_time_fpath) \
               < os.path.getmtime(self.asm_seqannot_fpath)
    # end def

    def _rewrite_asm_seqannot_timestamp(self):
        if os.path.exists(self.asm_seqannot_check_time_fpath):
            os.unlink(self.asm_seqannot_check_time_fpath)
        # end if
        with open(self.asm_seqannot_check_time_fpath, 'wt') as file:
            pass
        # end with
    # end def

    def _test_gb_gz_file(self, fpath):
        with gzip.open(fpath, 'rt') as infile:
            seq_records = tuple(SeqIO.parse(infile, 'gb'))
            if len(seq_records) == 0:
                raise NoSeqsError()
            # end if
        # end with
    # end def

    def _seqannot_file_is_404(self):
        try:
            with gzip.open(self.asm_seqannot_fpath, 'rt') as infile:
                if 'Error 404' in infile.read():
                    return True
                # end if
            # end with
        except gzip.BadGzipFile:
            return True
        # end try
        return False
    # end def

    def _empty_genome_dir(self):
        files_to_rm = glob.iglob(
            os.path.join(self.genome_dirpath, '*')
        )
        for fpath in files_to_rm:
            try:
                os.unlink(fpath)
            except OSError as err:
                err_msg = 'Error: cannot remove file `{}`: {}' \
                    .format(fpath, err)
                print(err_msg)
                raise DownloadError(err_msg)
        # end def
    # end def

    def _download_genome_files(self):
        if not self._genome_dir_exists():
            self._make_genome_dir()
        # end if
        self._download_asm_report()
        self._download_seqannot()
    # end def

    def _make_genome_dir(self):
        try:
            os.makedirs(self.genome_dirpath)
        except OSError as err:
            err_msg = 'Error: cannot create directory `{}`: {}' \
                .format(genome_dirpath, err)
            print(err_msg)
            raise DownloadError(err_msg)
        # end try
    # end def

    def _download_asm_report(self):
        download_url = self._make_asm_report_url()
        save_fpath = self.asm_report_fpath
        self._download_file(download_url, save_fpath)
    # end def

    def _make_asm_report_url(self):
        # E.g.
        # Dir url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/762/265/GCF_000762265.1_ASM76226v1
        # File name: GCF_000762265.1_ASM76226v1_assembly_report.txt
        # Or
        # Dir url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/012/635/GCF_900012635.1_Pyrococcus_chitonophagus_genome_sequence
        # File name: GCF_900012635.1_Pyrococcus_chitonophagus_genome_sequence_assembly_report.txt
        return '{}/{}_{}_assembly_report.txt'.format(
            self.asm_basedir_url,
            self.assembly_accession,
            self.asm_name
        )
    # end def

    def _download_file(self, download_url, save_fpath):
        command = ' '.join(
            [
                'curl',
                '"{}"'.format(download_url),
                '-o "{}"'.format(save_fpath),
            ]
        )
        pipe = sp.Popen(
            command,
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            encoding='utf-8'
        )

        out, err = pipe.communicate()
        if pipe.returncode != 0:
            print(err)
            raise DownloadError(err)
        # end if
    # end def

    def _download_seqannot(self):
        download_url = self._make_seqannot_url()
        save_fpath = self.asm_seqannot_fpath
        self._download_file(download_url, save_fpath)
    # end def

    def _make_seqannot_url(self):
        # E.g.
        # Dir url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/762/265/GCF_000762265.1_ASM76226v1
        # File name: GCF_000762265.1_ASM76226v1_genomic.gbff.gz
        # Or
        # Dir url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/012/635/GCF_900012635.1_Pyrococcus_chitonophagus_genome_sequence
        # File name: GCF_900012635.1_Pyrococcus_chitonophagus_genome_sequence_genomic.gbff.gz
        return '{}/{}_{}_genomic.gbff.gz'.format(
            self.asm_basedir_url,
            self.assembly_accession,
            self.asm_name
        )
    # end def
# end class


class DownloadStatus:

    ALREADY_HERE = 0
    DOWNLOADED   = 1
    FAILED       = 2

    def __init__(self, status_code, err_msg=None):
        self.status_code = status_code
        self.err_msg     = err_msg
    # end def
# end class


class DownloadError(Exception):
    pass
# end class

class FilesMissingError(DownloadError):
    pass
# end class

class NoSeqsError(DownloadError):
    pass
# end class

class Error404(DownloadError):
    pass
# end class
