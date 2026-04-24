import hashlib
import threading
import time
from typing import Optional
from ftplib import FTP, all_errors
import re
import os
import concurrent.futures
import logging

logger = logging.getLogger(__name__)
MAX_RETRIES = 3


def get_filelist(ftp, pattern):
    files_to_process = []
    re_pattern = re.compile(pattern)

    def filter_ftp_files(line):
        match = re_pattern.search(line)
        if match is not None:
            files_to_process.append(line)

    ftp.retrlines("NLST", filter_ftp_files)
    logger.info(f"found {len(files_to_process)} files matching the required name pattern")
    return files_to_process


def ftp_connect(host, path, config=None):
    if config is None:
        omarc_file = os.getenv("DARWIN_OMA_RC")
        if not omarc_file:
            omarc_file = os.path.abspath(os.path.expanduser("~/.omarc"))
    else:
        omarc_file = os.path.abspath(os.path.expanduser(config))
    if os.path.exists(omarc_file):
        logger.info("reading ftp config from {}".format(omarc_file))
        configs = {}
        with open(omarc_file, "r") as configurations:
            for config_line in configurations:
                if re.search("^(#.*|\n)", config_line):
                    continue
                pair = config_line.strip().split("=")
                configs[pair[0]] = pair[1]
        ftp = FTP(configs["firewall-host"])
        ftp.login(f'anonymous@{configs["firewall-user"]}@{host}', f'pwd@{configs["firewall-password"]}')
    else:
        ftp = FTP(host)
        ftp.login()
    ftp.cwd(path)
    logger.info(f"connected to FTP host {host} on {path}")
    return ftp


def compute_md5_checksum(fname):
    with open(fname, "rb") as f:
        file_hash = hashlib.md5()
        chunk = f.read(8192)
        while chunk:
            file_hash.update(chunk)
            chunk = f.read(8192)
    return file_hash.hexdigest()


_threadLocal = threading.local()


def fetch_file(filename, outdir, host, path, config, crcsums=None):
    if crcsums is None:
        crcsums = {}
    ftp = getattr(_threadLocal, "ftp", None)
    if ftp is None:
        _threadLocal.ftp = ftp = ftp_connect(host, path, config)
    n_retries = 0
    successful_download = False
    while not successful_download:
        try:
            with open("{}/{}".format(outdir, filename), "wb") as fp:
                ftp.retrbinary("RETR {}".format(filename), fp.write)
            if filename in crcsums:
                cksum = compute_md5_checksum(os.path.join(outdir, filename))
                if cksum != crcsums[filename]:
                    logger.error(
                        "wrong md5 checksum for {}: expected {}, computed {}".format(filename, crcsums[filename], cksum)
                    )
                    raise ValueError("wrong checksum for " + filename)
                else:
                    logger.info("checked crc checksum for {}".format(filename))
            successful_download = True
        except all_errors + (ValueError,):
            n_retries += 1
            if n_retries > MAX_RETRIES:
                logger.exception(f"persistent error for {filename}")
                logger.error(f"failed retrieving {filename} {n_retries} times, giving up")
                try:
                    os.remove(os.path.join(outdir, filename))
                except OSError:
                    pass
                return False

            logger.warning(f"failed retrieving {filename}, reconnecting and trying again")
            try:
                ftp.quit()
                delattr(_threadLocal, "ftp")
            except all_errors:
                pass
            time.sleep(2)
            _threadLocal.ftp = ftp_connect(host, path, config)
    return successful_download


def load_crc(host, path, pattern, config):
    ftp = ftp_connect(host, path, config)
    crc_file = get_filelist(ftp, pattern).pop()
    ftp.quit()

    suc = fetch_file(crc_file, "/tmp", host, path, config)
    assert suc
    crc_path = os.path.join("/tmp", crc_file)
    with open(crc_path, "rt") as fh:
        crcs = {}
        for line in fh:
            crc, fname = line.strip().split("\t")
            crcs[fname] = crc
    os.remove(crc_path)
    return crcs


def fetch(
    host: str,
    directory: str,
    pattern: str,
    out_dir: str = "./",
    checksum_ftp_path: Optional[str] = None,
    ftp_config=None,
    nr_cpu: int = 1,
) -> None:
    ftp = ftp_connect(host, directory, ftp_config)
    files_to_process = get_filelist(ftp, pattern)
    if checksum_ftp_path is not None:
        crc_dir, crc_file = os.path.split(checksum_ftp_path)
        crcsums = load_crc(host, crc_dir, crc_file, ftp_config)
    else:
        crcsums = {}
    succ = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=5 * nr_cpu) as executor:
        future_to_fname = {
            executor.submit(fetch_file, fname, out_dir, host, directory, ftp_config, crcsums): fname
            for fname in files_to_process
        }
        for future in concurrent.futures.as_completed(future_to_fname):
            fname = future_to_fname[future]
            try:
                succ += 1 if future.result() else 0
            except Exception as exc:
                logger.warning(f"{fname} generated exception: {exc}")
    logger.info(f"Downloaded {succ} of {len(files_to_process)} files successfully")
