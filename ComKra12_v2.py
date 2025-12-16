#!/usr/bin/env python3
"""
kraken_kmer_compare.py

Purpose
-------
Compare per-read k-mer annotation results between two Kraken output files
and evaluate whether sequences unclassified in the 2nd run map to references
in the provided df1 reference file.

Input files (command-line)
- out1: Kraken output file from first annotation (per-read output)
- out2: Kraken output file from second annotation (per-read output)
- df1 : two-column tab file mapping taxid -> species string (e.g. "3120256\tGilliamella sp. wkB108")

Output
------
Tab-delimited file with columns:
  title	out1_kmer_info	out2_kmer_info	ref_taxid	ref_species	composite_kmer_counts	final_judgement

Notes
-----
- The script expects Kraken "--output" style files. By default it treats
  columns as: 1=status (C/U), 2=title, 3=taxid, 5=kmer-annotation-field.
  If your files use different column layout, you may change indices via
  the command-line options --title-col/--taxid-col/--kmer-col (1-based).
- Final judgement rules (default multiplier = 3):
  * If the top taxid in the composite kmer counts equals ref_taxid -> T
  * Else if taxid '0' is top, ref_taxid is second, and count(0) >= multiplier * count(ref_taxid) -> F
  * Else if ref_taxid ranks second and count(0) < multiplier * count(ref_taxid) -> T
  * Else -> F

Usage example
-------------
python3 kraken_kmer_compare.py --out1 out1.kraken --out2 out2.kraken --df1 df1.txt --out-prefix result

Notes on this modified version
------------------------------
- Added --threads to parallelize processing of unannotated titles.
- Use a SQLite index for out1 to avoid loading the whole file into memory.
- Verbose logging (-v/--verbose) is ON by default; use --no-verbose to turn it off.
"""

import argparse
import csv
import re
import sys
from collections import Counter
import logging
import concurrent.futures
import sqlite3
import threading
import os

# Thread-local storage for per-thread sqlite connections
_thread_local = threading.local()


def parse_args():
    p = argparse.ArgumentParser(description='Compare k-mer annotation between two Kraken outputs and evaluate targets using df1.')
    p.add_argument('--out1', required=True, help='Kraken output file from first annotation')
    p.add_argument('--out2', required=True, help='Kraken output file from second annotation')
    p.add_argument('--df1', required=True, help='Two-column tab file: taxid\tspecies_name')
    p.add_argument('--out-prefix', required=True, help='Output file prefix')
    p.add_argument('--multiplier', type=float, default=3.0, help='Multiplier used in judgement rule (default: 3)')
    p.add_argument('--title-col', type=int, default=2, help='1-based column index for sequence title in kraken outputs (default: 2)')
    p.add_argument('--taxid-col', type=int, default=3, help='1-based column index for taxid in kraken outputs (default: 3)')
    p.add_argument('--kmer-col', type=int, default=5, help='1-based column index for per-read kmer annotation field (default: 5)')
    p.add_argument('--min-k', type=int, default=0, help='Minimum kmer count to consider (not enforced automatically, provided for future use)')
    # -v default ON per user request; add --no-verbose to allow disabling
    p.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=True, help='Verbose logging (default ON)')
    p.add_argument('--no-verbose', dest='verbose', action='store_false', help='Disable verbose logging')
    p.add_argument('--threads', type=int, default=1, help='Number of worker threads to use for processing unannotated titles (default: 1)')
    # Optionally allow forcing index rebuild (handy)
    p.add_argument('--rebuild-index', action='store_true', help='Force rebuild the out1 sqlite index even if index file exists')
    return p.parse_args()


def setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='[%(levelname)s] %(message)s')


def load_df1(path):
    """Load df1 mapping file: taxid -> species string.
    Return dict taxid(str) -> species(str-without-spaces-replaced-by-hyphen)
    """
    d = {}
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                # tolerate space-separated
                parts = line.split()
                if len(parts) < 2:
                    logging.warning('Skipping bad df1 line: %s', line)
                    continue
            taxid = parts[0].strip()
            species = '\t'.join(parts[1:]).strip()
            # remove internal whitespace by replacing with hyphen
            species_clean = re.sub(r'\s+', '_', species)
            d[taxid] = species_clean
    logging.info('Loaded %d df1 entries from %s', len(d), path)
    return d


def parse_kmer_field(field_text: str):
    """Parse a kmer-annotation text and return a Counter{taxid: count}.

    Expected common formats:
      - "taxid:count taxid:count ..."
      - "taxid:count;taxid:count"
      - other tokens containing taxid:count pairs
    Fallback: try to extract any digits:num pairs with regex.
    """
    if not field_text:
        return Counter()
    # try to find all id:count pairs
    pairs = re.findall(r"(\d+):(\d+)", field_text)
    if pairs:
        cnt = Counter()
        for tid, c in pairs:
            cnt[tid] += int(c)
        return cnt
    # fallback: split on whitespace or semicolon and look for token like tid:count
    tokens = re.split(r'[\s;|,]+', field_text)
    cnt = Counter()
    for tok in tokens:
        if ':' in tok:
            a, b = tok.split(':', 1)
            if a.isdigit() and b.isdigit():
                cnt[a] += int(b)
    return cnt


def parse_kraken_out(path, title_col=2, taxid_col=3, kmer_col=5):
    """Parse kraken per-read output into dicts keyed by title.

    Returns:
      records: dict title -> dict(status, taxid, kmer_field, kmer_counter, raw_line)
    Note: this function is still used to parse out2 into memory (out2 is used to find
    unannotated titles and quick checks). For out1 we will instead build/use a sqlite index.
    """
    records = {}
    with open(path, 'r', encoding='utf-8') as fh:
        for ln, line in enumerate(fh, start=1):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            # safe access with 1-based indices provided
            status = parts[0] if len(parts) >= 1 else ''
            title = parts[title_col - 1] if len(parts) >= title_col else ''
            taxid = parts[taxid_col - 1] if len(parts) >= taxid_col else '0'
            kmer_field = parts[kmer_col - 1] if len(parts) >= kmer_col else ''
            kmer_counter = parse_kmer_field(kmer_field)
            records[title] = {
                'status': status,
                'taxid': taxid,
                'kmer_field': kmer_field,
                'kmer_counter': kmer_counter,
                'raw_line': line,
                'line_no': ln
            }
    logging.info('Parsed %d records from %s', len(records), path)
    return records


def composite_counts(c1: Counter, c2: Counter):
    """Return composite Counter adding counts from two Counters."""
    c = Counter()
    c.update(c1)
    c.update(c2)
    return c


def sorted_composite_str(counter_obj: Counter):
    """Return string like taxid:count>taxid:count sorted desc by count."""
    if not counter_obj:
        return ''
    items = sorted(counter_obj.items(), key=lambda x: (-x[1], x[0]))
    return '>'.join(f"{tid}:{cnt}" for tid, cnt in items)


def determine_ref_taxid_from_out1(kmer_counter: Counter, df1_taxids_set):
    """From out1 kmer counter pick the taxid that exists in df1 (highest count among those).
    Return taxid or 'NA'."""
    if not kmer_counter:
        return 'NA'
    # filter keys that are in df1
    candidates = [(tid, cnt) for tid, cnt in kmer_counter.items() if tid in df1_taxids_set]
    if not candidates:
        return 'NA'
    # pick highest count
    candidates.sort(key=lambda x: -x[1])
    return candidates[0][0]


def final_judgement(composite_counter: Counter, ref_taxid: str, multiplier: float):
    """Apply the judgement rules described in the specification.

    Returns 'T' or 'F'.
    """
    if ref_taxid == 'NA' or not composite_counter:
        return 'F'
    items = sorted(composite_counter.items(), key=lambda x: (-x[1], x[0]))
    top_tid, top_cnt = items[0]
    # if top equals ref_taxid -> T
    if top_tid == ref_taxid:
        return 'T'
    # else if top is '0' and ref_taxid is second
    if top_tid == '0':
        # find rank of ref_taxid
        ref_rank = None
        ref_cnt = 0
        for idx, (tid, cnt) in enumerate(items):
            if tid == ref_taxid:
                ref_rank = idx + 1
                ref_cnt = cnt
                break
        if ref_rank is None:
            return 'F'
        # ref is second and top_cnt >= multiplier * ref_cnt => F
        if ref_rank == 2 and top_cnt >= multiplier * (ref_cnt if ref_cnt>0 else 1):
            return 'F'
        # ref is second and top_cnt < multiplier * ref_cnt => T
        if ref_rank == 2 and top_cnt < multiplier * ref_cnt:
            return 'T'
        # otherwise F
        return 'F'
    # other cases: if ref_taxid rank is 2 and top not 0, or ref rank 3+ => F
    # if ref is second but top not 0 -> F per spec
    return 'F'


######## SQLite index helpers for out1 (streaming approach) ########

def build_out1_sqlite_index(out1_path, db_path, title_col=2, taxid_col=3, kmer_col=5):
    """
    Build a sqlite index file for out1 so we don't need to hold all out1 records in memory.
    Table schema: out1(title PRIMARY KEY, status, taxid, kmer_field, raw_line, line_no)
    If db_path already exists, this function does nothing (to reuse existing index).
    """
    if os.path.exists(db_path):
        logging.info('SQLite index already exists at %s, skipping build.', db_path)
        return

    logging.info('Building sqlite index for out1 (%s) -> %s', out1_path, db_path)
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("PRAGMA synchronous = OFF;")
    cur.execute("PRAGMA journal_mode = MEMORY;")
    cur.execute("""
        CREATE TABLE IF NOT EXISTS out1 (
            title TEXT PRIMARY KEY,
            status TEXT,
            taxid TEXT,
            kmer_field TEXT,
            raw_line TEXT,
            line_no INTEGER
        );
    """)
    conn.commit()

    insert_sql = "INSERT OR REPLACE INTO out1 (title, status, taxid, kmer_field, raw_line, line_no) VALUES (?, ?, ?, ?, ?, ?)"
    batch = []
    batch_size = 2000
    with open(out1_path, 'r', encoding='utf-8') as fh:
        for ln, line in enumerate(fh, start=1):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            status = parts[0] if len(parts) >= 1 else ''
            title = parts[title_col - 1] if len(parts) >= title_col else ''
            taxid = parts[taxid_col - 1] if len(parts) >= taxid_col else '0'
            kmer_field = parts[kmer_col - 1] if len(parts) >= kmer_col else ''
            batch.append((title, status, taxid, kmer_field, line, ln))
            if len(batch) >= batch_size:
                cur.executemany(insert_sql, batch)
                conn.commit()
                batch = []
    if batch:
        cur.executemany(insert_sql, batch)
        conn.commit()
    conn.close()
    logging.info('Finished building sqlite index.')


def get_thread_db_conn(db_path):
    """
    Return a per-thread sqlite3 connection (stored in thread-local).
    Use check_same_thread=False to be safe; each thread gets its own connection.
    """
    conn = getattr(_thread_local, 'conn', None)
    if conn is None:
        conn = sqlite3.connect(db_path, check_same_thread=False)
        # return rows as tuples (we parse fields ourselves)
        _thread_local.conn = conn
    return conn


def fetch_out1_record_from_db(db_path, title):
    """
    Fetch a single record for title from the sqlite index.
    Returns a dict like parse_kraken_out produces, or None if not found.
    """
    conn = get_thread_db_conn(db_path)
    cur = conn.cursor()
    cur.execute("SELECT status, taxid, kmer_field, raw_line, line_no FROM out1 WHERE title = ?", (title,))
    row = cur.fetchone()
    if not row:
        return None
    status, taxid, kmer_field, raw_line, line_no = row
    return {
        'status': status,
        'taxid': taxid,
        'kmer_field': kmer_field,
        'kmer_counter': parse_kmer_field(kmer_field),
        'raw_line': raw_line,
        'line_no': line_no
    }


########################## Main ##########################

def main():
    args = parse_args()
    setup_logging(args.verbose)

    if args.threads is None or args.threads < 1:
        args.threads = 1

    df1 = load_df1(args.df1)
    df1_taxids = set(df1.keys())

    # parse out2 first and count annotated sequences existence in df1
    out2_records = parse_kraken_out(args.out2, title_col=args.title_col, taxid_col=args.taxid_col, kmer_col=args.kmer_col)

    total_checked = 0
    exist_in_df1 = 0
    not_exist_in_df1 = 0
    unannotated_titles = []

    for title, rec in out2_records.items():
        status = rec['status']
        taxid = rec['taxid'] if rec['taxid'] else '0'
        # treat 'C' as classified; also treat taxid != '0' as annotated
        if (status and status.startswith('C')) or (taxid != '0' and taxid != 'NA'):
            # this sequence participates in df1 check
            total_checked += 1
            if taxid in df1_taxids:
                exist_in_df1 += 1
            else:
                not_exist_in_df1 += 1
        else:
            # unannotated in out2, add to review set
            unannotated_titles.append(title)

    logging.info('Out2: total sequences participating in df1 check: %d', total_checked)
    logging.info('Out2: sequences with taxid present in df1: %d', exist_in_df1)
    logging.info('Out2: sequences with taxid NOT present in df1: %d', not_exist_in_df1)
    logging.info('Out2: unannotated sequences to review: %d', len(unannotated_titles))

    outpath = f"{args.out_prefix}.tsv"
    # write header early
    with open(outpath, 'w', encoding='utf-8') as outfh:
        writer = csv.writer(outfh, delimiter='\t', lineterminator='\n')
        header = ['title','out1_kmer_info','out2_kmer_info','ref_taxid','ref_species','composite_kmer_counts','final_judgement']
        writer.writerow(header)

    if not unannotated_titles:
        logging.info('No unannotated sequences in out2. Exiting.')
        logging.info('Wrote empty output %s', outpath)
        return

    # Prepare sqlite index for out1 (streaming)
    db_path = f"{args.out_prefix}.out1.sqlite"
    if args.rebuild_index and os.path.exists(db_path):
        try:
            os.remove(db_path)
            logging.info('Removed existing index due to --rebuild-index flag.')
        except Exception as e:
            logging.warning('Could not remove existing sqlite index: %s', e)

    build_out1_sqlite_index(args.out1, db_path, title_col=args.title_col, taxid_col=args.taxid_col, kmer_col=args.kmer_col)

    # Worker function to process one title and return a row list
    def worker(title):
        # fetch out1 record from sqlite index (streaming approach)
        rec1 = fetch_out1_record_from_db(db_path, title)
        rec2 = out2_records.get(title, None)
        out1_field = rec1['kmer_field'] if rec1 else ''
        out2_field = rec2['kmer_field'] if rec2 else ''
        out1_counter = rec1['kmer_counter'] if rec1 else Counter()
        out2_counter = rec2['kmer_counter'] if rec2 else Counter()
        # determine ref_taxid from out1 kmer results
        ref_taxid = determine_ref_taxid_from_out1(out1_counter, df1_taxids)
        ref_species = df1.get(ref_taxid, 'NA') if ref_taxid != 'NA' else 'NA'
        # composite counts (sum of out1 and out2 for this title)
        comp = composite_counts(out1_counter, out2_counter)
        comp_str = sorted_composite_str(comp)
        judgement = final_judgement(comp, ref_taxid, args.multiplier)
        return [title, out1_field, out2_field, ref_taxid, ref_species, comp_str, judgement]

    processed = 0
    # If threads == 1, do simple loop (slightly faster/no overhead)
    if args.threads == 1:
        logging.info('Processing %d unannotated titles sequentially...', len(unannotated_titles))
        with open(outpath, 'a', encoding='utf-8') as outfh:
            writer = csv.writer(outfh, delimiter='\t', lineterminator='\n')
            for title in unannotated_titles:
                row = worker(title)
                writer.writerow(row)
                processed += 1
                if processed % 1000 == 0:
                    logging.info('Processed %d / %d unannotated sequences', processed, len(unannotated_titles))
    else:
        max_workers = args.threads
        logging.info('Processing %d unannotated titles using %d threads...', len(unannotated_titles), max_workers)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            # executor.map preserves input order; iterate results and write sequentially
            with open(outpath, 'a', encoding='utf-8') as outfh:
                writer = csv.writer(outfh, delimiter='\t', lineterminator='\n')
                for i, row in enumerate(ex.map(worker, unannotated_titles), start=1):
                    writer.writerow(row)
                    if i % 1000 == 0:
                        logging.info('Processed %d / %d unannotated sequences', i, len(unannotated_titles))
                processed = i

    # Close any thread-local sqlite connections
    conn = getattr(_thread_local, 'conn', None)
    if conn:
        try:
            conn.close()
        except Exception:
            pass

    logging.info('Finished. Processed %d unannotated sequences. Results written to %s', processed, outpath)


if __name__ == '__main__':
    main()
