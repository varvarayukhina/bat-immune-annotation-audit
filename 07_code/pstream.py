#!/usr/bin/env python3
"""Stream the first `limit_bytes` of an HTTPS file to stdout using `n` parallel range requests
(ENA throttles a single connection to ~0.6 MB/s; 6-8 parallel ranges give 3-4 MB/s).
Chunks are fetched ahead of the write pointer and emitted strictly in order, so the output is a
plain prefix of the remote file (a gzip stream prefix decompresses up to the cut).
usage: pstream.py URL limit_bytes [n_parallel] [chunk_bytes] [tmp_dir]"""
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

url, limit = sys.argv[1], int(sys.argv[2])
npar = int(sys.argv[3]) if len(sys.argv) > 3 else 8
chunk = int(sys.argv[4]) if len(sys.argv) > 4 else 32 * 1024 * 1024
tmp = Path(sys.argv[5]) if len(sys.argv) > 5 else Path("/tmp")
tmp.mkdir(parents=True, exist_ok=True)


def fetch(i: int) -> Path:
    a, b = i * chunk, min(limit, (i + 1) * chunk) - 1
    p = tmp / f"pstream_{Path(url).name}_{i}.part"
    for attempt in range(6):
        r = subprocess.run(["curl", "-s", "-f", "-m", "900", "-r", f"{a}-{b}", "-o", str(p), url], check=False)
        if r.returncode == 0 and p.exists() and p.stat().st_size == b - a + 1:
            return p
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"chunk {i} ({a}-{b}) failed after retries")


nchunks = (limit + chunk - 1) // chunk
out = sys.stdout.buffer
t0 = time.time()
with ThreadPoolExecutor(max_workers=npar) as ex:
    futs = {i: ex.submit(fetch, i) for i in range(min(npar, nchunks))}
    nxt = len(futs)
    for i in range(nchunks):
        p = futs.pop(i).result()
        if nxt < nchunks:
            futs[nxt] = ex.submit(fetch, nxt)
            nxt += 1
        with open(p, "rb") as fh:
            out.write(fh.read())
        out.flush()
        p.unlink()
        done = min(limit, (i + 1) * chunk)
        if i % 8 == 7 or i == nchunks - 1:
            el = time.time() - t0
            print(f"pstream {done / 1e9:.2f} GB of {limit / 1e9:.2f} in {el / 60:.1f} min ({done / el / 1e6:.2f} MB/s)", file=sys.stderr, flush=True)
