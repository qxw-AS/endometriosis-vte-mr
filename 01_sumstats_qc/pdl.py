#!/usr/bin/env python3
"""Parallel HTTP range downloader (real data only)."""
import os, sys, time, threading
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

PROXIES = {"http": os.environ.get("HTTP_PROXY"), "https": os.environ.get("HTTPS_PROXY")}

def head_size(url):
    r = requests.head(url, proxies=PROXIES, timeout=60, allow_redirects=True)
    r.raise_for_status()
    n = r.headers.get("content-length")
    if n is None:
        raise RuntimeError("no content-length")
    return int(n)

def fetch_range(url, start, end, path, tries=6):
    for k in range(tries):
        try:
            r = requests.get(url, headers={"Range": f"bytes={start}-{end}"},
                             proxies=PROXIES, timeout=300)
            if r.status_code in (200, 206) and len(r.content) == (end - start + 1):
                with open(path, "wb") as f:
                    f.write(r.content)
                return len(r.content)
            time.sleep(2 + 2 * k)
        except Exception:
            time.sleep(2 + 2 * k)
    raise RuntimeError(f"failed range {start}-{end}")

def download(url, out, workers=24, chunk=4 * 1024 * 1024):
    total = head_size(url)
    if os.path.exists(out) and os.path.getsize(out) == total:
        print(f"[skip] {out} already complete ({total})", flush=True)
        return
    tmpd = out + ".parts"
    os.makedirs(tmpd, exist_ok=True)
    ranges = [(s, min(s + chunk - 1, total - 1)) for s in range(0, total, chunk)]
    done = [0]
    lock = threading.Lock()
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fetch_range, url, s, e, os.path.join(tmpd, f"{i:06d}")): i
                for i, (s, e) in enumerate(ranges)}
        for fu in as_completed(futs):
            n = fu.result()
            with lock:
                done[0] += n
                if done[0] % (100 * 1024 * 1024) < chunk:
                    el = time.time() - t0
                    print(f"  {out}: {done[0]/1e6:.0f}/{total/1e6:.0f} MB "
                          f"({done[0]/1e6/el:.2f} MB/s)", flush=True)
    with open(out, "wb") as fo:
        for i in range(len(ranges)):
            with open(os.path.join(tmpd, f"{i:06d}"), "rb") as fi:
                fo.write(fi.read())
    size = os.path.getsize(out)
    if size != total:
        raise RuntimeError(f"size mismatch {size} != {total}")
    for i in range(len(ranges)):
        os.remove(os.path.join(tmpd, f"{i:06d}"))
    os.rmdir(tmpd)
    print(f"[done] {out} {size} bytes in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    download(sys.argv[1], sys.argv[2], workers=int(sys.argv[3]) if len(sys.argv) > 3 else 24)
