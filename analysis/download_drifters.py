import re
from pathlib import Path
from urllib.parse import urljoin, unquote

import requests

BASE_URL = "https://observations.ipsl.fr/aeris/whirls/data/products/DRIFTERS/"
OUT = Path.home() / "drifters"          # -> /Users/emmagurcan/drifters

s = requests.Session()
s.headers.update({
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
                  "(KHTML, like Gecko) Chrome/124.0 Safari/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9,fr;q=0.8",
    "Referer": BASE_URL,
})


def get(url):
    r = s.get(url, timeout=60)
    if r.status_code != 200:
        print(f"\n{r.status_code} for {url}\nServer said:\n{r.text[:800]}\n")
        r.raise_for_status()
    return r


def list_files(url):
    """All file URLs under url (follows subfolders)."""
    files = []
    for href in re.findall(r'href="([^"?#]+)"', get(url).text):
        full = urljoin(url, href)
        if not full.startswith(url) or full == url:
            continue                     # parent dir / other sites
        files += list_files(full) if full.endswith("/") else [full]
    return files


files = list_files(BASE_URL)
print(f"{len(files)} files")
for url in files:
    dest = OUT / unquote(url[len(BASE_URL):])
    if dest.exists():
        continue                         # already downloaded
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(get(url).content)
    print("got", dest.relative_to(OUT))
print("done ->", OUT)