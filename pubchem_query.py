import time
import requests
import pandas as pd
import urllib.parse
from difflib import SequenceMatcher
import tqdm

# --------------------------------------------------------------------------- #
# Tenacity for automatic retries / exponential back-off
# --------------------------------------------------------------------------- #
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)

# --------------------------------------------------------------------------- #
# Small in-memory cache so we keep the properties we just fetched
# --------------------------------------------------------------------------- #
_compound_cache = {}          # cid → {'smiles': str | None, 'title': str | None}


# --------------------------------------------------------------------------- #
# Main public helper
# --------------------------------------------------------------------------- #
def best_pubchem_cid(name, *, min_similarity=0.75, suggestions=20, timeout=5):
    """Return the PubChem CID that best matches *name*, or None."""
    # 1 — quick exact/near-exact lookup
    try:
        cid = _name_to_cid(name, timeout)
    except Exception:          # network hiccup or JSON error after retries
        cid = None
    if cid:
        return cid                      # perfect hit ─ done!

    # 2 — gather autocomplete suggestions (ranked by PubChem)
    try:
        hits = _autocomplete(name, suggestions, timeout)
    except Exception:
        hits = []                       # give up gracefully if autocomplete fails
    for hit in hits:
        try:
            cid_hit = _name_to_cid(hit, timeout)
        except Exception:
            cid_hit = None
        if cid_hit:
            score = SequenceMatcher(None, name.lower(), hit.lower()).ratio()
            if score >= min_similarity:
                return cid_hit          # first suggestion above threshold

    return None                         # nothing good enough


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=0.5, min=0.5, max=4),
    retry=retry_if_exception_type(requests.RequestException),
    reraise=True,
)
def _name_to_cid(query, timeout):
    """
    Replacement for the old name→CID lookup that also records SMILES & title.
    """
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{urllib.parse.quote(query)}/property/CanonicalSMILES,Title/JSON"
    )
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    props = r.json().get("PropertyTable", {}).get("Properties", [])
    if props:
        cid = props[0].get("CID")
        smiles = props[0].get("CanonicalSMILES")
        title = props[0].get("Title") or query
        if cid:
            # Store for later use by get_pubchem_data
            _compound_cache[cid] = {"smiles": smiles, "title": title}
            return cid
    return None


@retry(
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=0.5, min=0.5, max=4),
    retry=retry_if_exception_type(requests.RequestException),
    reraise=True,
)
def _autocomplete(query, limit, timeout):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/"
        f"{urllib.parse.quote(query)}/JSON?limit={limit}"
    )
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        return r.json().get("dictionary_terms", {}).get("compound", [])
    except Exception:
        # network / API failure
        raise


def get_smiles(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url, timeout=5)
        response.raise_for_status()
        properties = response.json().get("PropertyTable", {}).get("Properties", [])
        if properties and "CanonicalSMILES" in properties[0]:
            return properties[0]["CanonicalSMILES"]
    except Exception:
        print(f"Failed to fetch SMILES for CID {cid}")
        return None
    return None


def get_pubchem_data(substance_names):
    data = []
    # for name in substance_names:
    for name in tqdm.tqdm(substance_names, desc="Fetching PubChem data"):
        # cid = attempt_get_cid(name)
        cid = best_pubchem_cid(
            name,
            # min_similarity=0.75,
            # suggestions=20,
            # timeout=5
        )
        if cid is None:
            continue

        # Fetch from cache; fall back to two-step SMILES call if needed
        record = _compound_cache.get(cid, {})
        smiles = record.get("smiles")
        resolved_name = record.get("title", name)
        if smiles is None:
            smiles = get_smiles(cid)
        if smiles is None:
            continue

        link = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
        data.append({"name": resolved_name, "SMILES": smiles, "PubChem": link})
        time.sleep(0.2)  # polite delay

    return pd.DataFrame(data)
