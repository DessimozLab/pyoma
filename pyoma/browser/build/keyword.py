import collections
from typing import List, Union, Set, Dict, Tuple
import logging
import re
import collections
import numpy
from tqdm import tqdm
import tables

logger = logging.getLogger(__name__)

SKIPWORDS = {
    "hmm",
    "ec",
    "aa",
    "kegg",
    "pfam",
    "tigerfam",
    "go_function",
    "go_process",
    "go_component",
    "acc",
    "pid",
    "fasta",
    "z-score",
    "goid",
    "uniprot",
    "blast",
    "waterman",
}


def get_descriptions_for_entries(h5db: tables.File, entries: numpy.ndarray) -> Dict[int, str]:
    """
    Get descriptions for a numpy array of protein entries.

    This function retrieves descriptions for each protein entry from the HDF5 database.
    It processes the descriptions to remove certain patterns and returns a dictionary
    mapping entry numbers to their cleaned descriptions.

    Parameters:
    h5db (tables.File): The HDF5 database file.
    entries (numpy.ndarray): A numpy array of protein entries.

    Returns:
    Dict[int, str]: A dictionary mapping entry numbers to their descriptions.
    """
    descriptions = {}
    desc_arr = h5db.get_node("/Protein/DescriptionBuffer")
    for entry in entries:
        desc = (
            desc_arr[entry["DescriptionOffset"] : entry["DescriptionOffset"] + entry["DescriptionLength"]]
            .tobytes()
            .decode()
        )
        desc = re.split(r"\[Source:|\(ec|GO:| COG]", desc)[0]
        desc = re.sub(r"\(ec[0-9.-]+\)", " #ec# ", desc)
        desc = re.sub(r"GO:[0-9]+", " #go# ", desc)
        desc = re.sub(r" COG[0-9]+", " #cog# ", desc)
        desc = re.sub(r"transcript_id=[^ ]+", " #transcript# ", desc)
        descriptions[entry["EntryNr"]] = desc
    return descriptions


def get_xrefs_for_entries(h5db: tables.File, entries: numpy.ndarray) -> Dict[int, Dict[str, str]]:
    """
    Get cross-references for a numpy array of protein entries.

    This function retrieves cross-references for each protein entry from the HDF5 database.
    It processes the cross-references and returns a dictionary mapping entry numbers to
    their cross-references.

    Parameters:
    h5db (tables.File): The HDF5 database file.
    entries (numpy.ndarray): A numpy array of protein entries.

    Returns:
    Dict[int, Dict[str, str]]: A dictionary mapping entry numbers to their cross-references.
    """
    xref_tab: tables.Table = h5db.get_node("/XRef")
    src_enum = xref_tab.get_enum("XRefSource")
    src_query = "|".join([f"(XRefSource == {src_enum[x]})" for x in ("SourceID", "Gene Name", "Protein Name")])

    def search_by_iter(en_nrs: numpy.array):
        en_set = set(en_nrs)
        src_set = set([src_enum[x] for x in ("SourceID", "Gene Name", "Protein Name")])
        it = xref_tab.where(
            f"(EntryNr >= en_min) & (EntryNr <= en_max)", condvars={"en_min": en_nrs.min(), "en_max": en_nrs.max()}
        )
        for row in it:
            if row["EntryNr"] in en_set and row["XRefSource"] in src_set:
                yield row.fetch_all_fields()

    def search_incr_iter(en_nrs):
        chunk = 31
        for start in tqdm(range(0, len(en_nrs), chunk)):
            en_query = "|".join([f"(EntryNr == {x})" for x in en_nrs[start : start + chunk]])
            it = xref_tab.where(f"({en_query}) & ({src_query})")
            for row in it:
                yield row.fetch_all_fields()

    getter = search_by_iter if len(entries) > 10000 else search_incr_iter
    data = numpy.fromiter(getter(entries["EntryNr"]), dtype=xref_tab.dtype)
    data.sort(order=["EntryNr", "XRefSource"])

    res, cur_en = {}, None
    for row in data:
        if row["EntryNr"] != cur_en:
            if cur_en is not None:
                res[cur_en] = cur_dict
            cur_en = row["EntryNr"]
            cur_dict = collections.defaultdict(str)
        src = src_enum(row["XRefSource"])
        if src in cur_dict:
            cur_dict[src] += "; " + row["XRefId"].decode()
        else:
            cur_dict[src] = row["XRefId"].decode()
    if cur_en is not None:
        res[cur_en] = cur_dict
    return res


def clean_descriptions(descriptions: Dict[int, str], xrefs: Dict[int, Dict[str, str]]) -> Tuple[List[str], List[str]]:
    """
    Clean descriptions and return a list of cleaned descriptions and a list of original descriptions.

    This function processes the descriptions and cross-references to clean and format them.
    It returns a tuple containing a list of cleaned descriptions and a list of original descriptions.

    Parameters:
    descriptions (Dict[int, str]): A dictionary mapping entry numbers to their descriptions.
    xrefs (Dict[int, Dict[str, str]]): A dictionary mapping entry numbers to their cross-references.

    Returns:
    Tuple[List[str], List[str]]: A tuple containing a list of cleaned descriptions and a list of original descriptions.
    """
    original_descs = []
    descs = []
    for en, desc in descriptions.items():
        z = xrefs[en]
        de = [z[key] for key in ("Gene Name", "Protein Name") if z[key]]
        de.append(desc)
        de = " #br# ".join(de)
        de = de.replace(";", " #br# ").replace(")(", " #br# ")
        de = de.replace("proposed based on presence of conserved amino acid motif", " #br# ")
        de = de.replace("identified by", " #br# ").replace("similar to", " #br# ")
        de = re.sub(r"Derived by automated computational analysis using gene prediction method:[\s\w,]+\.", "", de)
        de = re.sub(r"(?<!\d)[.,()\[\]:](?!\d)", "", de)
        de = re.sub(r"\s{2,}", " ", de)

        original_descs.append(de)
        de = de.lower()
        descs.append(de)
    return descs, original_descs


def generate_ngrams(text, n=2):
    """
    Generate n-grams from a text.

    This function generates n-grams from the given text, excluding certain skip words.
    It returns a list of n-grams.

    Parameters:
    :param str text: The input text.
    :param int n: The n-gram size. Default is 2.

    Returns:
    List[str]: A list of n-grams.
    """
    words = [word for word in text.split(" ") if word not in SKIPWORDS]
    is_break = [word.startswith("#") for word in words]

    ngrams = []
    for i in range(len(words) - n + 1):
        if any(is_break[i : i + n]):
            continue
        ngrams.append(" ".join(words[i : i + n]))
    return ngrams


def score_kw(n, occ):
    """Score a keyword.

    This function calculates a score for a keyword based on its n-gram size and the number of occurrences.
    It returns the calculated score.

    :param n: the n-gram size
    :param occ: the number of occurrences
    """
    return occ * n**4


def find_best_keyword(descs: List[str]) -> str:
    """
    Find the best keyword from a list of descriptions.

    This function identifies the best keyword from the given list of descriptions.
    It returns the best keyword.

    Parameters:
    descs (List[str]): A list of descriptions.

    Returns:
    str: The best keyword.
    """
    n, score, topkw = 2, 0, ""
    while True:
        print(n)
        kws = collections.Counter()
        for desc in descs:
            kws.update(set(generate_ngrams(desc, n)))
        print(kws.most_common(4))
        if len(kws) == 0 or score_kw(n, kws.most_common(1)[0][1]) < score:
            break
        score = score_kw(n, kws.most_common(1)[0][1])
        topkw = kws.most_common(1)[0][0]
        n += 1
    return topkw


def get_most_common_caps(kw, orig_desc) -> str:
    """
    Get the most common capitalization of a keyword.

    This function finds the most common capitalization of the given keyword in the original descriptions.
    It returns the most common capitalization.

    Parameters:
    kw (str): The keyword.
    orig_desc (List[str]): A list of original descriptions.

    Returns:
    str: The most common capitalization of the keyword.
    """
    caps = collections.Counter()
    for desc in orig_desc:
        if kw in desc.lower():
            i = desc.lower().index(kw)
            caps[desc[i : i + len(kw)]] += 1
    return caps.most_common(1)[0][0]


def collect_keywords_from_entries(h5db: tables.File, entries: numpy.ndarray) -> str:
    """
    Collect keywords from a numpy array of protein entries.

    :param h5db: The HDF5 database file.
    :type h5db: tables.File
    :param entries: Protein Entries forming a group for which one keyword is to be found
    :type entries: numpy.ndarray
    :return: keyword describing the group of protein entries best.
    :rtype: str
    """
    descriptions = get_descriptions_for_entries(h5db, entries)
    xrefs = get_xrefs_for_entries(h5db, entries)
    idc = collections.Counter(x["SourceID"] for x in xrefs.values() if len(x["SourceID"]) > 0)
    best_id = "-"
    if len(idc) > 0 and idc.most_common(1)[0][1] > 1:
        best_id = idc.most_common(1)[0][0]

    clean_desc, orig_desc = clean_descriptions(descriptions, xrefs)
    kw = find_best_keyword(clean_desc)
    if len(kw) < 3:
        kw = best_id
    else:
        kw = get_most_common_caps(kw, orig_desc)
    return kw
