import collections
from collections import Counter, defaultdict

import numpy
import pandas
import pathlib
from property_manager import lazy_property
from ..browser.db import Database
from ..browser.geneontology import GOAspect
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.manifold import MDS
from typing import Iterable, Union
import plotly
import plotly.graph_objects as go
import itertools
import pandas as pd
import logging

logger = logging.getLogger(__name__)
LONGONT = {"F": "MF", "P": "BP", "C": "CC"}


class GOEnrichmentAnalysis(object):
    def __init__(self, db, population, annots, most_specific=False, ensure_term=False):
        self._db = db
        self._go = db.gene_ontology

        self._annots = self._annot_coverage(annots, most_specific, ensure_term)
        self._with_annot = set(self._annots.keys())
        self._pop = set(population) & self._with_annot
        self._pop_counts = self._get_term_counts(self._pop)

    def _annot_coverage(self, annots, most_specific, ensure_term):
        if not most_specific or ensure_term:
            #  first get the most specific set for each entry
            for g in annots:
                if ensure_term:
                    annots[g] = set(map(self._go.term_by_id, annots[g]))

                z = set()
                for t in annots[g]:
                    if (len(annots[g] & self._go.get_subterms(t))) == 1:
                        z.add(t)
                annots[g] = z

        #  create vocabulary for goea
        vocab = set(itertools.chain.from_iterable(annots.values()))
        # now expand into the vocab for each term
        for g in annots:
            z = set()
            for t in annots[g]:
                z.update(self._go.get_superterms_incl_queryterm(t) & vocab)
            annots[g] = z

        return annots

    @lazy_property
    def _term_to_entry(self):
        z = defaultdict(set)
        for k, v in self._annots.items():
            for t in v:
                z[t].add(k)
        return z

    def _get_term_counts(self, ents):
        # counts a set of terms
        z = Counter()
        for t, t_ents in self._term_to_entry.items():
            n = len(t_ents & ents)
            if n > 0:
                z[t] = n
        return z

    def run_study(self, study, alpha=0.05, correction=None):
        """
        Run a GOEA with the study set as a foreground, with the already loaded background.
        """

        def do(study, study_counts):
            study_n = len(study)
            pop_n = len(self._pop)

            for t in study_counts:
                study_count = study_counts[t]
                pop_count = self._pop_counts[t]

                study_ratio = study_count / study_n
                pop_ratio = pop_count / pop_n

                # do we need to be able to compute the study_entries?
                study_entries = self._term_to_entry[t] & study

                if study_ratio >= pop_ratio:
                    # only test enrichment
                    a = study_count
                    b = study_n - study_count
                    c = pop_count - study_count
                    d = pop_n - pop_count - b

                    p = fisher_exact([[a, b], [c, d]], alternative="greater")[1]
                    yield (t, study_count, pop_count, study_entries, study_n, pop_n, p)

        if correction is None:
            correction = "fdr_bh"
        if correction not in ("fdr_bh", "bonferroni"):
            raise ValueError(f"correction parameter must be 'fdr_bh' or 'bonferroni': {correction}")
        study = set(study) & self._with_annot
        study_counts = self._get_term_counts(study)

        df = pd.DataFrame(
            do(study, study_counts),
            columns=[
                "GO_ID",
                "study_count",
                "population_count",
                "study_entries_with_go_term",
                "study_size",
                "pop_n",
                "p_uncorrected",
            ],
        )
        if len(df) == 0:
            return None

        # corrected p values
        df["p_bonferroni"] = multipletests(df.p_uncorrected.values, alpha=alpha, method="bonferroni")[1]
        df["p_fdr_bh"] = multipletests(df.p_uncorrected.values, alpha=alpha, method="fdr_bh")[1]
        if correction == "fdr_bh":
            df = df[df["p_fdr_bh"] <= alpha]
        elif correction == "bonferroni":
            df = df[df["p_bonferroni"] <= alpha]

        # results
        df["GO_name"] = df["GO_ID"].apply(lambda t: t.name)
        df["GO_depth"] = df["GO_ID"].apply(lambda t: t.min_depth)
        df["GO_aspect"] = df["GO_ID"].apply(lambda t: LONGONT[GOAspect.to_char(t.aspect)])
        df["GO_ID"] = df["GO_ID"].apply(lambda t: str(t))

        df["study_proportion"] = df.apply(lambda x: (x["study_count"] / x["study_size"]), axis=1)
        df["study_ratio"] = df.apply(lambda x: "{} / {}".format(x["study_count"], x["study_size"]), axis=1)
        df["study_entries_with_go_term"] = df["study_entries_with_go_term"].apply(
            lambda x: ", ".join(sorted(map(lambda e: e.decode() if isinstance(e, bytes) else str(e), x)))
        )

        df["population_proportion"] = df.apply(lambda x: (x["population_count"] / x["pop_n"]), axis=1)
        df["population_ratio"] = df.apply(lambda x: "{} / {}".format(x["population_count"], x["pop_n"]), axis=1)
        df["population_size"] = df["pop_n"]
        df["fold_change"] = df["study_proportion"] / df["population_proportion"]

        cols = [
            "GO_ID",
            "GO_name",
            "GO_aspect",
            "GO_depth",
            "p_uncorrected",
            "p_bonferroni",
            "p_fdr_bh",
            "study_count",
            "study_size",
            "study_ratio",
            "study_proportion",
            "population_count",
            "population_size",
            "population_ratio",
            "population_proportion",
            "fold_change",
            "study_entries_with_go_term",
        ]
        df = df[cols]

        return df.sort_values(
            [f"p_{correction}", "fold_change", "p_uncorrected"],
            ascending=[True, False, True],
        ).reset_index(drop=True)


def extant_species_go_enrichment(db: Database, foreground_entries: Iterable[int], alpha=0.05) -> pd.DataFrame:
    foreground = sorted(set(foreground_entries))
    omaid_mapper = db.id_mapper["OMA"]
    rng = omaid_mapper.genome_range(foreground[0])
    if rng[1] < foreground[-1]:
        fnd_species = collections.Counter(map(lambda enr: omaid_mapper.map_entry_nr(enr)[:5], foreground_entries))
        raise ValueError(
            f"foreground_entries must all be from the same species. "
            f"Found {len(fnd_species)} different ones (top 3: {fnd_species.most_common(3)}"
        )
    anno = defaultdict(set)
    for row in db.get_gene_ontology_annotations(rng[0], stop=rng[1] + 1):
        anno[omaid_mapper.map_entry_nr(int(row["EntryNr"]))].add(int(row["TermNr"]))
    main_iso = [omaid_mapper.map_entry_nr(int(row["EntryNr"])) for row in db._main_isoform_iter(rng[0], rng[1])]
    goea = GOEnrichmentAnalysis(db, main_iso, annots=anno, ensure_term=True)
    logger.debug("go_traverse_cache: %s", goea._go._traverseGraph.cache_info())
    return goea.run_study(map(omaid_mapper.map_entry_nr, foreground), alpha=alpha)


def ancestral_species_go_enrichment(
    db: Database,
    level: Union[str, bytes, int],
    foreground_hogs: Iterable[Union[str, bytes]],
    alpha=0.05,
    score_cutoff=0.1,
) -> pd.DataFrame:
    anc_node = db._ancestral_node(level)
    hog_ids = anc_node.Hogs.read(field="ID")
    foreground_hogs = numpy.fromiter(
        set(x.encode("utf-8") if isinstance(x, str) else x for x in foreground_hogs),
        dtype="S255",
    )
    print(foreground_hogs)
    if not numpy.isin(foreground_hogs, hog_ids).all():
        logger.warning(f"not all hogs found at the requested level: {numpy.isin(foreground_hogs, hog_ids)}")
        print(numpy.isin(foreground_hogs, hog_ids))
        raise ValueError("foreground_hogs must all be from the same level")
    annots = defaultdict(set)
    for row in anc_node.GeneOntology:
        if row["RawScore"] >= score_cutoff:
            annots[hog_ids[row["HogRow"]]].add(int(row["TermNr"]))
    goea = GOEnrichmentAnalysis(db, hog_ids, annots=annots, ensure_term=True)
    return goea.run_study(foreground_hogs, alpha=alpha)


class GO_MDS:
    def __init__(self, db):
        self.db = db

    def __call__(self, terms):
        return self.mds(terms)

    def mds(self, terms):
        """
        Perform MDS on the semantic similarity between the terms
        """
        # split the terms into different aspects first
        terms_by_aspect = defaultdict(set)
        for t in map(self.db.freq_aware_gene_ontology.term_by_id, terms):
            terms_by_aspect[t.aspect].add(t)

        # for each aspect, perform the
        mds = {}
        for aspect, subterms in terms_by_aspect.items():
            aspect_desc = GOAspect.to_string(aspect)

            subterms = sorted(map(lambda t: t.id, subterms))

            n = len(subterms)
            ss = numpy.ones((n, n), dtype=numpy.float64)
            for i in range(n):
                for j in range(i + 1, n):
                    ss[i, j] = ss[j, i] = self.db.freq_aware_gene_ontology.semantic_similarity(subterms[i], subterms[j])

            #  now perform a 2d MDS
            xy = MDS(dissimilarity="precomputed", normalized_stress="auto").fit(10 * (1 - ss)).embedding_
            mds[aspect_desc] = {t: tuple(xy[i]) for (i, t) in enumerate(subterms)}
        return mds


def generate_plots(enr_df, db, folder=None, pval_col="p_bonferroni"):
    """generate plots MDS of similarities of GO terms with significance and information content.

    The dataframes with the resulting MDS embedings is also returned."""
    if folder is None:
        folder = ""
    folder = pathlib.Path(folder)
    mds = GO_MDS(db)
    xy = mds.mds(enr_df["GO_ID"])
    per_aspect = {}
    for aspect, terms in xy.items():
        labels = sorted(terms.keys())
        df = pandas.DataFrame(
            {
                "id": labels,
                "MDS1": [terms[x][0] for x in labels],
                "MDS2": [terms[x][1] for x in labels],
            }
        )
        df["GO_term"] = df["id"].apply(db.freq_aware_gene_ontology.ensure_term)
        df["id"] = df["GO_term"].apply(str)
        df["name"] = df["GO_term"].apply(lambda t: t.name)
        df["ic"] = df["GO_term"].apply(lambda t: db.freq_aware_gene_ontology.ic(t))
        df = pandas.merge(df, enr_df[["GO_ID", pval_col]], left_on="id", right_on="GO_ID")
        df.drop(columns=["id", "GO_term"], inplace=True)

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=df["MDS1"],
                y=df["MDS2"],
                mode="markers",
                marker={
                    "size": 200 / (df["ic"] + 3),
                    "color": numpy.log10(df[pval_col]),
                    "colorbar": {"title": f"p-value ({pval_col}) (log10 scale)"},
                },
                opacity=0.9,
                customdata=numpy.stack((df["GO_ID"], df["name"], df["ic"], df[pval_col]), axis=-1),
                hovertemplate=(
                    "<b>ID</b>: %{customdata[0]}<br>"
                    + "<b>Name</b>: %{customdata[1]}<br>"
                    + "<b>Information Content</b>: %{customdata[2]:,.3f}<br>"
                    + "<b>p-value</b>: %{customdata[3]:,.3e}<br>"
                    + "<extra></extra>"
                ),
            )
        )
        fig.update_layout(
            xaxis={"title": "Semantic Space X"},
            yaxis={"title": "Semantic Space Y"},
            title=f"{aspect.title()}",
        )
        plotly.offline.plot(fig, filename=str(folder / f"{aspect}.html"))
        per_aspect[aspect] = df
    return per_aspect
