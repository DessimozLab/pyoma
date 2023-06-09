from collections import Counter, defaultdict

import numpy
from property_manager import lazy_property
from ..browser.db import Database
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from typing import Iterable, Union
import itertools
import pandas as pd
import logging

logger = logging.getLogger(__name__)


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

    def run_study(self, study, alpha=0.05):
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

        study = set(study) & self._with_annot
        study_counts = self._get_term_counts(study)

        df = pd.DataFrame(
            do(study, study_counts),
            columns=[
                "GO_ID",
                "study_count",
                "pop_count",
                "Study_Entries",
                "study_n",
                "pop_n",
                "p_uncorrected",
            ],
        )
        if len(df) == 0:
            return None

        # corrected p values
        df["p_bonferroni"] = multipletests(
            df.p_uncorrected.values, alpha=alpha, method="bonferroni"
        )[1]
        df["p_fdr_bh"] = multipletests(
            df.p_uncorrected.values, alpha=alpha, method="fdr_bh"
        )[1]

        # results
        df["GO_Name"] = df["GO_ID"].apply(lambda t: t.name)
        # df["GO_Depth"] = df["GO_ID"].apply(lambda t: self._go[t].min_depth)
        df["GO_Aspect"] = df["GO_ID"].apply(lambda t: t.aspect)
        df["GO_ID"] = df["GO_ID"].apply(lambda t: str(t))

        df["Study_Proportion"] = df.apply(
            lambda x: (x["study_count"] / x["study_n"]), axis=1
        )
        df["Study_Ratio"] = df.apply(
            lambda x: "{} / {}".format(x["study_count"], x["study_n"]), axis=1
        )
        df["Study_Entries"] = df["Study_Entries"].apply(
            lambda x: ", ".join(sorted(map(str, x)))
        )

        df["Population_Proportion"] = df.apply(
            lambda x: (x["pop_count"] / x["pop_n"]), axis=1
        )
        df["Population_Ratio"] = df.apply(
            lambda x: "{} / {}".format(x["pop_count"], x["pop_n"]), axis=1
        )
        df["Fold_Change"] = df["Study_Proportion"] / df["Population_Proportion"]

        return df.sort_values(
            ["p_bonferroni", "Fold_Change", "p_uncorrected"],
            ascending=[True, False, True],
        ).reset_index(drop=True)


def extent_species_go_enrichment(
    db: Database, foreground_entries: Iterable[int], alpha=0.05
) -> pd.DataFrame:
    foreground = sorted(set(foreground_entries))
    rng = db.id_mapper["OMA"].genome_range(foreground[0])
    if rng[1] < foreground[-1]:
        raise ValueError("foreground_entries must all be from the same speices")
    anno = defaultdict(set)
    for row in db.get_gene_ontology_annotations(rng[0], stop=rng[1]):
        anno[int(row["EntryNr"])].add(int(row["TermNr"]))
    main_iso = [int(row["EntryNr"]) for row in db._main_isoform_iter(rng[0], rng[1])]
    goea = GOEnrichmentAnalysis(db, main_iso, annots=anno, ensure_term=True)
    logger.debug(f"go_traverse_cache: {goea._go._traverseGraph.cache_info()}")
    return goea.run_study(foreground, alpha=alpha)


def ancestral_species_go_enrichment(
    db: Database,
    level: Union[str, bytes, int],
    foreground_hogs: Iterable[Union[str, bytes]],
    alpha=0.05,
) -> pd.DataFrame:
    anc_node = db._ancestral_node(level)
    hog_ids = anc_node.Hogs.read(field="ID")
    foreground_hogs = set(
        x.encode("utf-8") if isinstance(x, str) else x for x in foreground_hogs
    )
    if not numpy.isin(foreground_hogs, hog_ids).all():
        raise ValueError("foreground_hogs must all be from the same level")
    annots = defaultdict(set)
    for row in anc_node.GeneOntology:
        annots[hog_ids[row.nrow]].add(int(row["TermNr"]))
    goea = GOEnrichmentAnalysis(db, hog_ids, annots=annots, ensure_term=True)
    return goea.run_study(foreground_hogs, alpha=alpha)
