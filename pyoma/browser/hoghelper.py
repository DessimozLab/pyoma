import collections
import logging
import re
import numpy

logger = logging.getLogger(__name__)


class HogLevelFilter(object):
    """Filter sub-hogs at given level (or below)

    This class allows to identify  at a given level or - if a hog
    does not go deep enough to cover the desired level, the whole hog will be returned.
    """

    def __init__(self, db, level):
        self.db = db
        self.level = level
        tax_below_level = self.db.tax.get_subtaxonomy_rooted_at(level)
        self.sub_level_distances = self._build_distance_to_query_level_lookup(
            tax_below_level
        )
        self.parent_levels = set(
            db.tax.get_parent_taxa(
                db.tax.get_taxnode_from_name_or_taxid(level)[0]["NCBITaxonId"]
            )["Name"]
        )
        # cache all the families nrs that are defined in the clade of interest
        self.fams_in_clade = set(
            int(roothog["Fam"])
            for roothog in db.get_hdf5_handle().root.HogLevel.where(
                '~contains(ID, b".") & (IsRoot == True)'
            )
            if roothog["Level"] in self.sub_level_distances
            or roothog["Level"] in self.parent_levels
        )

    def _build_distance_to_query_level_lookup(self, tax):
        nodes = {}

        def traverse(node, dist=0):
            nodes[node["name"].encode("utf-8")] = dist
            try:
                for c in node["children"]:
                    traverse(c, dist + 1)
            except KeyError:
                pass

        traverse(tax.as_dict())
        return nodes

    def find_closest_children_level(self, levs):
        logger.debug(levs)
        min_dist, closest = 1000, None

        for lev in levs:
            try:
                if self.sub_level_distances[lev] < min_dist:
                    closest = lev.item().decode()
                    min_dist = self.sub_level_distances[lev]
            except KeyError:
                pass
        logger.debug("closest level {} at {}".format(closest, min_dist))
        return closest

    def analyse_families(self, it):
        for fam in it:
            if fam not in self.fams_in_clade:
                continue
            subhog_ids = self.db.get_subhogids_at_level(fam_nr=fam, level=self.level)
            if len(subhog_ids) == 0:
                # let's check if we have any level in the family that is a sub clade
                levs = self.db.hog_levels_of_fam(fam_nr=fam)
                closest_level = self.find_closest_children_level(levs)
                if closest_level is not None:
                    yield (self.db.format_hogid(fam), closest_level)
            else:
                for subhog_id in subhog_ids:
                    yield (subhog_id, self.level)


def build_hog_to_og_map(db):
    hog_to_og_cnts = collections.defaultdict(collections.Counter)
    entries = db.get_hdf5_handle().get_node("/Protein/Entries")
    for protein in entries.where("(OmaHOG != b'') & (OmaGroup > 0)"):
        hog_id = protein["OmaHOG"]
        og = int(protein["OmaGroup"])
        for p in re.finditer(br"\.", hog_id):
            cur_id = hog_id[: p.start()]
            hog_to_og_cnts[cur_id].update((og,))
        hog_to_og_cnts[hog_id].update((og,))
    for hog_id, cnts in hog_to_og_cnts.items():
        yield hog_id, cnts.most_common(1)


def compare_levels(
    parent_level_hogs: numpy.array,
    children_level_hogs: numpy.array,
    return_duplication_events=False,
):
    """compares hogs at two levels and returns duplicated/lost/gained/identical states

    The function requires that the hogs are sorted numpy arrays according to their HOGid"""
    annotated = []
    dups = 0
    i = j = 0
    while i < len(parent_level_hogs) and j < len(children_level_hogs):
        if parent_level_hogs[i]["Fam"] < children_level_hogs[j]["Fam"]:
            annotated.append(tuple(parent_level_hogs[i]) + (b"lost",))
            i += 1
        elif parent_level_hogs[i]["Fam"] > children_level_hogs[j]["Fam"]:
            annotated.append(tuple(children_level_hogs[j]) + (b"gained",))
            j += 1
        else:
            if parent_level_hogs[i]["ID"] == children_level_hogs[j]["ID"]:
                annotated.append(tuple(children_level_hogs[j]) + (b"retained",))
                i += 1
                j += 1
            else:
                dupl_fnd = 0
                while j < len(children_level_hogs) and children_level_hogs[j][
                    "ID"
                ].startswith(parent_level_hogs[i]["ID"]):
                    annotated.append(tuple(children_level_hogs[j]) + (b"duplicated",))
                    j += 1
                    dupl_fnd += 1
                if dupl_fnd > 0:
                    # parent duplicated into at least one subhog
                    dups += 1
                else:
                    # parent subhog does not exist, it is lost
                    annotated.append(tuple(parent_level_hogs[i]) + (b"lost",))
                i += 1

    while i < len(parent_level_hogs):
        annotated.append(tuple(parent_level_hogs[i]) + (b"lost",))
        i += 1
    while j < len(children_level_hogs):
        annotated.append(tuple(children_level_hogs[j]) + (b"gained",))
        j += 1
    dt = children_level_hogs.dtype.descr + [("Event", "S10")]
    diff = numpy.array(annotated, dtype=dt)
    if return_duplication_events:
        return diff, dups
    return diff
