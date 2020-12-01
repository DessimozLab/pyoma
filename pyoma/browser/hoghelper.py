import logging

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
        # cache all the families nrs that are defined in the clade of interest
        self.fams_in_clade = set(
            int(roothog["Fam"])
            for roothog in db.get_hdf5_handle().root.HogLevel.where(
                '~contains(ID, b".") & (IsRoot == True)'
            )
            if roothog["Level"] in self.sub_level_distances
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
