from .convert import *


class StandaloneExporter(DarwinExporter):

    def __init__(self, root, name):
        os.environ['DARWIN_BROWSERDATA_PATH'] = root
        super(StandaloneExporter, self).__init__(name)


    def add_homologs(self):
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'homologs' not in rel_node_for_genome:
                cache_file = os.path

    def get_version(self):
        #TODO: obtain real version
        return "OmaStandalone; 1.0.x"

