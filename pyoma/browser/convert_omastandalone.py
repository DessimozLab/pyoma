from .convert import *


class StandaloneExporter(DarwinExporter):
    DRW_CONVERT_FILE = os.path.abspath(os.path.splitext(__file__)[0] + ".drw")

    def __init__(self, root, name, **kwargs):
        os.environ['DARWIN_BROWSERDATA_PATH'] = root
        super(StandaloneExporter, self).__init__(name, **kwargs)
        self.transformed = False

    def add_homologs(self):
        self.assert_cached_results()
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'homologs' not in rel_node_for_genome:
                pass

    def get_version(self):
        #TODO: obtain real version
        return "OmaStandalone; 1.0.x"

    def assert_cached_results(self):
        if not self.transformed:
            cache_dir = os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), 'pyoma')
            res = self.call_darwin_export("TransformDataToCache('{}');".format(
                    cache_dir))
            if res != 'success':
                raise DarwinException('could not transform data from darwin')
            self.transformed = True
            os.environ['DARWIN_NETWORK_SCRATCH_PATH'] = os.getenv('DARWIN_BROWSERDATA_PATH')

    def add_orthologs(self):
        self.assert_cached_results()
        for gs in self.h5.root.Genome.iterrows():
            genome = gs['UniProtSpeciesCode'].decode()
            rel_node_for_genome = self._get_or_create_node('/PairwiseRelation/{}'.format(genome))
            if 'VPairs' not in rel_node_for_genome:
                cache_file = os.path.join(
                    os.getenv('DARWIN_NETWORK_SCRATCH_PATH', ''),
                    'pyoma', 'vps', '{}.txt.gz'.format(genome))
                if os.path.exists(cache_file):
                    data = load_tsv_to_numpy((cache_file, 0, 0, False,))
                else:
                    # fallback to read from VPsDB
                    data = self.call_darwin_export('GetVPsForGenome({})'.format(genome))

                vp_tab = self.h5.create_table(rel_node_for_genome, 'VPairs', tablefmt.PairwiseRelationTable,
                                              expectedrows=len(data))
                if isinstance(data, list):
                    data = self._convert_to_numpyarray(data, vp_tab)
                self._write_to_table(vp_tab, data)
                vp_tab.cols.EntryNr1.create_csindex()


def import_oma_run(path, outfile):
    log = getDebugLogger()
    x = StandaloneExporter(path, outfile, logger=log)
    x.add_version()
    x.add_species_data()
    x.add_orthologs()
    x.add_proteins()


if __name__ == "__main__":
    import_oma_run('~/Repositories/OmaStandalone', 'oma.h5')