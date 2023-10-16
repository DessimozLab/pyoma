import omataxonomy
import collections
import tables
import numpy
import csv
import os
import rapidfuzz
import pickle


t = omataxonomy.EnvReleaseTaxonomy()
gtdb2ncbi = collections.defaultdict(list)
c = t.db.execute('SELECT * from synonym where spname like "ncbi%";')
for gtdb, ncbi_text in c.fetchall():
    ncbi = int(ncbi_text[ncbi_text.index(':')+1:])
    gtdb2ncbi[gtdb].append(ncbi)

for g, n in gtdb2ncbi.items():
    if len(n) > 1: print(g, n)

ncbi2gtdb = collections.defaultdict(list)
for g, n in gtdb2ncbi.items():
    for k in n:
        ncbi2gtdb[k].append(g)
for n, g in ncbi2gtdb.items():
    if len(g)>1: print(n,g)

with open(os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), '..', 'downloads', 'oma-species.txt'), 'rt') as fh:
    r = csv.reader((l for l in fh if not l.startswith('#')), dialect='excel-tab')
    omagenomes = [row for row in r]

taxids = [int(d[1]) for d in omagenomes]
scinames = t.translate_to_names(taxids)
newdata = []
for g, sci in zip(omagenomes, scinames):
    taxid = int(g[1])
    ncbi = gtdb2ncbi[taxid][0] if taxid<0 else taxid
    newdata.append((g[0], taxid, ncbi, sci if taxid<0 else 'n/a', *g[2:]))

oma_gtdb = list(int(x[1]) for x in omagenomes if int(x[1])<0)
tree = t.get_topology(oma_gtdb, intermediate_nodes=True, annotate=True)
oma_ncbi = list(x[2] for x in newdata if x[1]<0)
oma_ncbi_all = list(x[2] for x in newdata)
tree2 = t.get_topology(oma_ncbi, intermediate_nodes=True, annotate=True)
tree2_all = t.get_topology(oma_ncbi_all, intermediate_nodes=True, annotate=True)

taxid2code = {x[1]: x[0] for x in newdata}
taxid2code.update({x[2]: x[0] for x in newdata})

def get_node_to_speciesSet(tree, name2sp):
    res = {}
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            try:
                spp = set([name2sp[int(n.name)]])
            except KeyError:
                spp = set([])
            n.add_feature('spset', spp)
            res[int(n.name)] = n.spset
        else:
            spset = set()
            for x in n.get_children(): 
                spset |= x.spset
            n.add_feature('spset', spset)
            res[int(n.name)] = spset
    return res
gtdb_nodes_to_speciesSet = get_node_to_speciesSet(tree, taxid2code)
ncbi_nodes_to_speciesSet = get_node_to_speciesSet(tree2, taxid2code)
ncbi_all_nodes_to_speciesSet = get_node_to_speciesSet(tree2_all, taxid2code)

maxpairs = {}
for gtdb, set1 in gtdb_nodes_to_speciesSet.items():
    res = []
    for ncbitaxid, set2 in ncbi_nodes_to_speciesSet.items():
        jac = 2*len(set2.intersection(set1)) / (len(set2) + len(set1))
        res.append((ncbitaxid, jac))
    res.sort(key=lambda x: -x[1])
    maxpairs[gtdb] = res[:50]

maxpairs_ncbi = {}
for ncbi, set1 in ncbi_all_nodes_to_speciesSet.items():
    res = []
    for gtdb_taxid, set2 in gtdb_nodes_to_speciesSet.items():
        jac = 2*len(set2.intersection(set1)) / (len(set2) + len(set1))
        res.append((gtdb_taxid, jac))
    res.sort(key=lambda x: -x[1])
    maxpairs_ncbi[ncbi] = res[:50]


with tables.open_file(os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), 'OmaServer.h5'), 'r') as h5:
    pyomatax = h5.root.Taxonomy.read()

code2taxids = collections.defaultdict(list)
for k, v in taxid2code.items():
    code2taxids[v].append(k)

ancestral, extant = [], []
for gtdbid, cand in maxpairs.items():
    if gtdbid in taxid2code:
        this_one = {'gtdb_taxid': gtdbid,'gtdb_acc': t.translate_to_names([gtdbid])[0]}
        ncbi_this  = [xx for xx in code2taxids[taxid2code[gtdbid]] if xx>0]
        if len(ncbi_this) > 0:
            ncbi_this = ncbi_this[0]
            this_one['ncbi_taxid'] = ncbi_this
            this_one['sciname'] = t.translate_to_names([ncbi_this])[0]
        extant.append(this_one)
    else:
        this_one = {'gtdb_taxid': gtdbid, 'gtdb_name': t.translate_to_names([gtdbid])[0]}
        if this_one['gtdb_name'] in ('d__Bacteria', 'd__Archaea'):
            this_one['valid_oma_tax'] = False
            if this_one['gtdb_name'] == 'd__Bacteria':
                this_one['representative'] = 2
            elif this_one['gtdb_name'] == 'd__Archaea':
                this_one['representative'] = 2157
            ancestral.append(this_one)
            continue

        if gtdbid not in pyomatax['NCBITaxonId']:
            this_one['valid_oma_tax'] = False
            tnode = tree.search_nodes(taxid=gtdbid)[0]
            while tnode.taxid not in pyomatax['NCBITaxonId']:
                tnode = tnode.children[0]
            this_one['representative'] = tnode.taxid
            print(f"{gtdbid} - {this_one['gtdb_name']} not in pyomatax -> {tnode.taxid} - {t.translate_to_names([tnode.taxid])[0]}")
        else:
            this_one['valid_oma_tax'] = True
        for k in range(len(cand)):
            if cand[k][1] < 0.4*cand[0][1]: break
        rel_cand = [(x[0], t.translate_to_names([x[0]])[0], x[1], rapidfuzz.fuzz.ratio(t.translate_to_names([x[0]])[0], this_one['gtdb_name'])) for x in cand[:k]]
        rel_cand.sort(key=lambda e: (-e[2], -e[3]))
        this_one['ncbi_candidates'] = rel_cand
        ancestral.append(this_one)
        

for ncbiid, cand in maxpairs_ncbi.items():
    if ncbiid in taxid2code:
        this_one = {'ncbi_taxid': ncbiid, 'sciname': t.translate_to_names([ncbiid])[0]}
        gtdb_this  = [xx for xx in code2taxids[taxid2code[ncbiid]] if xx<0]
        if len(gtdb_this) > 0:
            gtdbid = gtdb_this[0]
            this_one['gtdb_taxid'] = gtdbid
            this_one['gtdb_acc'] = t.translate_to_names([gtdbid])[0]
        extant.append(this_one)
    else:
        this_one = {'ncbi_taxid': ncbiid, 'sciname': t.translate_to_names([ncbiid])[0]}
        
        if 2759 in t.get_lineage(ncbiid) or ncbiid in (2, 2157):
            if ncbiid not in pyomatax['NCBITaxonId']:
                this_one['valid_oma_tax'] = False
                tnode = tree2_all.search_nodes(taxid=ncbiid)[0]
                while tnode.taxid not in pyomatax['NCBITaxonId']:
                    try:
                        tnode = tnode.children[0]
                    except IndexError:
                        print(f"Cannot find a genome for {ncbiid} - {this_one['sciname']}")
                        break
                if tnode.taxid not in pyomatax['NCBITaxonId']: 
                    continue
                this_one['representative'] = tnode.taxid
                print(f"{ncbiid} - {this_one['sciname']} not in pyomatax -> {tnode.taxid} - {t.translate_to_names([tnode.taxid])[0]}")
            else:
                this_one['valid_oma_tax'] = True
            ancestral.append(this_one)
            continue  # eukaryotic genomes do not have a gtdb mapping

        for k in range(len(cand)):
            if cand[k][1] < 0.4*cand[0][1]: break
        rel_cand = [(x[0], t.translate_to_names([x[0]])[0], x[1], rapidfuzz.fuzz.ratio(t.translate_to_names([x[0]])[0], this_one['sciname'])) for x in cand[:k]]
        rel_cand.sort(key=lambda e: (-e[2], -e[3]))
        this_one['gtdb_candidates'] = rel_cand
        ancestral.append(this_one)

def unique(lst):
    seen = set()
    for e in lst:
        t = (e.get('ncbi_taxid'), e.get('gtdb_taxid'))
        if t not in seen:
            seen.add(t)
            yield e

extant = list(unique(extant))
pickle_data = {'ancestral': ancestral, 'extant': extant}

with open(os.path.join(os.getenv('DARWIN_BROWSERDATA_PATH'), 'taxonomymap.pkl'), 'wb') as fout:
    pickle.dump(pickle_data, fout)


