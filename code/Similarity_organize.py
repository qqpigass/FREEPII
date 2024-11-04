from goatools.semantic import TermCounts, get_info_content, common_parent_go_ids, min_branch_length, semantic_distance, semantic_similarity, deepest_common_ancestor
from goatools.semantic import resnik_sim, lin_sim
from goatools.associations import read_associations, dnld_assc, read_gaf
from goatools.base import get_godag
from goatools.obo_parser import GODag


#### Define Relationship Semantic Contribution 
is_a = 0.8
part_of = 0.6
regulates = 0.7
negatively_regulates = 0.7
positively_regulates = 0.7
reg = 'reg0.7'
c = 0.67 # for GOGO


#### load Gene Ontology
# go = get_godag("go-basic.obo", optional_attrs={'relationship'})


#### Load association (i.e., assc): gene to GO id mapping
# ex:
# annotations = [('human', 'goa_human.gaf'), ('yeast', 'sgd.gaf'),]
# godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)

# for species, assc_name in annotations:
#     fin_assc = os.path.join(REPO, assc_name)
#     assc_gene2gos = dnld_assc(fin_assc, godag, namespace='MF', prt=None) # Get all the MF annotations for the current species
#     ## def dnld_assc(assc_name, go2obj=None, namespace='BP', prt=sys.stdout): Download association from http://geneontology.org/gene-associations.
#     ## Example assc_name: "tair.gaf"
#     ## Download the Association
#     ## if go2obj is None: return assc_orig
#     ## If a GO DAG is provided, use only GO IDs present in the GO DAG
#     termcounts = TermCounts(godag, association)
    
#     for go_a, go_b in itertools.combinations(sorted(goids), 2):
#         ## Resnik's similarity measure is defined as the information content of the most informative common ancestor. That is, the most specific common parent-term in the GO.
#         sim_r = resnik_sim(go_a, go_b, godag, termcounts)


# # dnld_assc('goa_human.gaf', go, namespace='MF', prt=None) # download annotation from website. prt: output dir

# # associations = read_associations("association_IEA_MF.txt")
# associations = read_gaf("goa_human.gaf", namespace="MF", godag=go) # go can be None if no obo file exist
# termcounts = TermCounts(go, associations)
# print(termcounts.aspect_counts)
# print(' ')
# print(list(termcounts.annots.items())[0])
# print(' ')
# print(list(termcounts.gene2gos.items())[0])
# print(' ')
# print(list(termcounts.go2genes.items())[1000])
# print(' ')
# print(list(termcounts.gocnts.items())[1000])
# print(' ')
# print(list(termcounts.goids)[:5])


#### helper functions
def all_common_parent_go_ids(goids, godag):
    ''' This function finds the common ancestors in the GO tree of the list of goids in the input. '''
    # Find candidates from first
    rec = godag[goids[0]]
    candidates = rec.get_all_upper()
    candidates.update({goids[0]})

    # Find intersection with second to nth goid
    for goid in goids[1:]:
        rec = godag[goid]
        parents = rec.get_all_upper()
        parents.update({goid})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)
    return candidates

def lowest_common_ancestor(goterms, godag):
    ''' This function gets the nearest common ancestor using the above function. Only returns single most specific - assumes unique exists. '''
    # Take the element at maximum depth.
    return max(all_common_parent_go_ids(goterms, godag), key=lambda t: godag[t].depth)

def all_paths_to_top(term, godag):
    # inputs: term_id and Go dag with 'relationship' as optional attributes
    """ Returns all possible paths to the root node"""
    if term not in godag:
        sys.stderr.write("Term %s not found!\n" % term)
        return
    
    def _all_paths_to_top_recursive(rec):
        if rec.level == 0:
            return [[rec]]
        paths = []
        parents = rec.get_goterms_upper()
        for parent in parents:
            top_paths = _all_paths_to_top_recursive(parent)
            for top_path in top_paths:
                top_path.append(rec)
                paths.append(top_path)
        return paths
    
    go_term = godag[term]
    return _all_paths_to_top_recursive(go_term)

def all_paths_to_top_wang(term, godag, optional_relationships):
    # inputs: term_id and Go dag with 'relationship' as optional attributes
    """ Returns all possible paths to the root node"""
    if term not in godag:
        sys.stderr.write("Term %s not found!\n" % term)
        return
    def _all_paths_to_top_recursive(rec):
        if rec.level == 0:
            return [[rec]]
        paths = []
        parents = rec.get_goterms_upper_rels(optional_relationships)
        #parents = rec.get_goterms_upper()
        for parent in parents:
            #parent = go[parent1]
            top_paths = _all_paths_to_top_recursive(parent)
            for top_path in top_paths:
                top_path.append(rec)
                paths.append(top_path)
        return paths
    go_term = godag[term]
    return _all_paths_to_top_recursive(go_term)


#### Calculate semantic value
def Semantic_Value(go_id, go, method):
    ''' input: goterm_id returns all of the weighted (all relatinships in go-basic) paths to root relationship types are global variables with appropriate weights '''
    #### Wang and GOGO
    if method == 'wang':
        # calculates all paths to top (with all relationships)
        optional_relationships = {'part_of'}
        all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
        S_values = list()
        for index, path in enumerate(all_all_paths):
            S_values.append([])
            path.reverse()
            for idx, term in enumerate(path):
                if idx == 0:
                    S_values[index].append((go_id, 1)) # self's s_value = 1
                if idx < len(path)-1:
                    if term.relationship != {}:
                        if 'part_of' in term.relationship:
                            if path[idx+1] in term.relationship['part_of']:
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
                            else: 
                                S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
                        else: 
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                    else: 
                        S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
        return final_values(S_values,'max')
    
    #### GOGO: Almost Same as Wang. Only difference is that weight consider children number for each ancestor
    elif method == 'GOGO':
        optional_relationships = {'part_of'}
        all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
        S_values = list()
        for index, path in enumerate(all_all_paths):
            S_values.append([])
            path.reverse()
            for idx, term in enumerate(path):
                if idx == 0:
                    S_values[index].append((go_id, 1)) # self's s_value = 1
                if idx < len(path)-1:
                    if term.relationship != {}:
                        if 'part_of' in term.relationship:
                            if path[idx+1] in term.relationship['part_of']:
                                weight = (1 / (c + len(go[path[idx+1].item_id].children))) + part_of
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight))
                            else: 
                                weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
                        else: 
                            weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
                    else: 
                        weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
                        S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
        return final_values(S_values,'max')
    
    #### Baseline: Almost Same as Wang. Only difference is the zero weight of root ancestor and consider regulated realtionships
    else:
        all_all_paths = all_paths_to_top(go_id, go)
        S_values = list()
        for index, path in enumerate(all_all_paths):
            S_values.append([])
            path.reverse()
            for idx, term in enumerate(path):
                if idx == 0:
                    S_values[index].append((go_id, 1)) # self's s_value = 1
                if idx < len(path)-1 and path[idx+1].item_id != 'GO:0003674' and path[idx+1].item_id != 'GO:0005575' and path[idx+1].item_id != 'GO:0008150':
                    if term.relationship != {}: 
                        if 'part_of' in term.relationship:
                            if path[idx+1] in term.relationship['part_of']:
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
                            else: 
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                        elif 'regulates' in term.relationship:
                            if path[idx+1] in term.relationship['regulates']:
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
                            else: 
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                        elif 'negatively_regulates' in term.relationship:
                            if path[idx+1] in term.relationship['negatively_regulates']:
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * negatively_regulates))
                            else: 
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                        elif 'positively_regulates' in term.relationship:
                            if path[idx+1] in term.relationship['positively_regulates']:
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * positively_regulates))
                            else: 
                                S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                        else: 
                            S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
                    else: 
                        S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
                if (term.item_id == 'GO:0003674' or term.item_id == 'GO:0005575' or term.item_id == 'GO:0008150'):
                        S_values[index].append((term.item_id, 0)) # root's s_value = 0 in GOntoSim
        if method == 'Baseline':
            return final_values(S_values, 'max')
        
        #### Baseline_LCA: 
        if method == 'Baseline_LCA' :
            svaluesfinal = final_values(S_values, 'max') # edges with max weight between specified node and root
            S_values_Modified = list()
            for index, path in enumerate(S_values):
                S_values_Modified.append([])
                for idx, term1 in enumerate(path):
                    ind = [value for (term, value) in svaluesfinal if term == term1[0]]
                    # for each path in S_values, assign the max weight to the nodes (consider its weight in all paths) in path
                    S_values_Modified[index].append((path[idx][0],ind[0]))
            
            SumOfNodesOnEachPath = list()
            for index, path in enumerate(S_values_Modified):
                SumOfNodesOnEachPath.append(sum(x[1] for x in path)) # sum all the weight in path
            
            maxPath = SumOfNodesOnEachPath.index(max(SumOfNodesOnEachPath))
            return S_values_Modified[maxPath], SumOfNodesOnEachPath[maxPath] # return path with the max sum of weight and the sum

def final_values(S_values, isMax):
    ''' helper function to assign the max of the weights assigned to each term '''
    unique_terms_s_values = []
    for path in S_values:
        for term in path:
            unique_terms_s_values.append(term) # collect all (terms, weight) among paths
    unique_terms_s_values = sorted(unique_terms_s_values, key=lambda x: x[0]) # sorted based on term
    _s_values = {}
    for y, x in unique_terms_s_values: 
        if y in _s_values: 
            _s_values[y].append((y, x)) # if term in dictionary, append (terms, weight)
        else: 
            _s_values[y] = [(y, x)] # if term not in dictionary, key = term, values = (terms, weight)
    final_s_values = []
    if(isMax == 'max'):
        for node in _s_values:
            final_s_values.append(max(_s_values[node])) # append (term, max weight of this term) to final_s_values
    elif (isMax == 'min'):
        for node in _s_values:
            final_s_values.append(min(_s_values[node]))
    return final_s_values

def intersection(lst1, lst2): 
    ''' Helper Function to find intersecting terms from the two input lists of (term, s_Value) '''
    da = {v:k for v, k in lst1}
    db = {v:k for v, k in lst2} 
    return [(da[k],db[k]) for k in da.keys() & db.keys()] # return tuple composed of weight of comman node in goid_set1 and goid_set2, respectively 


#### Downward Graph
def common_children_go_ids(goids, godag):
    ''' This function finds the common children in the GO tree of the list of goids in the input. '''
    rec = godag[goids[0]]
    candidates = rec.get_all_lower() # Set the child of first goid as candidate
    candidates.update({goids[0]})
    
    for goid in goids[1:]: # for child of second to nth goid, find the intersection with the candidates, and update candidate.
        rec = godag[goid]
        children = rec.get_all_lower()
        children.update({goid})
        candidates.intersection_update(children) # intersection_update: keep if intersect else remove
    return candidates

def highest_common_descendant(goterms, godag):
    ''' This function gets the nearest common descendant using the above function. Only returns single most specific - assumes unique exists. '''
    # Take the element at minimum depth (highest).
    common_children = common_children_go_ids(goterms, godag)
    if len(common_children) != 0:
        return min(common_children, key=lambda t: godag[t].depth) # take depth attribute instead of depth to accomodate all relationships
    else:
        return 0

def all_paths_to_bottom(term, godag, x):
    ''' Inputs: Goterm_id and Go dag with 'relationship' as optional attributes, returns all possible paths to the root node '''
    if term not in godag:
        sys.stderr.write("Term %s not found!\n" % term)
        return
    
    def _all_paths_to_bottom_recursive(rec):
        if rec.depth == godag[term].depth+x:
            return [[rec]]
        else:
            paths = []
            children = rec.get_goterms_lower()
            for child in children:
                bottom_paths = _all_paths_to_bottom_recursive(child)
                for bottom_path in bottom_paths:
                    bottom_path.append(rec)
                    paths.append(bottom_path)
            return paths
    
    go_term = godag[term]
    return _all_paths_to_bottom_recursive(go_term)

def Downward_Semantic_Value(go_id, go, x):
    ''' Input: Goterm_id, returns all of the weighted nodes in path to the Goterm from the children at the xth level below.
        Relationship types are global variables with appropriate weights '''
    all_all_paths = all_paths_to_bottom(go_id, go, x)
    S_values = list()
    for index, path in enumerate(all_all_paths):
        S_values.append([])
        path.reverse()
        for idx, term in enumerate(path):
            if idx == 0:
                S_values[index].append((go_id, 1)) # self's s_value = 1
            if idx < len(path)-1:
                if term.relationship != {}:
                    if 'part_of' in term.relationship:
                        if path[idx+1] in term.relationship['part_of']:
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
                        else: 
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                    elif 'regulates' in term.relationship:
                        if path[idx+1] in term.relationship['regulates']:
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
                        else: 
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                    elif 'negatively_regulates' in term.relationship:
                        if path[idx+1] in term.relationship['negatively_regulates']:
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * negatively_regulates))
                        else: 
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                    elif 'positively_regulates' in term.relationship:
                        if path[idx+1] in term.relationship['positively_regulates']:
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * positively_regulates))
                        else: 
                            S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
                else: 
                    S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
    return final_values(S_values, 'max')


### Calculating Similarity: path-based method (Wang, GOGO, GOntoSim)
def Similarity_of_Two_GOTerms(go_id1, go_id2, go, method):
    if method == 'Baseline_LCA':
        lca = lowest_common_ancestor((go_id1, go_id2), go)
        sim_lca = Semantic_Value(lca, go, method)
        sim1 = Semantic_Value(go_id1, go, method)
        sim2 = Semantic_Value(go_id2, go, method) # sim1[1] and sim2[1] are the sums of all the nodes on the path with the max s-values. 
        sum_sim1_sim2 = (sim1[1] + sim2[1])
        return ((sim_lca[1]*2)/sum_sim1_sim2)
    
    elif method == 'GOntoSim':
        hcd = highest_common_descendant((go_id1, go_id2), go)
        
        if hcd != 0:
            hcd_depth = go[hcd].depth
            go1_depth = go[go_id1].depth
            go2_depth = go[go_id2].depth
            x = hcd_depth - go1_depth # distance to the common descendant of goid1
            y = hcd_depth - go2_depth # distance to the common descendant of goid2
            sv_a = Downward_Semantic_Value(go_id1, go, x) # path to the specified depth
            sv_b = Downward_Semantic_Value(go_id2, go, y)
            intersecting_terms = intersection(sv_a, sv_b)
            numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
            denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) # (where sv_a has 2 values for each term, second being the SV)
            sim_down = (numerator/denominator)
        else:
            sim_down = 0
        
        sim_upper = Similarity_of_Two_GOTerms(go_id1, go_id2, go, 'Baseline_LCA')
        sim = (sim_down*0.5) + (sim_upper*0.5)
        #sim = (sim_down*0.3) + (sim_upper*0.7)
        return sim
    
    elif method == 'wang' or method == 'Baseline' or method == 'GOGO':
        sv_a = Semantic_Value(go_id1, go, method)
        sv_b = Semantic_Value(go_id2, go,  method)
        intersecting_terms = intersection(sv_a, sv_b)
        numerator = sum([x for t in intersecting_terms for x in t])
        denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) # (where sv_a has 2 values for each term, second being the SV)
        Similarity = (numerator/denominator)
        return Similarity
    
    elif method == 'Baseline_Desc':
        hcd = highest_common_descendant((go_id1, go_id2), go)
        
        if hcd != 0:
            hcd_depth = go[hcd].depth
            go1_depth = go[go_id1].depth
            go2_depth = go[go_id2].depth
            x = hcd_depth - go1_depth
            y = hcd_depth - go2_depth
            sv_a = Downward_Semantic_Value(go_id1, go, x)
            sv_b = Downward_Semantic_Value(go_id2, go, y)
            intersecting_terms = intersection(sv_a,sv_b)
            numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
            denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
            sim_down = (numerator/denominator)
        else:
            sim_down = 0
        
        sim_upper = Similarity_of_Two_GOTerms(go_id1, go_id2, go, 'Baseline')
        sim = (sim_down*0.5) + (sim_upper*0.5)
        #sim = (sim_down*0.3) + (sim_upper*0.7)
        return sim
    
    elif method == 'Baseline_Desc_only':
        hcd = highest_common_descendant((go_id1, go_id2), go)
        
        if hcd != 0:
            hcd_depth = go[hcd].depth
            go1_depth = go[go_id1].depth
            go2_depth = go[go_id2].depth
            x = hcd_depth - go1_depth
            y = hcd_depth - go2_depth
            sv_a = Downward_Semantic_Value(go_id1, go, x)
            sv_b = Downward_Semantic_Value(go_id2, go, y)
            intersecting_terms = intersection(sv_a,sv_b)
            numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
            denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
            sim_down = (numerator/denominator)
        else:
            sim_down = 0
        return sim_down

def Similarity_of_Set_of_GOTerms(set1, set2, method):
    Sim1 = []
    Sim2 = []
    for idx, goterm in enumerate(set1):
        Sim1.append([])
        for goid in set2:
            Sim1[idx].append((goterm, goid,(Similarity_of_Two_GOTerms(goterm, goid, go, method))))
    for idx, goterm in enumerate(set2):
        Sim2.append([])
        for goid in set1:
            Sim2[idx].append((goterm, goid,(Similarity_of_Two_GOTerms(goterm, goid, go, method))))
    
    sem1 = []
    sem2 = []
    for index, goterm in enumerate(Sim1):
        sem1.append((max(Sim1[index], key=lambda x: x[2])))
    for index, goterm in enumerate(Sim2):
        sem2.append((max(Sim2[index], key=lambda x: x[2])))
    
    similarity = (sum(x[2] for x in sem1)+ sum(x[2] for x in sem2) )/(len(set1) + len(set2))
    return round(similarity, 3)


#### Calculating Similarity: IC-based method (Resnik, Lin)
def Similarity_of_Two_GOTerms_IC(go_id1, go_id2, go, method, termcounts):
    if method == 'Resnik':
        Similarity = resnik_sim(go_id1, go_id2, go, termcounts)
        return Similarity
    elif method == 'Lin':
        Similarity = lin_sim(go_id1, go_id2, go, termcounts)
        return  Similarity

def Similarity_of_Set_of_GOTerms_IC(set1, set2, method, termcounts):
    Sim1 = []
    Sim2 = []
    
    for idx, goterm in enumerate(set1):
        Sim1.append([])
        for goid in set2:
            if method == 'Resnik':
                Sim1[idx].append((goterm, goid,(resnik_sim(goterm, goid, go, termcounts))))
            elif method == 'Lin':
                Sim1[idx].append((goterm, goid,(lin_sim(goterm, goid, go, termcounts))))
    for idx, goterm in enumerate(set2):
        Sim2.append([])
        for goid in set1:
            if method == 'Resnik':
                Sim2[idx].append((goterm, goid,(resnik_sim(goterm, goid, go, termcounts))))
            elif method == 'Lin':
                Sim2[idx].append((goterm, goid,(lin_sim(goterm, goid, go, termcounts))))

    sem1 = []
    sem2 = []
    for index, goterm in enumerate(Sim1):
        sem1.append((max(Sim1[index], key=lambda x: x[2])))
    for index, goterm in enumerate(Sim2):
        sem2.append((max(Sim2[index], key=lambda x: x[2])))
    
    similarity = (sum(x[2] for x in sem1)+ sum(x[2] for x in sem2) )/(len(set1) + len(set2))
    return round(similarity, 3)





