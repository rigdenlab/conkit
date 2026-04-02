def construct_seq_from_chain(chain, return_borders = True, place_holder = '?',alphabet = 'RNA'):
    #takes in a biopython chain and returns its sequence accoring to the original numbering
    from Bio.PDB.Selection import unfold_entities
    import conkit.misc.rescodes

    #set the right alphabet for translation to 1 letter codes, only load standard AA/bases as modifications are unlikly to be reliably anotated in sequence file
    if alphabet == 'RNA':
        from conkit.misc.rescodes import RNA_standard_rescodes as rescodes
    if alphabet == 'Protein':
        from conkit.misc.rescodes import PROTEIN_standard_rescodes as rescodes
    if alphabet == 'DNA':
        from conkit.misc.rescodes import DNA_standard_rescodes as rescodes


    #read in biopython chain
    residues = unfold_entities(chain, "R")

    #make dict of residues with keys the resnumber and items the rescode
    res_seq_map = {}
    for res in residues:
        res_seq_map[res.get_id()[1]] = res.get_resname()
        
    #find limits of the numbering
    residues = res_seq_map.keys()
    first_res = min(residues)
    last_res = max(residues)

    #initialize sequence
    seq = ''

    #loop form the lowest key to the highest key
    #for every number add the rescode to the sequence, if there is no item add place holder
    for r in range(first_res, last_res+1):
        if r in residues: 
            if res_seq_map[r] in rescodes.keys():  #Check if residue name is a standard AA or nucleobase
                seq += rescodes[res_seq_map[r]]
            else: seq += place_holder
        else: seq += place_holder

    #output
    if return_borders: 
        return seq, first_res, last_res
    else: 
        return seq

def get_alignment_map_dict(moving, static, return_score = False, return_both_directions = True, place_holder = '?', mode='global', open_gap_score=-11, extend_gap_score=-2):
    #make a dict mapping positions in the moving seq to ones in the dynamic seq

    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.match_score = 1.0
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    aligner.wildcard = place_holder

    try:
        alignments = list(aligner.align(static, moving))
    except ValueError as e:
        print('Needleman-Wunsch alignment failed due to:\n{}'.format(e))
        return None

    alignments.sort(key=lambda x: x.score, reverse=True)
    aligned_indices = alignments[0].aligned
    score = alignments[0].score
    alignment_dict = {}
    reverse_alignment_dict = {}

    for static_chunk, moving_chunk in zip(*aligned_indices):
        for static_index, moving_index in zip(range(*static_chunk), range(*moving_chunk)):
            alignment_dict[moving_index] = static_index
            reverse_alignment_dict[static_index] = moving_index

    if return_both_directions and return_score:
        return alignment_dict, reverse_alignment_dict, score
    elif return_both_directions:
        return alignment_dict, reverse_alignment_dict
    elif return_score:
        return alignment_dict, score
    else:
        return alignment_dict

def write_renumbered_version_of_chain_in_struct(struct_file,file_type,seq,selected_chain='',moltype='RNA'):
    # identify a single chain in structure file, adapt the numbereing to match the given sequence and write out a structure file of that chain isolated and renumbered, usable for further analysis 
    import os.path
    from Bio.PDB.Selection import unfold_entities
    from Bio.PDB.PDBIO import Select
    
    # determine path and file name of struct file to construct output file names
    loc, base_fn = os.path.split(struct_file)
    try:
        outprefix, ext = base_fn.split('.')
    except: 
        print(f'cannot split {struct_file} into name and extension')
        return []

    # read in the structure
    if file_type == 'pdb':
        from Bio.PDB.PDBParser import PDBParser
        from Bio.PDB.PDBIO import PDBIO 
        from Bio.PDB.PDBIO import Select
        parser = PDBParser()
        structure = parser.get_structure(outprefix,struct_file)
        io = PDBIO()
        io.set_structure(structure)
    elif file_type == 'mmcif':
        from Bio.PDB.MMCIFParser import MMCIFParser
        from Bio.PDB.mmcifio  import MMCIFIO
        from Bio.PDB.mmcifio import Select
        parser = MMCIFParser()
        structure = parser.get_structure(outprefix,struct_file)
        io = MMCIFIO()
        io.set_structure(structure)
    else:
        print('type of structure file was not recognized, will not renumber, will most likely crash')
        return

    sequence = seq.seq  # take only the sequence of the input seq
    model = structure[0]
    chainlist = unfold_entities(model, "C")

    # if a chain was pre selected, find that chain in the pdbfile and make the sequence alignment

    if selected_chain != '':
        for chain in chainlist:
            if chain.id == selected_chain:
                chain_seq, selected_start, selected_stop = construct_seq_from_chain(chain, place_holder = '?', alphabet=moltype)
                break 
        alignment_dict, reverse_alignment_dict = get_alignment_map_dict(chain_seq, sequence, place_holder = '?')
    # if no chain was preselected, choose the chain that best aligns to the input sequence
    else:
        score_old = 0
        alignment_dict = {}
        reverse_alignment_dict = {}
        for chain in chainlist:
            chain_seq, start, stop = construct_seq_from_chain(chain, place_holder = '?', alphabet=moltype)
            alignment_dict_new, reverse_alignment_dict_new, score = get_alignment_map_dict(chain_seq, sequence, return_score = True, place_holder = '?')
            if score >= score_old:
                score_old = score
                selected_chain = chain.id
                selected_start = start
                alignment_dict = alignment_dict_new
                reverse_alignment_dict = reverse_alignment_dict_new
        if score_old <= -100000:
            print(f'no chain in {struct_file} has sufficient sequence similarity to input, aborting')
            #return 0
        
        print(f'in {struct_file} chain {selected_chain} best aligns to the provided sequence with an alignment score of {score_old}, will continue with this one')

    # renumber each residue in the selected chain based on the alignment
    chain = model[selected_chain]
    reslist = unfold_entities(chain, "R")

    unusable_residues = []
    for res in reslist:
        resid = res.get_id()
        if resid[1] - selected_start in alignment_dict.keys():
            resid = (resid[0], alignment_dict[resid[1] - selected_start]+1, resid[2])
            res.id = resid
        else: unusable_residues.append(resid)
    
    if file_type == 'pdb':
        io = PDBIO()
    elif file_type == 'mmcif':
        io = MMCIFIO()

    io.set_structure(structure)
    out_name = os.path.join(loc,f'renumbered_{selected_chain}_{outprefix}.{ext}')

    class ChainSelect(Select):

        def accept_residue(self, residue):
            if residue.get_full_id()[-2] == selected_chain and not(residue.get_id() in unusable_residues):
                return True
            else:
                return False

    print(f'Writing out isolated and renumbered chain to {loc}/renumbered_{selected_chain}_{outprefix}.{ext}')
    io.save(out_name,ChainSelect())

    return out_name, alignment_dict, reverse_alignment_dict

