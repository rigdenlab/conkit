import logging

logger = logging.getLogger(__name__)


def construct_seq_from_chain(chain, return_borders=True, place_holder='?', alphabet='RNA'):
    # Takes a Bio.PDB chain and returns its sequence according to the original numbering.
    # Only canonical residues (blank insertion code) contribute to the sequence string and
    # to the chain-position index used for alignment.  Insertion-code residues are returned
    # separately so the caller can decide how to handle them after alignment.
    #
    # Returns (seq, first_res, last_res, canonical_by_chain_pos, insertion_residues)
    # where canonical_by_chain_pos maps 0-based position in `seq` → residue object,
    # and insertion_residues is a list of residue objects whose id[2] (icode) is non-blank.
    from Bio.PDB.Selection import unfold_entities
    import conkit.misc.rescodes

    if alphabet == 'RNA':
        from conkit.misc.rescodes import RNA_standard_rescodes as rescodes
        from conkit.misc.rescodes import mod_nuclist as mod_rescodes
    elif alphabet == 'Protein':
        from conkit.misc.rescodes import PROTEIN_standard_rescodes as rescodes
        from conkit.misc.rescodes import mod_reslist as mod_rescodes
    elif alphabet == 'DNA':
        from conkit.misc.rescodes import DNA_standard_rescodes as rescodes
        from conkit.misc.rescodes import mod_nuclist as mod_rescodes
    else:
        from conkit.misc.rescodes import RNA_standard_rescodes as rescodes
        from conkit.misc.rescodes import mod_nuclist as mod_rescodes

    from conkit.misc.rescodes import IRRELEVANT_res_codes as rescodes_to_ignore

    residues = unfold_entities(chain, "R")
    # Sort by (integer seq_id, insertion code) so 92 comes before 92A before 92B.
    residues.sort(key=lambda r: (r.get_id()[1], r.get_id()[2]))

    if not residues:
        if return_borders:
            return '', 0, 0, {}, []
        return '', {}, []

    # Separate canonical residues (blank icode) from insertion-code residues.
    canonical = [r for r in residues if not r.get_id()[2].strip()]
    insertion_residues = [r for r in residues if r.get_id()[2].strip()]

    if not canonical:
        if return_borders:
            return '', 0, 0, {}, insertion_residues
        return '', {}, insertion_residues

    first_res = canonical[0].get_id()[1]
    last_res = canonical[-1].get_id()[1]

    seq = ''
    canonical_by_chain_pos = {}   # chain_pos (0-based in seq) → residue
    prev_seq_id = first_res - 1
    chain_pos = 0

    for res in canonical:
        seq_id = res.get_id()[1]
        # Fill integer-numbered gaps with placeholder characters.
        for _ in range(prev_seq_id + 1, seq_id):
            seq += place_holder
            chain_pos += 1

        resname = res.get_resname()
        if resname in rescodes:
            seq += rescodes[resname]
        elif resname in mod_rescodes:
            seq += mod_rescodes[resname]
        else:
            seq += place_holder

        canonical_by_chain_pos[chain_pos] = res
        chain_pos += 1
        prev_seq_id = seq_id

    if return_borders:
        return seq, first_res, last_res, canonical_by_chain_pos, insertion_residues
    else:
        return seq, canonical_by_chain_pos, insertion_residues


def get_alignment_map_dict(moving, static, return_score=False, return_both_directions=True, place_holder='?', mode='global', open_gap_score=-11, extend_gap_score=-2):
    # Make a dict mapping positions in the moving seq to ones in the static seq.

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
        logger.warning("Needleman-Wunsch alignment failed: %s", e)
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


def write_renumbered_version_of_chain_in_struct(struct_file, file_type, seq, selected_chain='', moltype='RNA'):
    # Identify a single chain in structure file, adapt the numbering to match the given
    # sequence and write out a structure file of that chain isolated and renumbered,
    # usable for further analysis.
    import os.path
    from Bio.PDB.Selection import unfold_entities
    from Bio.PDB.PDBIO import Select

    loc, base_fn = os.path.split(struct_file)
    try:
        outprefix, ext = base_fn.split('.')
    except Exception:
        logger.warning("Cannot split %r into name and extension; aborting renumbering.", struct_file)
        return []

    if file_type == 'pdb':
        from Bio.PDB.PDBParser import PDBParser
        from Bio.PDB.PDBIO import PDBIO
        from Bio.PDB.PDBIO import Select
        parser = PDBParser()
        structure = parser.get_structure(outprefix, struct_file)
        io = PDBIO()
        io.set_structure(structure)
    elif file_type == 'mmcif':
        from Bio.PDB.MMCIFParser import MMCIFParser
        from Bio.PDB.mmcifio import MMCIFIO
        from Bio.PDB.mmcifio import Select
        parser = MMCIFParser()
        structure = parser.get_structure(outprefix, struct_file)
        io = MMCIFIO()
        io.set_structure(structure)
    else:
        logger.warning("Unrecognised structure file type %r; cannot renumber.", file_type)
        return

    sequence = seq.seq
    model = structure[0]
    chainlist = unfold_entities(model, "C")

    if selected_chain != '':
        for chain in chainlist:
            if chain.id == selected_chain:
                chain_seq, selected_start, selected_stop, canonical_by_chain_pos, insertion_residues = \
                    construct_seq_from_chain(chain, place_holder='?', alphabet=moltype)
                break
        alignment_dict, reverse_alignment_dict = get_alignment_map_dict(
            chain_seq, sequence, place_holder='?')
    else:
        score_old = -1000
        alignment_dict = {}
        reverse_alignment_dict = {}
        canonical_by_chain_pos = {}
        insertion_residues = []
        for chain in chainlist:
            chain_seq, start, stop, cmap, icode_res = construct_seq_from_chain(
                chain, place_holder='?', alphabet=moltype)
            result = get_alignment_map_dict(chain_seq, sequence, return_score=True, place_holder='?')
            if result is None:
                continue
            alignment_dict_new, reverse_alignment_dict_new, score = result
            if score >= score_old:
                score_old = score
                selected_chain = chain.id
                selected_start = start
                alignment_dict = alignment_dict_new
                reverse_alignment_dict = reverse_alignment_dict_new
                canonical_by_chain_pos = cmap
                insertion_residues = icode_res
        if score_old == -1000:
            raise ValueError(
                f"No chain in {struct_file!r} could be aligned to the input sequence. "
                "Check that the right files are being used, or specify --chain to select a chain explicitly."
            )

        logger.info("In %r, chain %r best aligns to the input sequence (score %.1f).", struct_file, selected_chain, score_old)
        logger.debug("Chain sequence: %s", chain_seq)
        logger.debug("Input sequence: %s", sequence)

    # --- Step 1: renumber canonical (no-icode) residues ---
    # Use the chain-position index (key in canonical_by_chain_pos) rather than
    # auth_seq_id arithmetic, so gaps in author numbering are handled correctly.
    # Clear the insertion code to ' ' so downstream readers see clean integer ids.
    #
    # Canonical residues that have no FASTA counterpart (extra structural residues
    # the author gave their own seq_id) are reassigned an insertion code relative
    # to the last aligned FASTA position rather than excluded from the output.
    #
    # original_map records the provenance of every renumbered residue:
    #   (new_seq_id, new_icode) → (orig_seq_id, orig_icode, resname)
    # Callers can derive numbering anomalies from this map without needing
    # the renumbering function to classify them explicitly.
    unusable_residues = []
    claimed_fasta_positions = set()
    last_fasta_pos = -1
    icode_counter = {}  # base_seq_num → next icode ordinal (A=65, B=66, ...)
    original_map = {}

    for chain_pos, res in sorted(canonical_by_chain_pos.items()):
        original_id = res.get_id()
        resname = res.get_resname()
        if chain_pos in alignment_dict:
            fasta_pos = alignment_dict[chain_pos]
            new_seq_id = fasta_pos + 1
            res.id = (original_id[0], new_seq_id, ' ')
            claimed_fasta_positions.add(fasta_pos)
            last_fasta_pos = fasta_pos
            original_map[(new_seq_id, ' ')] = (original_id[1], original_id[2], resname)
        else:
            if last_fasta_pos >= 0:
                base_seq_num = last_fasta_pos + 1
                next_ord = icode_counter.get(base_seq_num, ord('A'))
                new_icode = chr(next_ord)
                res.id = (original_id[0], base_seq_num, new_icode)
                icode_counter[base_seq_num] = next_ord + 1
                original_map[(base_seq_num, new_icode)] = (original_id[1], original_id[2], resname)
                logger.warning(
                    "Residue %s has no FASTA counterpart; reassigned to (%d, %r) as insertion code. "
                    "This may indicate a numbering anomaly in the structure file.",
                    original_id, base_seq_num, new_icode)
            else:
                unusable_residues.append(original_id)

    # --- Step 2: handle insertion-code residues ---
    # For each insertion-code residue, look for the first unclaimed FASTA position
    # immediately after its base residue's FASTA position.
    #
    #   Free slot found  → author "misused" insertion codes; FASTA contains this
    #                      residue → assign a clean sequential number (no icode).
    #
    #   No free slot     → insertion code is used correctly (extra structural residue
    #                      not in the FASTA) → keep in output with the base's new
    #                      seq_num and the original insertion code preserved, so the
    #                      structural detail is not lost.  pdb.py will skip these
    #                      residues when computing contacts (they have no FASTA entry).
    seq_id_to_chain_pos = {res.get_id()[1]: cp for cp, res in canonical_by_chain_pos.items()}

    for res in sorted(insertion_residues, key=lambda r: (r.get_id()[1], r.get_id()[2])):
        original_id = res.get_id()
        resname = res.get_resname()
        base_seq_id = original_id[1]
        base_chain_pos = seq_id_to_chain_pos.get(base_seq_id)

        if base_chain_pos is None or base_chain_pos not in alignment_dict:
            unusable_residues.append(original_id)
            continue

        base_fasta_pos = alignment_dict[base_chain_pos]
        # Walk forward from base+1 to find the first FASTA position not yet assigned.
        candidate = base_fasta_pos + 1
        while candidate in claimed_fasta_positions:
            candidate += 1

        if candidate < len(sequence):
            # Unclaimed slot within FASTA → misused insertion code; give it a clean id.
            new_seq_id = candidate + 1
            res.id = (original_id[0], new_seq_id, ' ')
            claimed_fasta_positions.add(candidate)
            original_map[(new_seq_id, ' ')] = (original_id[1], original_id[2], resname)
            logger.debug(
                "Insertion-code residue %s reassigned to FASTA position %d (misused icode).",
                original_id, new_seq_id)
        else:
            # No free slot → correctly used insertion code; keep it with base seq_num
            # and original icode so structural detail is preserved in the output.
            base_new_seq_num = alignment_dict[base_chain_pos] + 1
            res.id = (original_id[0], base_new_seq_num, original_id[2])
            original_map[(base_new_seq_num, original_id[2])] = (original_id[1], original_id[2], resname)
            logger.debug(
                "Insertion-code residue %s kept in output as (%d, %r) (correctly used icode).",
                original_id, base_new_seq_num, original_id[2])

    if file_type == 'pdb':
        io = PDBIO()
    elif file_type == 'mmcif':
        io = MMCIFIO()

    io.set_structure(structure)
    out_name = os.path.join(loc, f'renumbered_{selected_chain}_{outprefix}.{ext}')

    class ChainSelect(Select):
        def accept_residue(self, residue):
            if residue.get_full_id()[-2] == selected_chain and residue.get_id() not in unusable_residues:
                return True
            return False

    logger.info("Writing renumbered chain to %s/renumbered_%s_%s.%s", loc, selected_chain, outprefix, ext)
    io.save(out_name, ChainSelect())

    return out_name, alignment_dict, reverse_alignment_dict, original_map


def detect_numbering_anomalies(original_map):
    """Classify renumbering events that indicate numbering anomalies in the input file.

    Parameters
    ----------
    original_map : dict
        {(new_seq_id, new_icode) → (orig_seq_id, orig_icode, resname)} as returned
        by write_renumbered_version_of_chain_in_struct.

    Returns
    -------
    list of (new_seq_id, new_icode, orig_seq_id, orig_icode, resname, anomaly_type)
    sorted by (new_seq_id, new_icode).  anomaly_type is one of:

        'MISUSED_ICODE'   — residue carried an insertion code in the input file but
                            the FASTA has a slot for it; assigned a clean sequential number.
        'EXTRA_CANONICAL' — residue had its own integer seq_id in the input file but
                            has no FASTA counterpart; reassigned an insertion code.
                            pdb.py skips these residues during contact extraction.
    """
    anomalies = []
    for (new_seq_id, new_icode), (orig_seq_id, orig_icode, resname) in original_map.items():
        if orig_icode.strip() and not new_icode.strip():
            anomalies.append((new_seq_id, new_icode, orig_seq_id, orig_icode, resname, 'MISUSED_ICODE'))
        elif not orig_icode.strip() and new_icode.strip():
            anomalies.append((new_seq_id, new_icode, orig_seq_id, orig_icode, resname, 'EXTRA_CANONICAL'))
    return sorted(anomalies, key=lambda x: (x[0], x[1]))
