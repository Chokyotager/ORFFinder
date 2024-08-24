from Bio.Seq import Seq

def __reformSequence (sequence):

    if isinstance(sequence, str):

        sequence = Seq(sequence)

    return sequence

def getORFs (sequence, minimum_length=75, start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"], remove_nested=False, trim_trailing=False):

    """
    Returns the loci of discovered ORFs in a dictionary format.
    sequence: sequence in Biopython Seq or String format.
    minimum_length: minimum size of ORF in nucleotides.
    start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
    stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
    remove_nested: remove all ORFs completely encased in another. Default: False
    trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
    """

    sequence = __reformSequence(sequence)

    def findSense (sequence, sense="+", start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"]):

        start_codon_positions = list()
        stop_codon_positions = list()

        # Iterate through frames
        for frame in range(3):

            for i in range(frame, len(sequence), 3):

                if sequence[i : i + 3] in start_codons:
                    start_codon_positions.append({"position": i + 1, "frame": frame + 1, "sense": sense})

                if sequence[i : i + 3] in stop_codons:
                    stop_codon_positions.append({"position": i + 4, "frame": frame + 1, "sense": sense})

        return start_codon_positions, stop_codon_positions

    sequence_length = len(sequence)
    forward = str(sequence.seq).upper()
    reverse = str(sequence.reverse_complement().seq).upper()

    forward_start, forward_stop = findSense(forward, "+", start_codons=start_codons, stop_codons=stop_codons)
    reverse_start, reverse_stop = findSense(reverse, "-", start_codons=start_codons, stop_codons=stop_codons)

    all_starts = forward_start + reverse_start
    all_stops = forward_stop + reverse_stop

    all_starts.sort(key=lambda x: x["position"], reverse=False)
    all_stops.sort(key=lambda x: x["position"], reverse=False)

    for stop_codon in all_stops:

        stop_codon["occupied"] = False

    orfs = list()

    # Corroborate search strategy
    for start_codon in all_starts:

        position = start_codon["position"]
        frame = start_codon["frame"]
        sense = start_codon["sense"]

        warp_case = True

        for stop_codon in all_stops:

            right_frame = stop_codon["frame"] == frame
            right_sense = stop_codon["sense"] == sense

            length = stop_codon["position"] - position
            right_length = length >= minimum_length

            if right_frame and right_sense and length > 0:

                warp_case = False

                if stop_codon["occupied"]:
                    break

                # Registered ORF
                if right_length:
                    orfs.append({"start": position, "end": stop_codon["position"], "frame": frame, "sense": sense, "length": length, "trailing": False})

                stop_codon["occupied"] = True
                break

        if warp_case and not trim_trailing:

            length = sequence_length - position
            right_length = length >= minimum_length

            if right_length:
                orfs.append({"start": position, "end": -1, "frame": frame, "sense": sense, "length": length, "trailing": True})

    # Reorder by length
    orfs.sort(key=lambda x: x["length"], reverse=True)

    # Remove nested
    if remove_nested:

        unnested_orfs = list()

        for orf_1 in orfs:

            appendable = True

            for orf_2 in orfs:

                if orf_2["start"] < orf_1["start"] and orf_2["end"] > orf_1["end"] and orf_1["end"] != -1:
                    appendable = False
                    break

            if appendable:

                unnested_orfs.append(orf_1)

        orfs = unnested_orfs

    for i in range(len(orfs)):

        orf = orfs[i]
        orf["index"] = i + 1

        if orf["sense"] == "-":
            orf["start"] = sequence_length - orf["start"] + 2

            if orf["end"] == -1:
                orf["end"] = 1

            else:
                orf["end"] = sequence_length - orf["end"] + 2

        elif orf["end"] == -1:
            orf["end"] = sequence_length

    return orfs

def getORFNucleotides (sequence, return_loci=False, **kwargs):

    """
    Returns a list of Biopython Seq objects or loci of discovered ORFs with Biopython Seq objects in a dictionary format.
    sequence: sequence in Biopython Seq or String format.
    return_loci: return the loci together with the nucleotide sequences. Default: False
    minimum_length: minimum size of ORF in nucleotides. Default: 75
    start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
    stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
    remove_nested: remove all ORFs completely encased in another. Default: False
    trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
    """

    sequence = __reformSequence(sequence)

    loci = getORFs(sequence, **kwargs)

    sequence_length = len(sequence)
    forward = str(sequence.seq).upper()
    reverse = str(sequence.reverse_complement().seq).upper()

    nucleotides = list()

    for locus in loci:

        if locus["sense"] == "+":
            locus["nucleotide"] = Seq(forward[locus["start"] - 1 : locus["end"] - 1])

        else:
            locus["nucleotide"] = Seq(reverse[sequence_length - locus["start"] + 1 : sequence_length - locus["end"] + 1])

        nucleotides.append(locus["nucleotide"])

    if return_loci:

        return loci

    else:

        return nucleotides

def getORFProteins (sequence, translation_table=1, return_loci=False, **kwargs):

    """
    Returns a list of Biopython Seq objects or loci of discovered ORFs with Biopython Seq objects in a dictionary format.
    sequence: sequence in Biopython Seq or String format.
    translation_table: translation table as per BioPython. Default: 1
    return_loci: return the loci together with the protein sequences. Default: False
    minimum_length: minimum size of ORF in nucleotides. Default: 75
    start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
    stop_codons: recognised 3-base pair condons for termination. Default: ["TAA", "TAG", "TGA"]
    remove_nested: remove all ORFs completely encased in another. Default: False
    trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
    """

    sequence = __reformSequence(sequence)

    loci = getORFs(sequence, **kwargs)

    sequence_length = len(sequence)
    forward = str(sequence.seq).upper()
    reverse = str(sequence.reverse_complement().seq).upper()

    proteins = list()

    for locus in loci:

        difference = locus["length"] % 3

        if locus["sense"] == "+":
            locus["protein"] = Seq(forward[locus["start"] - 1 : locus["end"] - 1 - difference]).translate(table=translation_table)

        else:
            locus["protein"] = Seq(reverse[sequence_length - locus["start"] + 1 : sequence_length - locus["end"] + 1  - difference]).translate(table=translation_table)

        proteins.append(locus["protein"])

    if return_loci:

        return loci

    else:

        return proteins
