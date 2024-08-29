from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def __reformSequence(sequence):
    if isinstance(sequence, str):  # if sequence is str, ...

        sequence = SeqRecord(Seq(sequence), id="seqence")

    return sequence


def getORFs(sequence, minimum_length=75, start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"], remove_nested=False,
            trim_trailing=False):
    """
    sequence: sequence in Biopython Seq or String format.
    minimum_length: minimum size of ORF in nucleotides, include stop_codon.
    start_codons: recognised 3-base-pair codons for initialisation. Default: ["ATG"]
    stop_codons: recognised 3-base pair codons for termination. Default: ["TAA", "TAG", "TGA"]
    remove_nested: remove all ORFs completely encased in another. Default: False
    trim_trailing: remove ORFs are the edge of the sequence that do not have a defined stop codon. Default: False
    Returns the loci of discovered ORFs in a dictionary format.
    - start: the first basepair position of the start_codon.
    - end: the last basepair position of the stop_codon. (If stop_codon's position is 58..60, the end value is 60)
    - frame: value should be 1,2,3.
    - sense: value should be "+","-".
    - length: end - start + 1, the length include stop_codon.
    - trailing: flag, indicate whether a stop codon is found at the end of sequence.
    - index: sort by length of ORFs.
    """

    # README
    # This is a description of my own logic. While it does *not* reflect the actual code implementation,
    # the *results* should be consistent. It is intended to help users gain insight into the expected results of the code.
    # For an input sequence, it is first divided into six-frame DNA sequences.
    # For each frame, an ORF scan is performed as follows:
    # 1. Look for a start codon.
    # 2. Once a start codon is found, look for a stop codon.
    # 3. Return to step 1 after finding a stop codon.
    # ps. If a stop codon is not found at the end of the sequence,
    #     the last three-letter codon of the sequence will be treated as a stop codon.
    #
    # Due to the process of the ORF scan, we can clarify some confusing behaviours:
    # 1. "remove_nested" parameter is used to handle the incorporation of ORFs that exist in different translation frames.
    #    In fact, because of the ORF scan progress, it's not possible for ORFs in the same frame to overlap each other.

    sequence = __reformSequence(sequence)

    def findSense(sequence, sense="+", start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"]):
        """
        get information about all start codons and stop codons in the forward frame of input sequence.
        """
        start_codon_positions = list()
        stop_codon_positions = list()

        # Iterate through frames
        for frame in range(3):

            for i in range(frame, len(sequence), 3):

                if sequence[i: i + 3] in start_codons:
                    start_codon_positions.append({"position": i + 1, "frame": frame + 1, "sense": sense})  # append potential start position {position, frame, sense}

                if sequence[i: i + 3] in stop_codons:
                    # (Change 'i + 3' to modify the "end" key of the ORF)
                    stop_codon_positions.append({"position": i + 3, "frame": frame + 1, "sense": sense})  # append potential stop position {position, frame, sense}

        return start_codon_positions, stop_codon_positions

    sequence_length = len(sequence)
    forward = str(sequence.seq).upper()
    reverse = str(sequence.reverse_complement().seq).upper()

    forward_start, forward_stop = findSense(forward, "+", start_codons=start_codons, stop_codons=stop_codons)
    reverse_start, reverse_stop = findSense(reverse, "-", start_codons=start_codons, stop_codons=stop_codons)

    all_starts = forward_start + reverse_start  # list( {position, frame, sense}, {}... )
    all_stops = forward_stop + reverse_stop

    # sort in ascending order to ensure that the preceding start_codon/stop_codon iterates first.
    all_starts.sort(key=lambda x: x["position"], reverse=False)
    all_stops.sort(key=lambda x: x["position"], reverse=False)

    # same stop_codon can be contained in an ORF only once.
    # Assuming a situation: ATG1.....ATG2......TGA
    # After sorting, ATG1..TGA will be appended to orf_list. Then, since TGA is "occupied", ATG2..TGA will not be appended to orf_list.
    # Initialize the 'occupied' key for stop codons
    for stop_codon in all_stops:
        stop_codon["occupied"] = False
    # Initialize the 'occupied' key for end of each frame in the sequence
    forward_occupied = [0, 0, 0]  # +1, +2, +3 frame
    reverse_occupied = [0, 0, 0]  # -1, -2, -3 frame

    orfs = list()

    # Iterate through all start codons
    for start_codon in all_starts:

        position = start_codon["position"]
        frame = start_codon["frame"]
        sense = start_codon["sense"]

        # flag, indicate whether a stop_codon is found at the end of seqence.
        # if found, warp_case==False.
        warp_case = True

        # Iterate through all stop codons
        for stop_codon in all_stops:

            right_frame = stop_codon["frame"] == frame  # flag, indicate whether the start and stop codon have a same translation frame.
            right_sense = stop_codon["sense"] == sense  # flag, indicate whether the start and stop codon have a same sense.

            length = stop_codon["position"] - position + 1
            right_length = length >= minimum_length  # flag, indicate whether the length is sufficient.

            # if length > 0, the stop_codon is valid, warp_case is set to False
            if right_frame and right_sense and length > 0:

                warp_case = False  # found a stop_codon before the end of seqence

                if stop_codon["occupied"]:
                    break

                # Registered ORF
                if right_length:
                    orfs.append({"start": position, "end": stop_codon["position"], "frame": frame, "sense": sense,
                                 "length": length, "trailing": False})

                stop_codon["occupied"] = True
                break

        # if a stop_codon isn't found at the end of sequence, and trim_trailing==False, do...
        if warp_case and not trim_trailing:

            # Calculate the true length, after trimming non-intact codon
            length = (sequence_length - position + 1) // 3 * 3
            right_length = length >= minimum_length

            if right_length:
                # the "end position" of this orf is set to "-1", so that reset the position according to sense and frame below.
                if sense == '+':
                    if forward_occupied[frame - 1] == 0:  # frame +1 -> forward_occupied[0] ..., if != 0, it is occupied.
                        forward_occupied[frame - 1] += 1
                        orfs.append({"start": position, "end": -1, "frame": frame, "sense": sense, "length": length, "trailing": True})
                elif sense == '-':
                    if reverse_occupied[frame - 1] == 0:  # frame -1 -> reverse_occupied[0] ..., if != 0, it is occupied.
                        reverse_occupied[frame - 1] += 1
                        orfs.append({"start": position, "end": -1, "frame": frame, "sense": sense, "length": length, "trailing": True})

    # Reorder the ORFs by length to assign "index" keys.
    orfs.sort(key=lambda x: x["length"], reverse=True)

    # Transform the "start/end" keys of reverse_frame or trailing ORFs.
    for orf in orfs:

        # Transform the "start" and "end" keys of reverse_frame or trailing ORFs.
        if orf["sense"] == "-":
            orf["start"] = sequence_length - orf["start"] + 1

            if orf["end"] == -1:
                orf["end"] = orf["start"] - orf["length"] + 1

            else:
                orf["end"] = sequence_length - orf["end"] + 1

        elif orf["end"] == -1:
            orf["end"] = orf["start"] + orf["length"] - 1

    # Remove nested AFTER transforming the "start/end" keys of reverse_frame or trailing ORFs.
    if remove_nested:

        unnested_orfs = list()

        for orf_1 in orfs:

            # Because the "end" keys of trailing ORFs have been set to a specified value according to translation frame,
            # an additional conversion process is needed to recognize trailing ORFs.
            # Convert the "end" key of trailing ORF to the boundary of the sequence.
            orf_1_start = orf_1["start"]
            orf_1_end = orf_1["end"]
            if orf_1["trailing"]:
                if orf_1_end <= 3:
                    orf_1_end = 1
                elif orf_1_end >= len(sequence) - 2:
                    orf_1_end = len(sequence)

            # Because the 'start' of a sense'+' ORF corresponds to the 'end' of a sense'-' ORF,
            # the 'start' and 'end' keys of ORFs with a '-' sense will be swapped.
            if orf_1["sense"] == '-':
                (orf_1_start, orf_1_end) = (orf_1_end, orf_1_start)

            appendable = True  # flag, indicate whether orf_1 has been fully contained

            for orf_2 in orfs:
                orf_2_start = orf_2["start"]
                orf_2_end = orf_2["end"]
                if orf_2["trailing"]:
                    if orf_2_end <= 3:
                        orf_2_end = 1
                    elif orf_2_end >= len(sequence) - 2:
                        orf_2_end = len(sequence)

                # skip the same ORF, to resolve code conflict below
                if orf_2["start"] == orf_1['start'] and orf_2["end"] == orf_1['end']:
                    continue

                # Because the 'start' of a sense'+' ORF corresponds to the 'end' of a sense'-' ORF,
                # the 'start' and 'end' keys of ORFs with a '-' sense will be swapped.
                if orf_2["sense"] == '-':
                    (orf_2_start, orf_2_end) = (orf_2_end, orf_2_start)

                # if an orf_2 is found that can fully contain orf_1, orf_1's appendable will be set False.
                if orf_2_start <= orf_1_start and orf_1_end <= orf_2_end:
                    appendable = False
                    break

            # After iterating through all orf, if no orf can fully contain orf_1, then orf_1 will be appended to unnested_orfs.
            if appendable:
                unnested_orfs.append(orf_1)

        orfs = unnested_orfs

    # Add ORF "index" key.
    for i in range(len(orfs)):
        orf = orfs[i]
        orf["index"] = i + 1

    # Return results
    return orfs

def getORFNucleotides(sequence, return_loci=False, **kwargs):
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
        start = locus["start"]
        end = locus["end"]

        if locus["sense"] == "+":
            locus["nucleotide"] = Seq(forward[start - 1: end])

        else:
            locus["nucleotide"] = Seq(reverse[sequence_length - start: sequence_length - end + 1])

        nucleotides.append(locus["nucleotide"])

    if return_loci:

        return loci

    else:

        return nucleotides


def getORFProteins(sequence, translation_table=1, return_loci=False, **kwargs):
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

    # Use the result of the "getORFNucleotides" function directly to reduce errors and make modifications easier.
    loci = getORFNucleotides(sequence, return_loci=True, **kwargs)

    proteins = list()

    for locus in loci:

        protein = locus["nucleotide"].translate(table=translation_table)

        locus["protein"] = protein
        proteins.append(protein)

    if return_loci:

        return loci

    else:

        return proteins
