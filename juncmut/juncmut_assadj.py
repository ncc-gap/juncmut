#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import gzip


def _gtf_key_value_split(gencode_gene_file):
    if gencode_gene_file.endswith("gtf.gz"):
        return " "
    elif gencode_gene_file.endswith("gff3.gz"):
        return "="
    else:
        raise Exception("unexpected gencode gene file")


def _load_exon_lengths(gencode_gene_file, wanted_transcripts):
    """Return {transcript_id: total exonic length} for the wanted transcripts.

    The total exonic length is the spliced transcript length (sum of exon
    lengths), matching the value used to pick the representative transcript
    (longest_cds_tx) in sjclass_transcript.py.
    """
    exon_len = {}
    if not wanted_transcripts:
        return exon_len

    key_value_split = _gtf_key_value_split(gencode_gene_file)
    with gzip.open(gencode_gene_file, 'rt') as hin:
        for line in hin:
            if line.startswith("#"):
                continue
            F = line.rstrip('\n').split('\t')
            if F[2] != "exon":
                continue

            transcript_id = ""
            for item in F[8].rstrip(";").split(";"):
                kv = item.strip(" ").replace('"', '').split(key_value_split)
                if len(kv) < 2:
                    continue
                if kv[0] == "transcript_id":
                    transcript_id = kv[1]
                    break

            if transcript_id not in wanted_transcripts:
                continue

            # GTF coordinates are 1-based inclusive: length = end - start + 1
            exon_len[transcript_id] = exon_len.get(transcript_id, 0) + (int(F[4]) - int(F[3]) + 1)

    return exon_len


def _mane_rank(mane_value):
    # Lower rank wins. Mirrors the MANE_Select > MANE_Plus_Clinical > others
    # preference of sjclass_transcript.define_transcript.
    if mane_value == "MANE_Select":
        return 0
    if mane_value == "MANE_Plus_Clinical":
        return 1
    return 2


def _select_index(gene_list, mane_list, is_boundary_list, exon_len):
    """Pick the representative transcript among the boundary-matching candidates.

    Selection rule (aligned with sjclass_transcript.py):
      1. MANE tier (MANE_Select > MANE_Plus_Clinical > others)
      2. longest total exonic length
      3. transcript id (ascending) as the final tie-break
    Returns the index into the lists, or None if there is no boundary match.
    """
    best_key = None
    best_index = None
    for i, is_boundary in enumerate(is_boundary_list):
        if is_boundary not in ("s", "e"):
            continue
        transcript = gene_list[i]
        key = (_mane_rank(mane_list[i]), -exon_len.get(transcript, 0), transcript)
        if best_key is None or key < best_key:
            best_key = key
            best_index = i
    return best_index


def _collect_candidate_transcripts(input_file):
    wanted = set()
    with open(input_file) as hin:
        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            if not "alternative" in csvobj["Splicing_Class"].lower():
                continue
            if not ("5'SS" in csvobj["Splicing_Class"] or "3'SS" in csvobj["Splicing_Class"]):
                continue
            for gene_col, boundary_col in (("Gene_1", "Is_Boundary_1"), ("Gene_2", "Is_Boundary_2")):
                gene_list = csvobj[gene_col].split(';')
                is_boundary_list = csvobj[boundary_col].split(';')
                for transcript, is_boundary in zip(gene_list, is_boundary_list):
                    if is_boundary in ("s", "e"):
                        wanted.add(transcript)
    return wanted


def juncmut_assadj(input_file, output_file, gencode_gene_file):
    # Pre-pass: collect every candidate transcript so the GTF is parsed once.
    exon_len = _load_exon_lengths(gencode_gene_file, _collect_candidate_transcripts(input_file))

    with open(input_file) as hin, open(output_file, 'w') as hout:
        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', fieldnames=[
            "SJ_key", "Start_ori", "End_ori", "Created_motif", "SJ_strand", "Transcript", "Gene", "MANE", "SJ_read_count"
        ])
        csvwriter.writeheader()

        csvreader = csv.DictReader(hin, delimiter='\t')
        for csvobj in csvreader:
            if not "alternative" in csvobj["Splicing_Class"].lower():
                continue

            created_motif = ""
            if "5'SS" in csvobj["Splicing_Class"]:
                created_motif = "Donor"
            elif "3'SS" in csvobj["Splicing_Class"]:
                created_motif = "Acceptor"
            else:
                continue

            transcript = None
            symbol = None
            offset = None
            mane = None
            strand = None

            # intron pos1
            mane_list_1 = csvobj["MANE_1"].split(';')
            is_boundary_list_1 = csvobj["Is_Boundary_1"].split(';')
            gene_list_1 = csvobj["Gene_1"].split(';')
            symbol_list_1 = csvobj["Symbol_1"].split(';')
            offset_list_1 = csvobj["Offset_1"].split(';')

            index_1 = _select_index(gene_list_1, mane_list_1, is_boundary_list_1, exon_len)
            if index_1 is not None:
                transcript = gene_list_1[index_1]
                symbol = symbol_list_1[index_1]
                offset = offset_list_1[index_1]
                mane = "MANE_Select" if mane_list_1[index_1] == "MANE_Select" else "NA"
                if is_boundary_list_1[index_1] == "s":
                    strand = "-"
                else:
                    strand = "+"

            # intron pos2
            if index_1 is None:
                mane_list_2 = csvobj["MANE_2"].split(';')
                is_boundary_list_2 = csvobj["Is_Boundary_2"].split(';')
                gene_list_2 = csvobj["Gene_2"].split(';')
                symbol_list_2 = csvobj["Symbol_2"].split(';')
                offset_list_2 = csvobj["Offset_2"].split(';')

                index_2 = _select_index(gene_list_2, mane_list_2, is_boundary_list_2, exon_len)
                if index_2 is not None:
                    transcript = gene_list_2[index_2]
                    symbol = symbol_list_2[index_2]
                    offset = offset_list_2[index_2]
                    mane = "MANE_Select" if mane_list_2[index_2] == "MANE_Select" else "NA"
                    if is_boundary_list_2[index_2] == "s":
                        strand = "+"
                    else:
                        strand = "-"

            if transcript is None:
                print(csvobj)
                raise Exception("juncmut_assadj.py: Unexpected data format")

            out_csvobj = {
                "SJ_key": "%s,%d,%d" % (csvobj["SJ_1"], int(csvobj["SJ_2"]) - 1 * int(offset), int(csvobj["SJ_3"]) - 1 * int(offset)),
                "Start_ori": csvobj["SJ_2"],
                "End_ori": csvobj["SJ_3"],
                "Created_motif": created_motif,
                "SJ_read_count": csvobj["SJ_7"],
                "SJ_strand": strand,
                "Transcript": transcript,
                "Gene": symbol,
                "MANE": mane,
            }
            csvwriter.writerow(out_csvobj)


if __name__ == "__main__":
    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gencode_gene_file = sys.argv[3]

    juncmut_assadj(input_file, output_file, gencode_gene_file)
