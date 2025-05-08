#!/usr/bin/env python3

import pysam
import gzip
import csv

def annot_clinvar(input_file, output_file, clinvar_file):
    
    """
    Annotates variants found within the ClinVar dataset.
    Specific clinical signifcance level of interest is
    filtered (See INFO field CLNSIG)
    
    Additional fields:
    ID            ClinVar Variation ID
    
    ClinVar specific INFO fields:
    AF_ESP              allele frequencies from GO-ESP
    AF_EXAC             allele frequencies from ExAC
    AF_TGP              allele frequencies from TGP
    *ALLELEID            The ClinVar Allele ID
    *CLNDN               ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB
    CLNDNINCL           For included Variant : ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB 
    *CLNDISDB            Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
    CLNDISDBINCL        For included Variant: Tag-value pairs of disease database name and identifier, e.g. OMIM:NNNNNN
    CLNHGVS             Top-level (primary assembly, alt, or patch) HGVS expression.
    CLNREVSTAT          ClinVar review status for the Variation ID
    *CLNSIG              Clinical significance for this single variant
    CLNSIGINCL          Clinical significance for a haplotype or genotype that includes this variant. Reported as pairs of VariationID:clinical significance.
    CLNVC               Variant type
    CLNVCSO             Sequence Ontology id for variant type
    CLNVI               The variant's clinical sources reported as tag-value pairs of database and variant identifier
    DBVARID             nsv accessions from dbVar for the variant
    GENEINFO            Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)
    *MC                  comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence
    ORIGIN              Allele origin. One or more of the following values may be added: 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other
    RS                  dbSNP ID (i.e. rs number)
    SSR                 Variant Suspect Reason Codes. One or more of the following values may be added: 0 - unspecified, 1 - Paralog, 2 - byEST, 4 - oldAlign, 8 - Para_EST, 16 - 1kg_failed, 1024 - other
    /(naomemo:INCL means included variants, case of haplotype or genotype )
    """

    #clinvar_file = f'{path_to_db}/clinvar.vcf.gz'

    NO_DATA = "---"
    gene_pathogenic_list=[]
    with gzip.open(clinvar_file, 'rb') as tin:
        for row in tin:
            if row.decode().startswith('#'): continue
            else:
                rec, gene = NO_DATA, NO_DATA
                R = row.decode().rstrip("\n").split('\t')
                # select pathogenic
                infos = R[7].split(';')
                for info in infos:
                    #GENEINFO=AGRN:375790
                    if info.startswith("GENEINFO="):
                        geneinfo = info.replace("GENEINFO=", '')
                        gene = geneinfo.split(":")[0]
                    elif info.startswith("CLNSIG="):
                        rec = info.replace("CLNSIG=", '')
                    
                if gene != NO_DATA and "pathogenic" in rec.lower():
                    gene_pathogenic_list.append(gene)

    gene_pathogenic_set = set(gene_pathogenic_list)
    
    tb = pysam.TabixFile(clinvar_file, encoding="utf-8")
    with open(input_file, 'r') as hin, open(output_file,'w') as hout:
        csvreader = csv.DictReader(hin, delimiter='\t')

        csvwriter = csv.DictWriter(hout, delimiter='\t', lineterminator='\n', quotechar='"', fieldnames=csvreader.fieldnames + [
            "ClinVar_allele", "ClinVar_ID", "ClinVar_DN", "ClinVar_DISDB", "ClinVar_SIG", "ClinVar_MC",
            "Nearby_ClinVar_pos", "Nearby_ClinVar_allele", "Nearby_ClinVar_ID", "Nearby_ClinVar_DN", "Nearby_ClinVar_DISDB", "Nearby_ClinVar_SIG", "Nearby_ClinVar_MC",
            "Norm_ClinVar_pos", "Norm_ClinVar_allele", "Norm_ClinVar_ID", "Norm_ClinVar_DN", "Norm_ClinVar_DISDB", "Norm_ClinVar_SIG", "Norm_ClinVar_MC", "Pathogenicity"
        ])
        csvwriter.writeheader()

        for csvobj in csvreader:
            gene = csvobj["Gene"]
            (mut_chr, mut_pos, mut_ref, mut_alt) = csvobj["Mut_key"].split(',')
            mut_pos = int(mut_pos)

            rec_allele, rec_allele_id, rec_cln_dn, rec_cln_disdb, rec_cln_sig, rec_mc = NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA
            m5_rec_pos, m5_rec_allele, m5_rec_allele_id, m5_rec_cln_dn, m5_rec_cln_disdb, m5_rec_cln_sig, m5_rec_mc = NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA
            norm_rec_pos, norm_rec_allele, norm_rec_allele_id, norm_rec_cln_dn, norm_rec_cln_disdb, norm_rec_cln_sig, norm_rec_mc = NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA, NO_DATA

            # around alt_pos
            tabix_records =tb.fetch(mut_chr.replace('chr',''), mut_pos - 5, mut_pos + 5)
            
            m5_pos_list = []
            m5_allele_list = []
            m5_allele_id_list = []
            m5_cln_dn_list = []
            m5_cln_disdb_list = []
            m5_cln_sig_list = []
            m5_mc_list = []
            
            if tabix_records is None:
                tabix_records = []
            # row : [0]chr[1]position[2].[3]REF[4]ALT
            for record in tabix_records:
                R = record.rstrip('\n').split('\t')
                record_pos = int(R[1])
                record_ref = R[3]
                record_alt = R[4]
                record_infos = R[7].split(';')

                # complete match allele
                if mut_pos == record_pos and mut_alt == record_alt:
                    rec_allele = record_ref + ">" + record_alt
                    for info in record_infos:
                        if info.startswith("ALLELEID="):
                            rec_allele_id = info.replace("ALLELEID=", '')
                        elif info.startswith("CLNDN="): 
                            rec_cln_dn = info.replace("CLNDN=", '')
                        elif info.startswith("CLNDISDB="):
                            rec_cln_disdb = info.replace("CLNDISDB=", '')
                        elif info.startswith("CLNSIG="): 
                            rec_cln_sig = info.replace("CLNSIG=", '')
                        #Molecular consequence
                        elif info.startswith("MC="):
                            rec_mc = info.split('|')[1].split(',')[0]

                # +/-5 bp around mutation
                elif len(record_ref) == 1 and len(record_alt) == 1:
                    m5_pos_list.append(record_pos)
                    m5_allele = record_ref + ">" + record_alt
                    for m5_info in record_infos:
                        if m5_info.startswith("ALLELEID="):
                            m5_allele_id = m5_info.replace("ALLELEID=", '')
                        if m5_info.startswith("CLNDN="): 
                            m5_cln_dn = m5_info.replace("CLNDN=", '')
                        if m5_info.startswith("CLNDISDB="):
                            m5_cln_disdb = m5_info.replace("CLNDISDB=", '')
                        if m5_info.startswith("CLNSIG="): 
                            m5_cln_sig = m5_info.replace("CLNSIG=", '')
                        #Molecular consequence
                        if m5_info.startswith("MC="):
                            m5_mc = m5_info.split('|')[1].split(',')[0] 
                    m5_allele_list.append(m5_allele)
                    m5_allele_id_list.append(m5_allele_id)
                    m5_cln_dn_list.append(m5_cln_dn)
                    m5_cln_disdb_list.append(m5_cln_disdb)
                    m5_cln_sig_list.append(m5_cln_sig)
                    m5_mc_list.append(m5_mc)

            if len(m5_pos_list) > 0:
                m5_rec_pos = ','.join(map(str,m5_pos_list)) + ","
                m5_rec_allele = ','.join(m5_allele_list) 
                m5_rec_allele_id = ','.join(m5_allele_id_list) + ","
                m5_rec_cln_dn = ','.join(m5_cln_dn_list)
                m5_rec_cln_disdb = ','.join(m5_cln_disdb_list)
                m5_rec_cln_sig = ','.join(m5_cln_sig_list)
                m5_rec_mc = ','.join(m5_mc_list)
            
            ## within gene 
            # when multiple genes are persent or no gene, skip
            gene_rec_cln_sig = NO_DATA
            if "," in gene or gene == NO_DATA:
                gene_rec_cln_sig = NO_DATA
            elif gene in gene_pathogenic_set:
                gene_rec_cln_sig = "Pathogenic"

            rec_cln_sig_mod = rec_cln_sig.replace('Conflicting_interpretations_of_pathogenicity', '')
            m5_rec_cln_sig_mod = m5_rec_cln_sig.replace('Conflicting_interpretations_of_pathogenicity', '')

            tier = NO_DATA
            if "pathogenic" in rec_cln_sig_mod.lower():
                tier = "Tier1"
            elif "pathogenic" in m5_rec_cln_sig_mod.lower():
                tier = "Tier2"
            elif gene_rec_cln_sig == "Pathogenic":
                tier = "Tier3"

            csvobj["ClinVar_allele"] = rec_allele
            csvobj["ClinVar_ID"] = rec_allele_id
            csvobj["ClinVar_DN"] = rec_cln_dn
            csvobj["ClinVar_DISDB"] = rec_cln_disdb
            csvobj["ClinVar_SIG"] = rec_cln_sig
            csvobj["ClinVar_MC"] = rec_mc
            csvobj["Nearby_ClinVar_pos"] = m5_rec_pos
            csvobj["Nearby_ClinVar_allele"] = m5_rec_allele
            csvobj["Nearby_ClinVar_ID"] = m5_rec_allele_id
            csvobj["Nearby_ClinVar_DN"] = m5_rec_cln_dn
            csvobj["Nearby_ClinVar_DISDB"] = m5_rec_cln_disdb
            csvobj["Nearby_ClinVar_SIG"] = m5_rec_cln_sig
            csvobj["Nearby_ClinVar_MC"] = m5_rec_mc
            csvobj["Norm_ClinVar_pos"] = norm_rec_pos
            csvobj["Norm_ClinVar_allele"] = norm_rec_allele
            csvobj["Norm_ClinVar_ID"] = norm_rec_allele_id
            csvobj["Norm_ClinVar_DN"] = norm_rec_cln_dn
            csvobj["Norm_ClinVar_DISDB"] = norm_rec_cln_disdb
            csvobj["Norm_ClinVar_SIG"] = norm_rec_cln_sig
            csvobj["Norm_ClinVar_MC"] = norm_rec_mc
            csvobj["Pathogenicity"] = tier

            csvwriter.writerow(csvobj)

if __name__== "__main__":
    import sys

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    clinvar_file = sys.argv[3]
    
    annot_clinvar(input_file, output_file, clinvar_file)
