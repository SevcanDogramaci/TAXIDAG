import pysam
import pathlib
import pandas as pd 

def write_variant_calls_to_file():
    f = open(f"{current_dir}/data_out/1st_sample/snpCalls.txt", "w")
    for call in snpCalls:
        f.write(call)
    f.close()

def filter_indels(variant_calls):
    import re

    variant_call_ids_to_filter = []
    alt_base_column = []
    for i, read_result in variant_calls['READ_RESULTS'].iteritems():

        # search indels
        ins_search = "\+[0-9]+[ACGTNacgtn]+"
        del_search = "-[0-9]+[ACGTNacgtn]+"
        ins = re.search(ins_search, read_result)
        dels = re.search(del_search, read_result)

        if ins == None and dels == None:
            bases = "[ACGTNacgtn]"

            # get alt base by searching
            alt = re.search(bases, read_result)
            if alt != None:
                # get first alt base ?
                alt = alt.group(0)
            else:
                alt = '.'
            alt_base_column.append(alt)
            print(i, alt, variant_calls["POS"][i])
        else:
            row_id_to_delete = variant_calls.index[i]
            variant_call_ids_to_filter.append(row_id_to_delete)
    
    print("\nLength of variant calls before: ", len(variant_calls), \
     " Length of variant calls to delete: ", len(variant_call_ids_to_filter))
    variant_calls = variant_calls.drop(variant_call_ids_to_filter,inplace=False)
    print("Length of variant calls after: ", len(variant_calls), "\n")
    variant_calls.insert(len(variant_calls.columns), 'ALT', alt_base_column)
    return variant_calls
    

def filter_columns(variant_calls, columns_to_filter):
    return variant_calls.drop(columns_to_filter, axis=1, inplace=False)

def filter_variant_calls(variant_calls):
    from io import StringIO

    # convert variant calls to dataframe for easier data manipulation
    variant_calls = StringIO(variant_calls)
    variant_call_columns = ['CHROM', 'POS', 'REF', 'READ_COUNT', 'READ_RESULTS', 'READ_QUALITY']
    variant_calls_dataframe = pd.read_table(variant_calls, names=variant_call_columns)

    print("\n\n>>> FILTER INSERTION & DELETIONS <<<")
    variant_calls_dataframe = filter_indels(variant_calls_dataframe)
    print(variant_calls_dataframe.head(20))
    
    print("\n\n>>> FILTER COLUMNS <<<")
    columns_to_delete = ['READ_COUNT', 'READ_QUALITY', 'READ_RESULTS']
    variant_calls_dataframe = filter_columns(variant_calls_dataframe, columns_to_delete)
    print(variant_calls_dataframe.head(20))

    return variant_calls_dataframe

def filter_post_mortem_transitions(variant_calls):
    
    postmortem_transitions = [['T', 'C'], ['A', 'G']]
    post_mortem_trans_row_ids = []

    for i, variant_call in variant_calls.iterrows():
        ref = variant_call['REF'].upper()
        alt = variant_call['ALT'].upper()

        print("Ref:", ref, " Alt:", alt)
    
        for transition in postmortem_transitions:
            if (ref == transition[0] and alt == transition[1]) or \
                (alt == transition[0] and ref == transition[1]):
                    post_mortem_trans_row_ids.append(i)
                    print("Postmortem transition found - ", "Pos:", variant_call["POS"], "i:", i)
        print()
    return variant_calls.drop(post_mortem_trans_row_ids, inplace=False)

# --- MAIN PROGRAM STARTS ---

# file paths
current_dir = pathlib.Path(__file__).parent.absolute()

aligned_sample_to_oar_file = f"{current_dir}/data/tmp/1st_sample/o_dp_tps062_MT.bam"
sorted_aligned_sample_to_oar_file = f"{current_dir}/data_out/1st_sample/sorted_o_dp_tps062_MT.bam"
oar_ref_file = f"{current_dir}/data/Ovis_aries.Oar_v3.1.dna_rm.chromosome.MT.fa"
transv_poly_oar_file = f"{current_dir}/data/tmp/trvposoar.bed"

# samtools faidx sequence1.fasta
pysam.faidx(oar_ref_file)

# sort and index the sample's aligned sequence to the sheep's mtDNA reference
pysam.sort("-o", sorted_aligned_sample_to_oar_file, aligned_sample_to_oar_file)
pysam.index(sorted_aligned_sample_to_oar_file)

# bcftools mpileup -B -f fastafile -R bedfile bamfile | 
# bcftools call -mV indels -A --ploidy 1 -o tr_samp_oar.vcf
print("Snp calling is starting ...")
variant_calls = pysam.mpileup("-f", oar_ref_file, "-B", "-l", transv_poly_oar_file, sorted_aligned_sample_to_oar_file)
print("Snp calling finished ...")

# Open comment if you want to see variant calls in a txt file
# write_variant_calls_to_file()

variant_calls = filter_variant_calls(variant_calls)

# awk '!( $3 == "T" && $4 == "C" || $3 == "C" && $4 == "T" || 
#         $3 == "A" && $4 == "G" || $3 == "G" && $4 == "A" )
print("\n\n>>> FILTER POSTMORTEM TRANSITIONS <<<")
variant_calls = filter_post_mortem_transitions(variant_calls)
print(variant_calls.head(20))