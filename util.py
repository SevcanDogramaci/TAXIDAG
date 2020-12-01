import pysam
import pathlib
import pandas as pd 
import pybedtools

from file_paths import current_dir, current_sample_dir

def write_to_file(data, path):
    f = open(path, "w")
    
    for item in data:
        f.write(str(item))
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
                # get first alt base 
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


def get_variants_info(data):

    new_data = [data["CHROM"], data["POS"].apply(lambda item: (item-1)), 
                data["POS"], data["ALT"]]
    new_data_headers = ["CHROM", "POS-1", "POS", "ALT"]
    dataframe = pd.concat(new_data, axis=1, keys=new_data_headers)

    print(dataframe.head(20))

    return dataframe


def find_shared_reads(chi_sample, oar_sample):
    headers = ["V1", "V2", "V3", "V4", "V5", "V6"]
    chi_btb = pd.read_table(chi_sample, sep = "\t", names = headers)
    oar_btb = pd.read_table(oar_sample,sep = "\t", names = headers)

    oar_btb["id"] = range(1, len(oar_btb)+1)
    print("\n\n>>> CHI <<<")
    print(chi_btb.head(20))
    print("\n\n>>> OAR <<<")
    print(oar_btb.head(20))

    shared2=chi_btb.merge(oar_btb, on="V4")
    shared2.sort_values(by=["id"])

    shared_chi_data = {
        "V1": shared2["V1_x"],
        "V2": shared2["V2_x"],
        "V3": shared2["V3_x"],
        "V4": shared2["V4"],
        "V5": shared2["V5_x"],
        "V6": shared2["V6_x"]
    }
    shared_chi = pd.DataFrame(shared_chi_data)

    shared_oar_data = {
        "V1": shared2["V1_y"],
        "V2": shared2["V2_y"],
        "V3": shared2["V3_y"],
        "V4": shared2["V4"],
        "V5": shared2["V5_y"],
        "V6": shared2["V6_y"]
    }
    shared_oar = pd.DataFrame(shared_oar_data)

    shared_chi.to_csv(f"{current_dir}/data_out/{current_sample_dir}/shared_chi.bed", 
                        header=False, index=False, sep="\t", mode="w")
    shared_oar.to_csv(f"{current_dir}/data_out/{current_sample_dir}/shared_oar.bed", 
                        header=False, index=False, sep="\t", mode="w")


def find_alt_freqs(all_intersects_chi, all_intersects_oar):
    all_intersects = all_intersects_chi.merge(all_intersects_oar, on="name")
    
    print("\n\n>>> ALL INTERSECTIONS <<<")
    print(all_intersects.head(5))

    all_intersects = all_intersects.query('not (Alt_y>Ref_x  &  Alt_x>Ref_y) | (Alt_y<Ref_x & Alt_x<Ref_y)')

    all_intersects.insert(    
        len(all_intersects.columns),
        "Oar_Alt_Freq", (all_intersects["Alt_y"] / all_intersects["Total_y"])
    )
    all_intersects.insert(    
        len(all_intersects.columns),
        "Chi_Alt_Freq", (all_intersects["Alt_x"] / all_intersects["Total_x"])
    )

    all_alt_freqs = pd.concat([all_intersects["Chi_Alt_Freq"], all_intersects["Oar_Alt_Freq"]], axis=1)
    return all_alt_freqs


def sort_and_index_aligned_file(sorted_aligned_sample_file, aligned_sample_file):
    pysam.sort("-o", sorted_aligned_sample_file, aligned_sample_file)
    pysam.index(sorted_aligned_sample_file)


def call_variants(ref_file, transv_poly_file, sorted_aligned_sample_file):
    # bcftools mpileup -B -f fastafile -R bedfile bamfile | 
    # bcftools call -mV indels -A --ploidy 1 -o tr_samp_oar.vcf
    print("Snp calling is starting ...")
    variant_calls = pysam.mpileup("-f", ref_file, "-B", "-l", transv_poly_file, sorted_aligned_sample_file)
    print("Snp calling finished ...")
    return variant_calls


def convert_bam_to_bed(bam_file, bam_to_bed_file):
    # bedtools bamtobed -i bamfile > samp_oarbtb.bed 
    sample_in_bam = pybedtools.example_bedtool(bam_file)
    sample_in_bed = sample_in_bam.bam_to_bed()
    write_to_file(sample_in_bed, bam_to_bed_file) 


def find_intersections(shared_file, transv_sample_file, uniq_intersections_sample_file):
    from re import sub
    
    # bedtools intersect -a sha.oar.bed -b transvoar_samp.bed -wb 
    shared = pybedtools.example_bedtool(shared_file)
    transv_poly = pybedtools.example_bedtool(transv_sample_file)
    intersects = shared.intersect(transv_poly, wb=True)
    intersects = intersects.to_dataframe()

    intersects = pd.concat([intersects["name"], intersects["blockCount"]], axis=1)

    # |awk '{gsub(/A|T|G|C/,"N",$2)}1' | awk '{gsub(/N,N/,"N",$2)}1'
    intersects["blockCount"] = intersects["blockCount"].apply(
                                        lambda item: sub("A|a|T|t|G|g|C|c|(N,N)", "N", item))

    # sort | uniq -c |awk '{ print $1,'\t',$2,'\t',$3 }'
    intersects = intersects.value_counts().reset_index(name='counts')
    intersects = intersects.sort_values(by=["name"])

    print("\n\n >>> INTERSECTIONS <<<")
    print(intersects.head(10))
    intersects.to_csv(uniq_intersections_sample_file, header=False, sep="\t", mode="w")
        
    return intersects


def insert_ref_and_alt_allele_numbers(intersects):
    intersects["id"] = range(1, len(intersects)+1)

    intersects_with_ref_allele = intersects[intersects["blockCount"] == '.']
    intersects_with_alt_allele = intersects[intersects["blockCount"] == 'N']

    intersects_with_ref_allele.insert( len(intersects_with_ref_allele.columns), "Ref", intersects_with_ref_allele["counts"])
    intersects_with_ref_allele.insert( len(intersects_with_ref_allele.columns), "Alt", [0]*len(intersects_with_ref_allele))

    intersects_with_alt_allele.insert( len(intersects_with_alt_allele.columns),"Alt", intersects_with_alt_allele["counts"])
    intersects_with_alt_allele.insert(len(intersects_with_alt_allele.columns), "Ref", [0]*len(intersects_with_alt_allele))

    all_intersects = pd.concat([intersects_with_ref_allele, intersects_with_alt_allele])
    all_intersects.insert(len(all_intersects.columns), "Total", all_intersects["Ref"] + all_intersects["Alt"])
    all_intersects = filter_columns(all_intersects, ["counts", "blockCount"])
    all_intersects = all_intersects.sort_values(by=["id"])
    
    print("\n\n>>> ALL INTERSECTIONS <<<")
    print(all_intersects.head(5))
    return all_intersects