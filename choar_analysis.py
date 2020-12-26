import sys
import pysam
import pybedtools

from util import *
import file_paths

if __name__ == "__main__":
    from os import path, mkdir
    
    for arg in sys.argv:
        if arg == "-d":
            file_paths.DEBUG = True
            
    if not path.exists(file_paths.data_out_dir):
        mkdir(data_out_dir)
    if not path.exists(file_paths.current_sample_data_out_dir):
        mkdir(file_paths.current_sample_data_out_dir)
    print("Debug Mode On\n" if file_paths.DEBUG else "Debug Mode Off\n")


# samtools faidx sequence1.fasta
pysam.faidx(file_paths.oar_ref_file)
pysam.faidx(file_paths.chi_ref_file)


# === STEP 1 ===
# sort and index the sample's aligned sequence to the sheep's or goat's mtDNA reference
sort_and_index_aligned_file(file_paths.sorted_aligned_sample_to_oar_file, file_paths.aligned_sample_to_oar_file)
sort_and_index_aligned_file(file_paths.sorted_aligned_sample_to_chi_file, file_paths.aligned_sample_to_chi_file)
# ==============


# === STEP 2 ===
variant_calls_oar = call_variants(file_paths.oar_ref_file, file_paths.transv_poly_oar_file, file_paths.sorted_aligned_sample_to_oar_file)
variant_calls_chi = call_variants(file_paths.chi_ref_file, file_paths.transv_poly_chi_file, file_paths.sorted_aligned_sample_to_chi_file)
# === ====== ===


# === STEP 3 ===
variant_calls_oar = filter_variant_calls_from_pileup_format(variant_calls_oar)
variant_calls_chi = filter_variant_calls_from_pileup_format(variant_calls_chi)


# awk '!( $3 == "T" && $4 == "C" || $3 == "C" && $4 == "T" || 
#         $3 == "A" && $4 == "G" || $3 == "G" && $4 == "A" )
variant_calls_oar = filter_post_mortem_transitions(variant_calls_oar)
variant_calls_chi = filter_post_mortem_transitions(variant_calls_chi)


# awk '{print $1,$2-1,$2,$4}'
variants_info_oar = get_variants_info(variant_calls_oar, file_paths.transv_sample_oar_file)
variants_info_chi = get_variants_info(variant_calls_chi, file_paths.transv_sample_chi_file)
# === ====== ===


# === STEP 4 ===
sample_oar_in_bed = convert_bam_to_bed(file_paths.sorted_aligned_sample_to_oar_file, file_paths.sample_oar_bamtobed_file)   
sample_chi_in_bed = convert_bam_to_bed(file_paths.sorted_aligned_sample_to_chi_file, file_paths.sample_chi_bamtobed_file)
# === ====== ===


# === STEP 5 ===
# sharead.R
shared_chi, shared_oar = find_shared_reads(sample_chi_in_bed, sample_oar_in_bed)
# === ====== ===


# === STEP 6 ===
intersect_oar = find_intersections(shared_oar, variants_info_oar, file_paths.uniq_intersections_sample_oar_file)
intersect_chi = find_intersections(shared_chi, variants_info_chi, file_paths.uniq_intersections_sample_chi_file)
# === ====== ===


# === STEP 7 ===
# create table of alternative, reference and total allele numbers
all_intersects_oar = insert_ref_and_alt_allele_numbers(intersect_oar)
all_intersects_chi = insert_ref_and_alt_allele_numbers(intersect_chi)
# === ====== ===


# === STEP 8 === 
all_alt_freqs = find_alt_freqs(all_intersects_chi, all_intersects_oar)
chi_alt_freq_with_0_or_1 = all_alt_freqs.query('Chi_Alt_Freq == 0 | Chi_Alt_Freq == 1')

total_read_numbers = len(chi_alt_freq_with_0_or_1)
chi_read_numbers = len(chi_alt_freq_with_0_or_1.query('Chi_Alt_Freq == 0'))
oar_read_numbers = len(chi_alt_freq_with_0_or_1.query('Oar_Alt_Freq == 0'))

print("Total read number: ", total_read_numbers)
print("Chi read number:", chi_read_numbers)
print("Oar read number:", oar_read_numbers)

from scipy.stats import binom_test
binom_test_result = binom_test(oar_read_numbers, total_read_numbers, p=0.5, alternative='two-sided')

read_numbers_with_binom_result = pd.DataFrame({
                                        "Total_reads": [total_read_numbers],
                                        "Oar_reads": [oar_read_numbers],
                                        "Chi_reads": [chi_read_numbers],
                                        "p_value": [binom_test_result]
                                    })
print("\b\n>>> RESULTS <<<")
print(read_numbers_with_binom_result)
# === ====== ===
