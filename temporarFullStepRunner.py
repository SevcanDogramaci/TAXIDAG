import pysam
import pybedtools
from util import *
from re import sub
from file_paths import *

# samtools faidx sequence1.fasta
pysam.faidx(oar_ref_file)
pysam.faidx(chi_ref_file)


# === STEP 1 ===
# sort and index the sample's aligned sequence to the sheep's mtDNA reference
pysam.sort("-o", sorted_aligned_sample_to_oar_file, aligned_sample_to_oar_file)
pysam.index(sorted_aligned_sample_to_oar_file)

# sort and index the sample's aligned sequence to the goat's mtDNA reference
pysam.sort("-o", sorted_aligned_sample_to_chi_file, aligned_sample_to_chi_file)
pysam.index(sorted_aligned_sample_to_chi_file)
# ==============


# === STEP 2 ===
# bcftools mpileup -B -f fastafile -R bedfile bamfile | 
# bcftools call -mV indels -A --ploidy 1 -o tr_samp_oar.vcf
print("Snp calling is starting for oar ...")
variant_calls_oar = pysam.mpileup("-f", oar_ref_file, "-B", "-l", transv_poly_oar_file, sorted_aligned_sample_to_oar_file)
print("Snp calling finished for oar...")

# bcftools mpileup -B -f fastafile -R bedfile bamfile | 
# bcftools call -mV indels -A --ploidy 1 -o tr_samp_oar.vcf
print("Snp calling is starting for chi ...")
variant_calls_chi = pysam.mpileup("-f", chi_ref_file, "-B", "-l", transv_poly_chi_file, sorted_aligned_sample_to_chi_file)
print("Snp calling finished for chi...")
# === ====== ===


# Open comment if you want to see variant calls in a txt file
# write_to_file(variant_calls_oar, f"{current_dir}/data_out/1st_sample/snpCalls.txt")


# === STEP 3 ===
variant_calls_oar = filter_variant_calls(variant_calls_oar)
variant_calls_chi = filter_variant_calls(variant_calls_chi)

# awk '!( $3 == "T" && $4 == "C" || $3 == "C" && $4 == "T" || 
#         $3 == "A" && $4 == "G" || $3 == "G" && $4 == "A" )
print("\n\n>>> FILTER POSTMORTEM TRANSITIONS FOR OAR<<<")
variant_calls_oar = filter_post_mortem_transitions(variant_calls_oar)
print(variant_calls_oar.head(20))
print("\n\n>>> FILTER POSTMORTEM TRANSITIONS FOR CHI<<<")
variant_calls_chi = filter_post_mortem_transitions(variant_calls_chi)
print(variant_calls_chi.head(20))

# awk '{print $1,$2-1,$2,$4}'
print("\n\n>>> CREATE BED FILE FOR OAR<<<")
variants_info_oar = get_variants_info(variant_calls_oar)
variants_info_oar.to_csv(transv_sample_oar_file, header=False, index=False, sep="\t", mode="w")
print("\n\n>>> CREATE BED FILE FOR CHI<<<")
variants_info_chi = get_variants_info(variant_calls_chi)
variants_info_chi.to_csv(transv_sample_chi_file, header=False, index=False, sep="\t", mode="w")
# === ====== ===


# === STEP 4 ===
# bedtools bamtobed -i bamfile > samp_oarbtb.bed 
sample_oar_in_bam = pybedtools.example_bedtool(sorted_aligned_sample_to_oar_file)
sample_oar_in_bed = sample_oar_in_bam.bam_to_bed()
write_to_file(sample_oar_in_bed, sample_oar_bamtobed_file)   

# bedtools bamtobed -i bamfile > samp_oarbtb.bed 
sample_chi_in_bam = pybedtools.example_bedtool(sorted_aligned_sample_to_chi_file)
sample_chi_in_bed = sample_chi_in_bam.bam_to_bed()
write_to_file(sample_chi_in_bed, sample_chi_bamtobed_file)
# === ====== ===


# === STEP 5 ===
# sharead.R
find_shared_reads(sample_chi_bamtobed_file, sample_oar_bamtobed_file)
# === ====== ===


# === STEP 6 ===
# bedtools intersect -a sha.oar.bed -b transvoar_samp.bed -wb 
shared_oar = pybedtools.example_bedtool(shared_oar_file)
transv_poly_oar = pybedtools.example_bedtool(transv_sample_oar_file)
intersect_oar = shared_oar.intersect(transv_poly_oar, wb=True)
intersect_oar = intersect_oar.to_dataframe()

shared_chi = pybedtools.example_bedtool(shared_chi_file)
transv_poly_chi = pybedtools.example_bedtool(transv_sample_chi_file)
intersect_chi = shared_chi.intersect(transv_poly_chi, wb=True)
intersect_chi = intersect_chi.to_dataframe()

intersect_oar = pd.concat([intersect_oar["name"], intersect_oar["blockCount"]], axis=1)
intersect_chi = pd.concat([intersect_chi["name"], intersect_chi["blockCount"]], axis=1)

# |awk '{gsub(/A|T|G|C/,"N",$2)}1' | awk '{gsub(/N,N/,"N",$2)}1'
intersect_oar["blockCount"] = intersect_oar["blockCount"].apply(
                                    lambda item: sub("A|a|T|t|G|g|C|c|(N,N)", "N", item))
intersect_chi["blockCount"] = intersect_chi["blockCount"].apply(
                                    lambda item: sub("A|a|T|t|G|g|C|c|(N,N)", "N", item))

# sort | uniq -c |awk '{ print $1,'\t',$2,'\t',$3 }'
intersect_oar = intersect_oar.value_counts().reset_index(name='counts')
intersect_oar = intersect_oar.sort_values(by=["name"])

intersect_chi = intersect_chi.value_counts().reset_index(name='counts')
intersect_chi = intersect_chi.sort_values(by=["name"])

print("\n\n >>> INTERSECTIONS FOR OAR <<<")
print(intersect_oar.head(10))
intersect_oar.to_csv(uniq_intersections_sample_oar_file, header=False, sep="\t", mode="w")
print("\n\n >>> INTERSECTIONS FOR CHI <<<")
print(intersect_chi.head(10))
intersect_chi.to_csv(uniq_intersections_sample_chi_file, header=False, sep="\t", mode="w")
# === ====== ===


# === STEP 7 ===
# create table of alternative, reference and total allele numbers
intersect_oar["id"] = range(1, len(intersect_oar)+1)
intersect_chi["id"] = range(1, len(intersect_chi)+1)

intersect_oar_with_ref_allele = intersect_oar[intersect_oar["blockCount"] == '.']
intersect_oar_with_alt_allele = intersect_oar[intersect_oar["blockCount"] == 'N']

intersect_chi_with_ref_allele = intersect_chi[intersect_chi["blockCount"] == '.']
intersect_chi_with_alt_allele = intersect_chi[intersect_chi["blockCount"] == 'N']

intersect_oar_with_ref_allele.insert(
    len(intersect_oar_with_ref_allele.columns),
    "Ref", intersect_oar_with_ref_allele["counts"]
)
intersect_oar_with_ref_allele.insert(
    len(intersect_oar_with_ref_allele.columns),
    "Alt", [0]*len(intersect_oar_with_ref_allele)
)

intersect_chi_with_ref_allele.insert(
    len(intersect_chi_with_ref_allele.columns),
    "Ref", intersect_chi_with_ref_allele["counts"]
)
intersect_chi_with_ref_allele.insert(
    len(intersect_chi_with_ref_allele.columns),
    "Alt", [0]*len(intersect_chi_with_ref_allele)
)

intersect_oar_with_alt_allele.insert(
    len(intersect_oar_with_alt_allele.columns),
    "Alt", intersect_oar_with_alt_allele["counts"]
)
intersect_oar_with_alt_allele.insert(
    len(intersect_oar_with_alt_allele.columns),
    "Ref", [0]*len(intersect_oar_with_alt_allele)
)

intersect_chi_with_alt_allele.insert(
    len(intersect_chi_with_alt_allele.columns),
    "Alt", intersect_chi_with_alt_allele["counts"]
)
intersect_chi_with_alt_allele.insert(
    len(intersect_chi_with_alt_allele.columns),
    "Ref", [0]*len(intersect_chi_with_alt_allele)
)

all_intersects_oar = pd.concat([intersect_oar_with_ref_allele, intersect_oar_with_alt_allele])
all_intersects_oar.insert(
    len(all_intersects_oar.columns),
    "Total", all_intersects_oar["Ref"] + all_intersects_oar["Alt"]
)
all_intersects_oar = filter_columns(all_intersects_oar, ["counts", "blockCount"])
all_intersects_oar = all_intersects_oar.sort_values(by=["id"])
print("\n\n>>> ALL INTERSECTIONS FOR OAR <<<")
print(all_intersects_oar.head(5))
all_intersects_oar.to_csv("intersections_oar.bed", header=True, index=False, sep="\t", mode="w")

all_intersects_chi = pd.concat([intersect_chi_with_ref_allele, intersect_chi_with_alt_allele])
all_intersects_chi.insert(
    len(all_intersects_chi.columns),
    "Total", all_intersects_chi["Ref"] + all_intersects_chi["Alt"]
)
all_intersects_chi = filter_columns(all_intersects_chi, ["counts", "blockCount"])
all_intersects_chi = all_intersects_chi.sort_values(by=["id"])
print("\n\n>>> ALL INTERSECTIONS FOR CHI <<<")
print(all_intersects_chi.head(5))
all_intersects_chi.to_csv("intersections_chi.bed", header=True, index=False, sep="\t", mode="w")
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
