import pysam
import pybedtools
from util import *
from re import sub
from file_paths import *

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
# write_to_file(variant_calls, f"{current_dir}/data_out/1st_sample/snpCalls.txt")

variant_calls = filter_variant_calls(variant_calls)

# awk '!( $3 == "T" && $4 == "C" || $3 == "C" && $4 == "T" || 
#         $3 == "A" && $4 == "G" || $3 == "G" && $4 == "A" )
print("\n\n>>> FILTER POSTMORTEM TRANSITIONS <<<")
variant_calls = filter_post_mortem_transitions(variant_calls)
print(variant_calls.head(20))

# awk '{print $1,$2-1,$2,$4}'
print("\n\n>>> CREATE BED FILE <<<")
variants_info = get_variants_info(variant_calls)
variants_info.to_csv(transv_sample_oar_file, header=False, index=False, sep="\t", mode="w")

# bedtools bamtobed -i bamfile > samp_oarbtb.bed 
sample_oar_in_bam = pybedtools.example_bedtool(sorted_aligned_sample_to_oar_file)
sample_oar_in_bed = sample_oar_in_bam.bam_to_bed()
write_to_file(sample_oar_in_bed, sample_oar_bamtobed_file)   

# ----- For Goat ------
# sort and index the sample's aligned sequence to the goat's mtDNA reference
pysam.sort("-o", sorted_aligned_sample_to_chi_file, aligned_sample_to_chi_file)
pysam.index(sorted_aligned_sample_to_chi_file)

# bedtools bamtobed -i bamfile > samp_oarbtb.bed 
sample_chi_in_bam = pybedtools.example_bedtool(sorted_aligned_sample_to_chi_file)
sample_chi_in_bed = sample_chi_in_bam.bam_to_bed()
write_to_file(sample_chi_in_bed, sample_chi_bamtobed_file)
# --------------------

# sharead.R
find_shared_reads(sample_chi_bamtobed_file, sample_oar_bamtobed_file)

# bedtools intersect -a sha.oar.bed -b transvoar_samp.bed -wb 
shared_oar = pybedtools.example_bedtool(shared_oar_file)
transv_poly_oar = pybedtools.example_bedtool(transv_sample_oar_file)
intersect_oar = shared_oar.intersect(transv_poly_oar, wb=True)
intersect_oar = intersect_oar.to_dataframe()

intersect_oar = pd.concat([intersect_oar["name"], intersect_oar["blockCount"]], axis=1)

# |awk '{gsub(/A|T|G|C/,"N",$2)}1' | awk '{gsub(/N,N/,"N",$2)}1'
intersect_oar["blockCount"] = intersect_oar["blockCount"].apply(
                                    lambda item: sub("A|T|G|C|(N,N)", "N", item))

# sort | uniq -c |awk '{ print $1,'\t',$2,'\t',$3 }'
intersect_oar = intersect_oar.value_counts().to_frame().sort_values(by=["name"])
print("\n\n >>> INTERSECTIONS <<<")
print(intersect_oar.head(10))
intersect_oar.to_csv(uniq_intersections_sample_oar_file, header=False, sep="\t", mode="w")
