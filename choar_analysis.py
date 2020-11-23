import pysam
import pathlib
import pybedtools
from util import *

# file paths
current_dir = pathlib.Path(__file__).parent.absolute()

aligned_sample_to_oar_file = f"{current_dir}/data/tmp/1st_sample/o_dp_tps062_MT.bam"
aligned_sample_to_chi_file = f"{current_dir}/data/tmp/1st_sample/c_dp_tps062_MT.bam"
sorted_aligned_sample_to_oar_file = f"{current_dir}/data_out/1st_sample/sorted_o_dp_tps062_MT.bam"
sorted_aligned_sample_to_chi_file = f"{current_dir}/data_out/1st_sample/sorted_c_dp_tps062_MT.bam"
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
variants_info.to_csv(f"{current_dir}/data_out/1st_sample/transv_oar_samp.bed", 
                    header=False, index=False, sep="\t", mode="w")

# bedtools bamtobed -i bamfile > samp_oarbtb.bed 
sample_oar_in_bam = pybedtools.example_bedtool(sorted_aligned_sample_to_oar_file)
sample_oar_in_bed = sample_oar_in_bam.bam_to_bed()
write_to_file(sample_oar_in_bed, 
                            f"{current_dir}/data_out/1st_sample/samp_oar_btb.bed")   

# ----- For Goat ------
# sort and index the sample's aligned sequence to the goat's mtDNA reference
pysam.sort("-o", sorted_aligned_sample_to_chi_file, aligned_sample_to_chi_file)
pysam.index(sorted_aligned_sample_to_chi_file)

# bedtools bamtobed -i bamfile > samp_oarbtb.bed 
sample_chi_in_bam = pybedtools.example_bedtool(sorted_aligned_sample_to_chi_file)
sample_chi_in_bed = sample_chi_in_bam.bam_to_bed()
write_to_file(sample_chi_in_bed, 
                            f"{current_dir}/data_out/1st_sample/samp_chi_btb.bed")
# --------------------

# sharead.R
oar_file = f"{current_dir}/data_out/1st_sample/samp_oar_btb.bed"
chi_file = f"{current_dir}/data_out/1st_sample/samp_chi_btb.bed"
find_shared_reads(chi_file, oar_file)