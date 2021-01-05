import pathlib

global DEBUG
DEBUG = False

# input file locations
# change file names according to the valid files to work with
current_sample_dir = "2nd_sample" # sample directory name to work with
current_dir = pathlib.Path(__file__).parent.absolute()

oar_ref_file = f"{current_dir}/data/Ovis_aries.Oar_v3.1.dna_rm.chromosome.MT.fa"
chi_ref_file = f"{current_dir}/data/Capra_hircus.ARS1.dna.chromosome.MT.fa"
transv_poly_oar_file = f"{current_dir}/data/tmp/trvposoar.bed"
transv_poly_chi_file = f"{current_dir}/data/tmp/trvposchi.bed"
aligned_sample_to_oar_file = f"{current_dir}/data/tmp/{current_sample_dir}/o_dp_tps083_MT.bam"
aligned_sample_to_chi_file = f"{current_dir}/data/tmp/{current_sample_dir}/c_dp_tps083_MT.bam"

# outputs file locations
data_out_dir = f"{current_dir}/data_out"
current_sample_data_out_dir = f"{data_out_dir}/{current_sample_dir}"

sorted_aligned_sample_to_oar_file = f"{current_sample_data_out_dir}/sorted_o_dp_MT.bam"
sorted_aligned_sample_to_chi_file = f"{current_sample_data_out_dir}/sorted_c_dp_MT.bam"
transv_sample_oar_file = f"{current_sample_data_out_dir}/transv_oar_samp.bed"
transv_sample_chi_file = f"{current_sample_data_out_dir}/transv_chi_samp.bed"
sample_oar_bamtobed_file = f"{current_sample_data_out_dir}/samp_oar_btb.bed"
sample_chi_bamtobed_file = f"{current_sample_data_out_dir}/samp_chi_btb.bed"
shared_oar_file = f"{current_sample_data_out_dir}/shared_oar.bed"
shared_chi_file = f"{current_sample_data_out_dir}/shared_chi.bed"
uniq_intersections_sample_oar_file = f"{current_sample_data_out_dir}/uniq_oar.bed"
uniq_intersections_sample_chi_file = f"{current_sample_data_out_dir}/uniq_chi.bed"