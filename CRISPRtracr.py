from modules.pipelines import CompleteTracrSearch
from modules.pipelines import CompleteTracrSearchWithModel

from run_inputs import folder_genomes, file_crispr_idetify, folder_temp, folder_output, output_file_name, hmm_model


def main():
    complete_tracr_search = CompleteTracrSearchWithModel(folder_output=folder_output,
                                                         folder_genome=folder_genomes,
                                                         crispr_identify_summary_file=file_crispr_idetify,
                                                         temp_folder_path=folder_temp,
                                                         output_file_name=output_file_name,
                                                         blast_similarity=0.7,
                                                         blast_coverage=0.6,
                                                         hmm_model=hmm_model)


if __name__ == "__main__":
    main()
