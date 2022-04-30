from modules.pipelines import CompleteTracrSearchWithModel
from modules.argument_parser import parse_input_arguments


def main():
    args = parse_input_arguments()
    run_type = args.run_type
    folder_input = args.input_folder
    folder_output = args.output_folder
    temp_folder_path = args.temp_folder_path
    output_file_name = args.output_summary_file
    model = "tools/cov_matrix/s2_b14/cov_matrix.cm" if args.model == "II" else "tools/cov_matrix/type_5/cov_matrix.cm"
    anti_repeat_similarity = args.anti_repeat_similarity
    anti_repeat_coverage = args.anti_repeat_coverage
    dict_weights = {"crispr_array_score": args.crispr_array_score,
                    "anti_repeat_similarity": args.weight_anti_repeat_sim,
                    "anti_repeat_coverage": args.weight_anti_repeat_coverage,
                    "anti_repeat_similarity_coverage_multiplication": args.weight_anti_sim_coverage,
                    "interaction_energy": args.weight_interaction_score * -1/10,
                    "tail_model_hit_score": args.weight_model_hit_score * 1/10,
                    "weight_terminator_hit_score": args.weight_terminator_hit_score,
                    "consistency_terminator_orientation": args.weight_consistency_orientation,
                    "consistency_anti_repeat_tail": args.weight_consistency_anti_repeat_tail,
                    "consistency_tail_terminator": args.weight_consistency_tail_terminator}

    if run_type == "complete_run":
        complete_tracr_search = CompleteTracrSearchWithModel(folder_output=folder_output,
                                                             folder_genome=folder_input,
                                                             temp_folder_path=temp_folder_path,
                                                             output_file_name=output_file_name,
                                                             blast_similarity=anti_repeat_similarity,
                                                             blast_coverage=anti_repeat_coverage,
                                                             hmm_model=model,
                                                             dict_weights=dict_weights)


if __name__ == "__main__":
    main()
