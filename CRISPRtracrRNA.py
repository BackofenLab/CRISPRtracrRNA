from modules.pipelines import CompleteTracrSearchWithModel
from modules.pipelines import TracrSearchWihtModelOnly
from modules.argument_parser import parse_input_arguments


def main():
    args = parse_input_arguments()
    run_type = args.run_type
    folder_input = args.input_folder
    folder_output = args.output_folder
    temp_folder_path = args.temp_folder_path
    output_file_name = args.output_summary_file
    if args.model_type == "II":
        model = ["tools/cov_matrix/s2_b14/cov_matrix.cm"]
    else:
        model = [
            "tools/cov_matrix/type_5/cluster1.top30.cm",
            "tools/cov_matrix/type_5/cluster2.top30.cm",
            "tools/cov_matrix/type_5/cluster3.top30.cm",
            "tools/cov_matrix/type_5/cluster4.top30.cm"
        ]
    anti_repeat_similarity = args.anti_repeat_similarity_threshold
    anti_repeat_coverage = args.anti_repeat_coverage_threshold
    dict_weights = {"crispr_array_score": args.weight_crispr_array_score,
                    "anti_repeat_similarity": args.weight_anti_repeat_sim,
                    "anti_repeat_coverage": args.weight_anti_repeat_coverage,
                    "anti_repeat_similarity_coverage_multiplication": args.weight_anti_sim_coverage,
                    "interaction_energy": args.weight_interaction_score * -1/10,
                    "tail_model_hit_score": args.weight_model_hit_score * 1/10,
                    "best_terminator_score": args.weight_terminator_hit_score,
                    "consistency_terminator_orientation": args.weight_consistency_orientation,
                    "consistency_anti_repeat_tail": args.weight_consistency_anti_repeat_tail,
                    "consistency_tail_terminator": args.weight_consistency_tail_terminator}
    perform_anti_repeat_search = args.perform_type_v_anti_repeat_analysis

    if run_type == "complete_run":
        complete_tracr_search = CompleteTracrSearchWithModel(folder_output=folder_output,
                                                             folder_genome=folder_input,
                                                             temp_folder_path=temp_folder_path,
                                                             output_file_name=output_file_name,
                                                             blast_similarity=anti_repeat_similarity,
                                                             blast_coverage=anti_repeat_coverage,
                                                             hmm_model=model,
                                                             dict_weights=dict_weights)
    else:
        model_tracr_search = TracrSearchWihtModelOnly(folder_input=folder_input,
                                                      folder_output=folder_output,
                                                      summary_file_name=output_file_name,
                                                      path_to_model=model,
                                                      flag_perform_anti_repeat_search=perform_anti_repeat_search)


if __name__ == "__main__":
    main()
