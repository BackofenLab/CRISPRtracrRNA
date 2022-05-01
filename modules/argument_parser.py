import argparse


def parse_input_arguments():
    parser = argparse.ArgumentParser(description='Run CRISPRtracrRNA ')
    parser.add_argument('--input_folder', type=str, default=None,
                        help='input folder with genomes (default: None)')
    parser.add_argument('--output_folder', type=str, default="CRISPRtracrRNA_output",
                        help='output folder (default: CRISPRtracrRNA_output)')
    parser.add_argument('--output_summary_file', type=str, default="CRISPRtracrRNA_result.csv",
                        help='output summary file (default: CRISPRtracrRNA_result.csv)')
    parser.add_argument('--temp_folder_path', type=str, default="temp",
                        help='temp folder path (default: temp)')
    parser.add_argument('--run_type', type=str, default="complete_run",
                        help='run type (default: complete_run )')
    parser.add_argument('--model_type', type=str, default="II",
                        help='model type (default: II)')
    parser.add_argument('--anti_repeat_similarity_threshold', type=float, default=0.7,
                        help='anti-repeat similarity threshold (default: 0.7)')
    parser.add_argument('--anti_repeat_coverage_threshold', type=float, default=0.6,
                        help='anti-repeat coverage threshold (default: 0.6)')
    parser.add_argument('--weight_crispr_array_score', type=float, default=0.5,
                        help='weight of crispr array score (default: 0.5)')
    parser.add_argument('--weight_anti_repeat_sim', type=float, default=0.5,
                        help='weight of anti-repeat similarity (default: 0.5)')
    parser.add_argument('--weight_anti_repeat_coverage', type=float, default=0.5,
                        help='weight of anti-repeat coverage (default: 0.5)')
    parser.add_argument('--weight_anti_sim_coverage', type=float, default=0.5,
                        help='weight of anti-similarity coverage (default: 0.5)')
    parser.add_argument('--weight_interaction_score', type=float, default=0.6,
                        help='weight of interaction score (default: 0.6)')
    parser.add_argument('--weight_model_hit_score', type=float, default=0.9,
                        help='weight of model hit score (default: 0.9)')
    parser.add_argument('--weight_terminator_hit_score', type=float, default=0.9,
                        help='weight of terminator hit score (default: 0.9)')
    parser.add_argument('--weight_consistency_orientation', type=float, default=0.1,
                        help='weight of consistency orientation (default: 0.1)')
    parser.add_argument('--weight_consistency_anti_repeat_tail', type=float, default=0.1,
                        help='weight of consistency anti-repeat tail (default: 0.1)')
    parser.add_argument('--weight_consistency_tail_terminator', type=float, default=0.1,
                        help='weight of consistency tail terminator (default: 0.1)')
    args = parser.parse_args()
    return args



