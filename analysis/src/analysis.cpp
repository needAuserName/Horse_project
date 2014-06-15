#include "analysis.h"
#include "parse_result.h"
#include "write_output.h"


int main(int argc, char* argv[])
{
	parse_annotationGTF("data/Equus.gtf");

	parse_cuffdiff_expr("data/cuff_anno_gene_exp.diff", "error/");
	parse_cuffdiff_trans("data/cuff_anno_isoform_exp.diff", "error/");

	parse_cuffdiff_expr_novel("data/cuff_no_anno_gene_exp.diff", "error/");
	parse_cuffdiff_trans_novel("data/cuff_no_anno_isoform_exp.diff", "error/");
	parse_cuffdiff_splicing_novel("data/cuff_no_anno_splicing.diff", "error/");

	parse_dfs_expr("data/dfs_differential_expression.txt", "error/");
	parse_dfs_trans("data/dfs_differential_transcription.txt", "error/");

	parse_mts_expr("data/mts_differential_expression.txt", "error/");
	parse_mts_trans("data/mts_differential_transcription.txt", "error/");

	write_output_anno("result/");
	write_output_novel("result/");

	return 0;
}

