/*    
 *    analysis.cpp		
 *
 *    Copyright (C) 2014 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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

