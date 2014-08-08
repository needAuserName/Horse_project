/*    
 *    write_output.cpp		
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

#include "write_output.h"

string bool2str(bool flag)
{
	return flag? "yes":"no";
}

double Log2(double n)  
{  
	return log(n) / log(2.);  
}

int write_output_anno(string dir_result)
{
	string filename = dir_result + "/summary_diff_gene_expression_anno.txt";
	ofstream outfile_expr(filename.c_str());
	filename = dir_result + "/summary_diff_transcription_anno.txt";
	ofstream outfile_trans(filename.c_str());
	
	if (!outfile_expr.is_open() || !outfile_trans.is_open())
	{
		cout << "Error: cannot open the output summary file." << endl;
		return 1;
	}

	outfile_expr << "Gene ID\tChromosome\tPosition\tAnnotated\t# of annotated isoforms\t"
		<< "Cuffdiff annotation\tCuffdiff annotation p-value\tCuffdiff annotation log2(fold change)\tCuffdiff annotation FPKM group1\tCuffdiff annotation FPKM group2\t"
		<< "Cuffdiff w/o annotation\tCuffdiff w/o annotation p-value\tCuffdiff w/o annotation log2(fold change)\tCuffdiff w/o annotation FPKM group1\tCuffdiff w/o annotation FPKM group2\t"
		<< "DiffSplice\tDiffSplice test statistic\tDiffSplice log2(fold change)\tDiffSplice coverage group1\tDiffSplice coverage group2\t"
		<< "MultiSplice\tMultiSplice test statistic\tMultiSplice log2(fold change)\tMultiSplice coverage group1\tMultiSplice coverage group2"
		<< endl;
	outfile_trans << "Gene ID\tChromosome\tPosition\tAnnotated\t# of annotated isoforms\t"
		<< "Cuffdiff annotation\tCuffdiff annotation p-value\tCuffdiff annotation changed isoforms and log2(fold change)\t"
		<< "Cuffdiff w/o annotation isoform\tCuffdiff w/o annotation isoform p-value\tCuffdiff w/o annotation changed isoforms and log2(fold change)\t"
		<< "Cuffdiff w/o annotation splicing\tCuffdiff w/o annotation splicing p-value\tCuffdiff w/o annotation splicing sqrtJS\t"
		<< "DiffSplice\tDiffSplice sqrt of JSD\tDiffSplice test statistic\tDiffSplice alter splicing category\t"
		<< "MultiSplice\tMultiSplice sqrt of JSD\tMultiSplice test statistic"
		<< endl;

	for (vector<gene>::iterator it_gene = global_list_gene_anno.begin() ; it_gene != global_list_gene_anno.end(); ++it_gene)
	{
		// differential expression
		outfile_expr << it_gene->geneID << "\t" << it_gene->chrNm << "\t" << it_gene->rangeLow << "-" << it_gene->rangeHigh << "\t" << bool2str(it_gene->is_annotated) << "\t" << (it_gene->anno_isoforms).size() << "\t";
		outfile_expr << bool2str(it_gene->cuff_anno_diffexpr) << "\t" << (it_gene->cuff_anno_expr_stat).p_value << "\t" << (it_gene->cuff_anno_expr_stat).foldchange << "\t" << (it_gene->cuff_anno_expr_stat).fpkm_g1 << "\t" << (it_gene->cuff_anno_expr_stat).fpkm_g2 << "\t";
		outfile_expr << bool2str(it_gene->cuff_no_anno_diffexpr) << "\t" << (it_gene->cuff_no_anno_expr_stat).p_value << "\t" << (it_gene->cuff_no_anno_expr_stat).foldchange << "\t" << (it_gene->cuff_no_anno_expr_stat).fpkm_g1 << "\t" << (it_gene->cuff_no_anno_expr_stat).fpkm_g2 << "\t";
		outfile_expr << bool2str(it_gene->dfs_diffexpr) << "\t" << it_gene->dfs_expr_stat << "\t";
		if (it_gene->dfs_fold_change > 0)
			outfile_expr << Log2(it_gene->dfs_fold_change);
		outfile_expr << "\t" << it_gene->dfs_cov1 << "\t" << it_gene->dfs_cov2;
		outfile_expr << "\t" << bool2str(it_gene->mts_diffexpr) << "\t" << it_gene->mts_expr_stat << "\t";
		if (it_gene->mts_fold_change > 0)
			outfile_expr << Log2(it_gene->mts_fold_change);
		outfile_expr << "\t" << it_gene->mts_cov1 << "\t" << it_gene->mts_cov2;
		outfile_expr << endl;

		// differential transcription
		outfile_trans << it_gene->geneID << "\t" << it_gene->chrNm << "\t" << it_gene->rangeLow << "-" << it_gene->rangeHigh << "\t" << bool2str(it_gene->is_annotated) << "\t" << (it_gene->anno_isoforms).size() << "\t";
		outfile_trans << bool2str(it_gene->cuff_anno_difftrans) << "\t" << it_gene->cuff_anno_difftrans_pvalue << "\t";
		for (vector<isoform>::iterator it_isoform = it_gene->cuff_anno_isoforms.begin(); it_isoform != it_gene->cuff_anno_isoforms.end(); ++it_isoform)
			if ((it_isoform->cuff).significant)
				outfile_trans << it_isoform->id << "(" << (it_isoform->cuff).foldchange << "), ";
		outfile_trans << "\t" << bool2str(it_gene->cuff_no_anno_difftrans) << "\t" << it_gene->cuff_no_anno_difftrans_pvalue << "\t";
		for (vector<isoform>::iterator it_isoform = it_gene->cuff_no_anno_isoforms.begin(); it_isoform != it_gene->cuff_no_anno_isoforms.end(); ++it_isoform)
			if ((it_isoform->cuff).significant)
				outfile_trans << it_isoform->id << "(" << (it_isoform->cuff).foldchange << "), ";
		outfile_trans << "\t" << bool2str(it_gene->cuff_no_anno_diffsplicing) << "\t" << it_gene->cuff_no_anno_diffsplicing_pvalue << "\t" << it_gene->cuff_no_anno_diffsplicing_sqrtJS;
		outfile_trans << "\t" << bool2str(it_gene->dfs_difftrans) << "\t" << it_gene->dfs_sqrtJSD << "\t" << it_gene->dfs_trans_stat << "\t" << it_gene->dfs_as_category;
		outfile_trans << "\t" << bool2str(it_gene->mts_difftrans) << "\t" << it_gene->mts_sqrtJSD << "\t" << it_gene->mts_trans_stat;
		outfile_trans << endl;
	}	

	outfile_expr.close();
	outfile_trans.close();

	return 0;
}


int write_output_novel(string dir_result)
{
	string filename = dir_result + "/summary_diff_gene_expression_novel.txt";
	ofstream outfile_expr(filename.c_str());
	filename = dir_result + "/summary_diff_transcription_novel.txt";
	ofstream outfile_trans(filename.c_str());

	if (!outfile_expr.is_open() || !outfile_trans.is_open())
	{
		cout << "Error: cannot open the output summary file." << endl;
		return 1;
	}

	outfile_expr << "Gene ID (Name)\tChromosome\tPosition\tAnnotated\t"
		<< "Cuffdiff w/o annotation\tCuffdiff w/o annotation p-value\tCuffdiff w/o annotation log2(fold change)\tCuffdiff w/o annotation FPKM group1\tCuffdiff w/o annotation FPKM group2\t"
		<< "DiffSplice\tDiffSplice test statistics\tDiffSplice log2(fold change)\tDiffSplice coverage group1\tDiffSplice coverage group2"
		<< endl;
	outfile_trans << "Gene ID (Name)\tChromosome\tPosition\tAnnotated\t"
		<< "Cuffdiff w/o annotation isoform\tCuffdiff w/o annotation isoform p-value\tCuffdiff w/o annotation changed isoforms and log2(fold change)\t"
		<< "Cuffdiff w/o annotation splicing\tCuffdiff w/o annotation splicing p-value\tCuffdiff w/o annotation splicing sqrtJS\t"
		<< "DiffSplice\tDiffSplice sqrt of JSD\tDiffSplice test statistics\tDiffSplice alter splicing category"
		<< endl;

	for (vector<gene>::iterator it_gene = global_list_gene_novel.begin() ; it_gene != global_list_gene_novel.end(); ++it_gene)
	{
		// differential expression
		outfile_expr << it_gene->geneID << " (" << it_gene->geneNm << ")\t" << it_gene->chrNm << "\t" << it_gene->rangeLow << "-" << it_gene->rangeHigh << "\t" << bool2str(it_gene->is_annotated) << "\t";
		outfile_expr << bool2str(it_gene->cuff_no_anno_diffexpr) << "\t" << (it_gene->cuff_no_anno_expr_stat).p_value << "\t" << (it_gene->cuff_no_anno_expr_stat).foldchange << "\t" << (it_gene->cuff_no_anno_expr_stat).fpkm_g1 << "\t" << (it_gene->cuff_no_anno_expr_stat).fpkm_g2 << "\t";
		outfile_expr << bool2str(it_gene->dfs_diffexpr) << "\t" << it_gene->dfs_expr_stat << "\t";
		if (it_gene->dfs_fold_change > 0)
			outfile_expr << Log2(it_gene->dfs_fold_change);
		outfile_expr << "\t" << it_gene->dfs_cov1 << "\t" << it_gene->dfs_cov2;
		outfile_expr << endl;

		// differential transcription
		outfile_trans << it_gene->geneID << " (" << it_gene->geneNm << ")\t" << it_gene->chrNm << "\t" << it_gene->rangeLow << "-" << it_gene->rangeHigh << "\t" << bool2str(it_gene->is_annotated) << "\t";
		outfile_trans << bool2str(it_gene->cuff_no_anno_difftrans) << "\t" << it_gene->cuff_no_anno_difftrans_pvalue << "\t";
		for (vector<isoform>::iterator it_isoform = it_gene->cuff_no_anno_isoforms.begin(); it_isoform != it_gene->cuff_no_anno_isoforms.end(); ++it_isoform)
			if ((it_isoform->cuff).significant)
				outfile_trans << it_isoform->id << "(" << (it_isoform->cuff).foldchange << "), ";
		outfile_trans << "\t" << bool2str(it_gene->cuff_no_anno_diffsplicing) << "\t" << it_gene->cuff_no_anno_diffsplicing_pvalue << "\t" << it_gene->cuff_no_anno_diffsplicing_sqrtJS;
		outfile_trans << "\t" << bool2str(it_gene->dfs_difftrans) << "\t" << it_gene->dfs_sqrtJSD << "\t" << it_gene->dfs_trans_stat << "\t" << it_gene->dfs_as_category;
		outfile_trans << endl;
	}	

	outfile_expr.close();
	outfile_trans.close();

	return 0;
}
