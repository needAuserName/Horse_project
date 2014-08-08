/*    
 *    gene_entry.h		
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

/************************************************************************/
/*						gene_entry.h                                    */
/* Declaration of entries used in the analysis of results.				*/
/************************************************************************/

#ifndef GENE_ENTRY
#define GENE_ENTRY

#include "analysis.h"

// make the fields as string because some can be inf
class test_cuff
{
public:
	double fpkm_g1; // FPKM of group 1
	double fpkm_g2; // FPKM of group 2
	string foldchange; // log2(fold_change), fold_change=fpkm_g2/fpkm_g1
	double p_value; // the q_value in Cuffdiff results, which is the corrected p_value with correction of multiple testing
	bool significant; // significant difference or not

	test_cuff();
};

class isoform
{
public:
	// basic information
	string name; // the annotated name or the Cufflinks name
	string id; // Ensembl ID

	// optional information: only available for Cuffdiff isoforms or Multisplice results, etc
	test_cuff cuff;

	isoform();
};


class gene
{
public:
	// basic information
	string geneNm; // the annotated name if this gene is in annotation, the cufflinks name otherwise
	string geneID; // Ensembl ID
	string chrNm; // chromosome
	long rangeLow; // start position in annotation or in cufflinks, depending on is_annotated
	long rangeHigh; // end position in annotation or in cufflinks, depending on is_annotated
	bool is_annotated; // in annotation or not

	// annotation info
	vector <isoform> anno_isoforms; // list of isoforms in annotation, only names will be sufficient

	// Cuffdiff with annotation info
	bool		cuff_anno_diffexpr; // differential gene expression, same as significant in gene_exp.diff
	test_cuff	cuff_anno_expr_stat; // information extracted from gene_exp.diff
	bool		cuff_anno_difftrans; // differential transcription, yes if at least one of the isoforms of this gene is significant in isoform_exp.diff
	double		cuff_anno_difftrans_pvalue; // minimum p-value of all isoforms in the gene, used to represent the p-value of the gene
	vector <isoform> cuff_anno_isoforms; // list of isoforms and their stats from isoform_exp.diff

	// Cuffdiff without annotation info
	string geneNm_cuff; // gene name from cufflinks no annotation, useful is is_annotation==true, get from isoform_exp.diff
	long rangeLow_cuff; // get from isoform_exp.diff
	long rangeHigh_cuff; // get from isoform_exp.diff
	bool		cuff_no_anno_diffexpr; // differential gene expression, same as significant in gene_exp.diff
	test_cuff	cuff_no_anno_expr_stat; // information extracted from gene_exp.diff
	bool		cuff_no_anno_difftrans; // differential transcription, yes if at least one of the isoforms of this gene is significant in isoform_exp.diff
	double		cuff_no_anno_difftrans_pvalue; // minimum p-value of all isoforms in the gene, used to represent the p-value of the gene
	vector <isoform> cuff_no_anno_isoforms; // list of isoforms and their stats from isoform_exp.diff
	bool		cuff_no_anno_diffsplicing; // differential splicing
	double		cuff_no_anno_diffsplicing_pvalue;
	double		cuff_no_anno_diffsplicing_sqrtJS;

	// diffsplice info
	bool		dfs_diffexpr; // differential gene expression of diffsplice, from differential_expression.txt
	double		dfs_fold_change; // fold_change, you can talk log2 in order to match cufflinks results
	double		dfs_expr_stat; // stat_diff_expr(|stat_d_expected-stat_d|)
	double		dfs_cov1; // diffsplice gene coverage in group1
	double		dfs_cov2; 
	bool		dfs_difftrans; // differential transcription of diffsplice, from differential_transcription.txt
	double		dfs_trans_stat; 
	string		dfs_as_category; // category of alternative splicing
	double		dfs_sqrtJSD; // square root of JSD

	// multisplice info
	bool		mts_diffexpr; // differential gene expression of diffsplice, from differential_expression.txt
	double		mts_fold_change; // fold_change, you can talk log2 in order to match cufflinks results
	double		mts_expr_stat;
	double		mts_cov1; // diffsplice gene coverage in group1
	double		mts_cov2; 
	bool		mts_difftrans; // differential transcription of diffsplice, from differential_transcription.txt
	double		mts_trans_stat;
	double		mts_sqrtJSD; // square root of JSD

	// astroid info

	gene();
};

extern vector <gene> global_list_gene_anno;
extern vector <gene> global_list_gene_novel;

#endif
