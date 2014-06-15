#include "gene_entry.h"

test_cuff::test_cuff()
{
	fpkm_g1 = 0;
	fpkm_g2 = 0;
	//foldchange = "n/a";
	p_value = 1.0;
	significant = false;
}

isoform::isoform()
{

}


gene::gene()
{
	rangeLow = 0;
	rangeHigh = 0;
	is_annotated = false;

	cuff_anno_diffexpr = false;
	cuff_anno_difftrans = false;
	cuff_anno_difftrans_pvalue = 1.0;

	rangeLow_cuff = 0;
	rangeHigh_cuff = 0;
	cuff_no_anno_diffexpr = false;
	cuff_no_anno_difftrans = false;
	cuff_no_anno_difftrans_pvalue = 1.0;
	cuff_no_anno_diffsplicing = false;
	cuff_no_anno_diffsplicing_pvalue = 1.0;
	cuff_no_anno_diffsplicing_sqrtJS = 0;

	dfs_diffexpr = false;
	dfs_expr_stat = 0;
	dfs_fold_change = 0;
	dfs_cov1 = 0;
	dfs_cov2 = 0;
	dfs_difftrans = false;
	dfs_trans_stat = 0;
	dfs_sqrtJSD = 0;
	//dfs_as_category = "n/a";

	mts_diffexpr = false;
	mts_expr_stat = 0;
	mts_fold_change = 0;
	mts_cov1 = 0;
	mts_cov2 = 0;
	mts_difftrans = false;
	mts_trans_stat = 0;
	mts_sqrtJSD = 0;
}


vector <gene> global_list_gene_anno;
vector <gene> global_list_gene_novel;


