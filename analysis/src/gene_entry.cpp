/*    
 *    gene_entry.cpp		
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


