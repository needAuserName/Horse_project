#include "parse_result.h"

bool str2bool(string str_logic)
{
	return str_logic.compare("yes") == 0 ? true : false;
}

int parse_annotationGTF(string filename)
{
	ifstream infile(filename.c_str());

	if (!infile.is_open())
	{
		cout << "Error: cannot open annotation GTF file." << endl;
		return 1;
	}

	string info_chr, info_feature, info_other;
	stringstream iss;
	char opt_field[500], opt_value[500];

	while (infile >> info_chr)
	{
		infile >> info_other;
		infile >> info_feature;
		if (info_feature.compare("gene") == 0)
		{
			gene new_gene;
			new_gene.chrNm = info_chr;
			new_gene.is_annotated = true;
			infile >> new_gene.rangeLow;
			infile >> new_gene.rangeHigh;
			for (int tmpcnt = 0; tmpcnt < 4; ++tmpcnt)
				getline(infile, info_other, '\t');
			getline(infile, info_other);

			if (!info_other.empty())
			{
				iss << info_other;
				while (iss.getline(opt_field, 500, ' '))
				{
					iss.getline(opt_value, 500, ' ');

					if (strcmp(opt_field, "gene_id") == 0)
					{
						new_gene.geneID = opt_value;
						new_gene.geneID = new_gene.geneID.substr(1, new_gene.geneID.size()-3);
					}
					else if (strcmp(opt_field, "gene_name") == 0)
					{
						new_gene.geneNm = opt_value;
						new_gene.geneNm = new_gene.geneNm.substr(1, new_gene.geneNm.size()-3);
					}
				}
				iss.clear();
				iss.str("");
			}
			else
				cout << "Warning: a gene entry with no info" << endl;
			
			global_list_gene_anno.push_back(new_gene);
		}
		else if (info_feature.compare("transcript") == 0)
		{
			isoform new_isoform;
			long tmp_rangeLow, tmp_rangeHigh;
			string tmp_geneID;
			infile >> tmp_rangeLow;
			infile >> tmp_rangeHigh;
			for (int tmpcnt = 0; tmpcnt < 4; ++tmpcnt)
				getline(infile, info_other, '\t');
			getline(infile, info_other);

			if (!info_other.empty())
			{
				iss << info_other;
				while (iss.getline(opt_field, 500, ' '))
				{
					iss.getline(opt_value, 500, ' ');

					if (strcmp(opt_field, "gene_id") == 0)
					{
						tmp_geneID = opt_value;
						tmp_geneID = tmp_geneID.substr(1, tmp_geneID.size()-3);
					}
					else if (strcmp(opt_field, "transcript_id") == 0)
					{
						new_isoform.id = opt_value;
						new_isoform.id = new_isoform.id.substr(1, new_isoform.id.size()-3);
					}
					else if (strcmp(opt_field, "transcript_name") == 0)
					{
						new_isoform.name = opt_value;
						new_isoform.name = new_isoform.name.substr(1, new_isoform.name.size()-3);
					}
				}
				iss.clear();
				iss.str("");

				vector<gene>::reverse_iterator it_gene;
				for (it_gene = global_list_gene_anno.rbegin(); it_gene != global_list_gene_anno.rend(); ++it_gene)
				{
					if ((it_gene->geneID).compare(tmp_geneID) == 0)
					{
						break;
					}
				}

				if (it_gene != global_list_gene_anno.rend())
				{
					// normal case, find the gene
					(it_gene->anno_isoforms).push_back(new_isoform);
				} 
				else
				{
					// abnormal case, cannot find the gene
					cout << "Warning: transcript " << new_isoform.id << " not belonging to any gene." << endl;
					gene new_gene;
					new_gene.chrNm = info_chr;
					new_gene.is_annotated = true;
					new_gene.rangeLow = tmp_rangeLow;
					new_gene.rangeHigh = tmp_rangeHigh;
					new_gene.geneNm = new_isoform.name;
					new_gene.geneID = new_isoform.id;
					new_gene.anno_isoforms.push_back(new_isoform);

					global_list_gene_anno.push_back(new_gene);
				}
				
			}
			else
				cout << "Warning: a transcript entry with no info" << endl;
		}
		else
		{
			getline(infile, info_other);
		}
	}

	infile.close();

	cout << "Recorded " << global_list_gene_anno.size() << " genes from annotation GTF." << endl;

	return 0;
}

// find a gene in the global_list_gene_anno that has the maximum overlap ratio (first the ratio of the incoming region, if the same then the ratio of the target region)
long locate_gene_by_overlap(const vector <gene> &list_gene, string chr, long rangeLow, long rangeHigh)
{
	long best_match = -1;
	double best_overlap_new = 0, best_overlap_target = 0;

	if (rangeHigh <= rangeLow)
		return best_match;

	for (long iLoop = 0; iLoop < list_gene.size(); ++iLoop)
	{
		if (list_gene[iLoop].chrNm.compare(chr) == 0 && list_gene[iLoop].rangeHigh > list_gene[iLoop].rangeLow)
		{
			long overlap = min(list_gene[iLoop].rangeHigh, rangeHigh) - max(list_gene[iLoop].rangeLow, rangeLow);
			double overlap_new = double(overlap) / (rangeHigh - rangeLow);
			double overlap_target = double(overlap) / (list_gene[iLoop].rangeHigh - list_gene[iLoop].rangeLow);
			if (overlap > 0 && (overlap_new > best_overlap_new || (overlap_new == best_overlap_new && overlap_target > best_overlap_target)))
			{
				best_match = iLoop;
				best_overlap_new = overlap_new;
				best_overlap_target = overlap_target;
			}
		}
	}

	return best_match;
}

int parse_dfs_expr(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open DiffSplice expression result." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/dfs_expr_discarded.txt").c_str() );

	string info_chr, info_other, titleline;
	long rangeLow, rangeHigh;
	
	getline(infile, titleline);

	while (infile >> info_chr)
	{
		infile >> rangeLow;
		infile >> rangeHigh;
		
		long index_match = locate_gene_by_overlap(global_list_gene_anno, info_chr, rangeLow, rangeHigh);
		if (index_match >= 0)
		{
			// find a match in annotated genes
			infile >> global_list_gene_anno.at(index_match).dfs_expr_stat; // stat_diff_expr
			infile >> global_list_gene_anno.at(index_match).dfs_fold_change;
			infile >> global_list_gene_anno.at(index_match).dfs_cov1; // coverage_group1
			infile >> global_list_gene_anno.at(index_match).dfs_cov2; // coverage_group2
			infile >> info_other; global_list_gene_anno.at(index_match).dfs_diffexpr = str2bool(info_other);
			getline(infile, info_other);
		} 
		else
		{
			// no match in annotated list
			long index_match = locate_gene_by_overlap(global_list_gene_novel, info_chr, rangeLow, rangeHigh);
			if (index_match >= 0)
			{
				// find a match in novel genes
				infile >> global_list_gene_novel.at(index_match).dfs_expr_stat; // stat_diff_expr
				infile >> global_list_gene_novel.at(index_match).dfs_fold_change;
				infile >> global_list_gene_novel.at(index_match).dfs_cov1; // coverage_group1
				infile >> global_list_gene_novel.at(index_match).dfs_cov2; // coverage_group2
				infile >> info_other; global_list_gene_novel.at(index_match).dfs_diffexpr = str2bool(info_other);
				getline(infile, info_other);
			} 
			else
			{
				getline(infile, info_other);
				outfile_error << info_chr << "\t" << rangeLow << "\t" << rangeHigh << info_other << endl;
			}
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}

int parse_dfs_trans(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open DiffSplice transcription result." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/dfs_trans_discarded.txt").c_str() );

	string info_chr, info_other, titleline;
	long rangeLow, rangeHigh;

	getline(infile, titleline);

	while (infile >> info_chr)
	{
		infile >> rangeLow;
		infile >> rangeHigh;

		long index_match_anno = locate_gene_by_overlap(global_list_gene_anno, info_chr, rangeLow, rangeHigh);
		if (index_match_anno >= 0)
		{
			// find a match in annotated genes
			infile >> global_list_gene_anno.at(index_match_anno).dfs_as_category;
			infile >> global_list_gene_anno.at(index_match_anno).dfs_trans_stat; // stat_diff_expr
			infile >> global_list_gene_anno.at(index_match_anno).dfs_sqrtJSD;
			infile >> info_other; // coverage_group1
			infile >> info_other; // coverage_group2
			infile >> info_other; global_list_gene_anno.at(index_match_anno).dfs_difftrans = str2bool(info_other);
			getline(infile, info_other);
		} 
		else
		{
			// no match in annotated list
			long index_match_novel = locate_gene_by_overlap(global_list_gene_novel, info_chr, rangeLow, rangeHigh);
			if (index_match_novel >= 0)
			{
				// find a match in novel genes
				infile >> global_list_gene_novel.at(index_match_novel).dfs_as_category;
				infile >> global_list_gene_novel.at(index_match_novel).dfs_trans_stat; // stat_diff_expr
				infile >> global_list_gene_novel.at(index_match_novel).dfs_sqrtJSD;
				infile >> info_other; // coverage_group1
				infile >> info_other; // coverage_group2
				infile >> info_other; global_list_gene_novel.at(index_match_novel).dfs_difftrans = str2bool(info_other);
				getline(infile, info_other);
			} 
			else
			{
				getline(infile, info_other);
				outfile_error << info_chr << "\t" << rangeLow << "\t" << rangeHigh << info_other << endl;
			}
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}

// find a gene in the global_list_gene_anno that has the maximum overlap ratio (first the ratio of the incoming region, if the same then the ratio of the target region)
long locate_gene_by_ID(const vector <gene> &list_gene, string geneID)
{
	for (long iLoop = 0; iLoop < list_gene.size(); ++iLoop)
	{
		if (list_gene[iLoop].geneID.compare(geneID) == 0)
		{
			return iLoop;
		}
	}

	return -1;
}

int parse_cuffdiff_expr(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open Cuffdiff gene expression result." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/cuffdiff_gene_discarded.txt").c_str() );

	string info_geneID, info_testID, info_other, titleline;
	
	getline(infile, titleline);

	while (infile >> info_testID)
	{
		infile >> info_geneID;

		long index_match = locate_gene_by_ID(global_list_gene_anno, info_geneID);
		if (index_match >= 0)
		{
			test_cuff new_cuff_res;

			// find a match in annotated genes
			infile >> info_other; // gene name
			infile >> info_other; // locus
			infile >> info_other; // sample_1
			infile >> info_other; // sample_2
			infile >> info_other; // status

			infile >> new_cuff_res.fpkm_g1;
			infile >> new_cuff_res.fpkm_g2;
			infile >> new_cuff_res.foldchange;
			if (new_cuff_res.foldchange.compare("inf") == 0 || new_cuff_res.foldchange.compare("-inf") == 0)
				new_cuff_res.foldchange = "";

			infile >> info_other; // test_stat
			infile >> info_other; // p_value before correction

			infile >> new_cuff_res.p_value;
			infile >> info_other; new_cuff_res.significant = str2bool(info_other);
			getline(infile, info_other);

			global_list_gene_anno.at(index_match).cuff_anno_diffexpr = new_cuff_res.significant;
			global_list_gene_anno.at(index_match).cuff_anno_expr_stat = new_cuff_res;
		} 
		else
		{
			// no match in annotated list
			getline(infile, info_other);
			outfile_error << info_testID << "\t" << info_geneID << info_other << endl;
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}


int parse_cuffdiff_trans(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open Cuffdiff isoform expression result." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/cuffdiff_isoform_discarded.txt").c_str() );

	string info_geneID, info_transID, info_other, titleline;

	getline(infile, titleline);

	while (infile >> info_transID)
	{
		infile >> info_geneID;

		long index_match = locate_gene_by_ID(global_list_gene_anno, info_geneID);
		if (index_match >= 0)
		{
			test_cuff new_cuff_res;

			// find a match in annotated genes
			infile >> info_other; // gene name
			infile >> info_other; // locus
			infile >> info_other; // sample_1
			infile >> info_other; // sample_2
			infile >> info_other; // status

			infile >> new_cuff_res.fpkm_g1;
			infile >> new_cuff_res.fpkm_g2;
			infile >> new_cuff_res.foldchange;
			if (new_cuff_res.foldchange.compare("inf") == 0 || new_cuff_res.foldchange.compare("-inf") == 0)
				new_cuff_res.foldchange = "";

			infile >> info_other; // test_stat
			infile >> info_other; // p_value before correction

			infile >> new_cuff_res.p_value;
			infile >> info_other; new_cuff_res.significant = str2bool(info_other);
			getline(infile, info_other);

			isoform new_isoform;
			new_isoform.id = info_transID;
			new_isoform.cuff = new_cuff_res;

			if (!global_list_gene_anno.at(index_match).cuff_anno_difftrans)
				global_list_gene_anno.at(index_match).cuff_anno_difftrans = new_cuff_res.significant;
			if (new_cuff_res.p_value < global_list_gene_anno.at(index_match).cuff_anno_difftrans_pvalue)
				global_list_gene_anno.at(index_match).cuff_anno_difftrans_pvalue = new_cuff_res.p_value;
			global_list_gene_anno.at(index_match).cuff_anno_isoforms.push_back(new_isoform);
		} 
		else
		{
			// no match in annotated list
			getline(infile, info_other);
			outfile_error << info_transID << "\t" << info_geneID << info_other << endl;
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}

// extract chromosome and range from a locus string
int extract_locus(const string locus, string &chr, long &rangeLow, long &rangeHigh)
{
	string info_locus = locus;
	
	size_t found = info_locus.find(':');
	if (found == string::npos)
		return 1;
	chr = info_locus.substr(0, found);
	info_locus = info_locus.substr(found+1);

	found = info_locus.find('-');
	if (found == string::npos)
		return 2;
	string tmpnumber = info_locus.substr(0, found);
	info_locus = info_locus.substr(found+1);
	rangeLow = atol(tmpnumber.c_str());
	
	rangeHigh = atol(info_locus.c_str());

	return 0;
}


int parse_cuffdiff_expr_novel(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open Cuffdiff gene expression result (no annotation)." << endl;
		return 1;
	}

	//ofstream outfile_error( (dir_error+"/cuffdiff_gene_no_annotation_discarded.txt").c_str() );

	string info_geneID, info_testID, info_locus, info_other, titleline, locus_chr;
	long locus_rangeLow, locus_rangeHigh;

	getline(infile, titleline);

	while (infile >> info_testID)
	{
		infile >> info_geneID;

		test_cuff new_cuff_res;

		// find a match in annotated genes
		infile >> info_other; // gene name
		infile >> info_locus; // locus
		infile >> info_other; // sample_1
		infile >> info_other; // sample_2
		infile >> info_other; // status

		infile >> new_cuff_res.fpkm_g1;
		infile >> new_cuff_res.fpkm_g2;
		infile >> new_cuff_res.foldchange;
		if (new_cuff_res.foldchange.compare("inf") == 0 || new_cuff_res.foldchange.compare("-inf") == 0)
			new_cuff_res.foldchange = "";

		infile >> info_other; // test_stat
		infile >> info_other; // p_value before correction

		infile >> new_cuff_res.p_value;
		infile >> info_other; new_cuff_res.significant = str2bool(info_other);
		getline(infile, info_other);
		
		extract_locus(info_locus, locus_chr, locus_rangeLow, locus_rangeHigh); 
		long index_match = locate_gene_by_overlap(global_list_gene_anno, locus_chr, locus_rangeLow, locus_rangeHigh);
		if (index_match >= 0)
		{
			global_list_gene_anno.at(index_match).cuff_no_anno_diffexpr = new_cuff_res.significant;
			global_list_gene_anno.at(index_match).cuff_no_anno_expr_stat = new_cuff_res;
		} 
		else
		{
			// no match in annotated list
			gene new_gene;
			new_gene.geneID = info_geneID;
			new_gene.is_annotated = false;
			new_gene.chrNm = locus_chr;
			new_gene.rangeLow = locus_rangeLow;
			new_gene.rangeHigh = locus_rangeHigh;
			new_gene.cuff_no_anno_diffexpr = new_cuff_res.significant;
			new_gene.cuff_no_anno_expr_stat = new_cuff_res;
			global_list_gene_novel.push_back(new_gene);
		}
	}

	infile.close();
	//outfile_error.close();

	return 0;
}


int parse_cuffdiff_trans_novel(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open Cuffdiff isoform expression result (no annotation)." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/cuffdiff_isoform_no_annotation_discarded.txt").c_str() );

	string info_geneID, info_transID, info_other, titleline, info_locus, locus_chr;
	long locus_rangeLow, locus_rangeHigh;

	getline(infile, titleline);

	while (infile >> info_transID)
	{
		infile >> info_geneID;

		test_cuff new_cuff_res;

		// find a match in annotated genes
		infile >> info_other; // gene name
		infile >> info_locus; // locus
		infile >> info_other; // sample_1
		infile >> info_other; // sample_2
		infile >> info_other; // status

		infile >> new_cuff_res.fpkm_g1;
		infile >> new_cuff_res.fpkm_g2;
		infile >> new_cuff_res.foldchange;
		if (new_cuff_res.foldchange.compare("inf") == 0 || new_cuff_res.foldchange.compare("-inf") == 0)
			new_cuff_res.foldchange = "";

		infile >> info_other; // test_stat
		infile >> info_other; // p_value before correction

		infile >> new_cuff_res.p_value;
		infile >> info_other; new_cuff_res.significant = str2bool(info_other);
		getline(infile, info_other);

		isoform new_isoform;
		new_isoform.id = info_transID;
		new_isoform.cuff = new_cuff_res;

		extract_locus(info_locus, locus_chr, locus_rangeLow, locus_rangeHigh);
		long index_match_anno = locate_gene_by_overlap(global_list_gene_anno, locus_chr, locus_rangeLow, locus_rangeHigh); // here match using the overlap!!
		if (index_match_anno >= 0)
		{
			if (!global_list_gene_anno.at(index_match_anno).cuff_no_anno_difftrans)
				global_list_gene_anno.at(index_match_anno).cuff_no_anno_difftrans = new_cuff_res.significant;
			if (new_cuff_res.p_value < global_list_gene_anno.at(index_match_anno).cuff_no_anno_difftrans_pvalue)
				global_list_gene_anno.at(index_match_anno).cuff_no_anno_difftrans_pvalue = new_cuff_res.p_value;
			global_list_gene_anno.at(index_match_anno).cuff_no_anno_isoforms.push_back(new_isoform);
		} 
		else
		{
			// no match in annotated list
			long index_match_novel = locate_gene_by_ID(global_list_gene_novel, info_geneID); // here match using gene ID because they are all Cufflinks ID
			if (index_match_novel >= 0)
			{
				if (!global_list_gene_novel.at(index_match_novel).cuff_no_anno_difftrans)
					global_list_gene_novel.at(index_match_novel).cuff_no_anno_difftrans = new_cuff_res.significant;
				if (new_cuff_res.p_value < global_list_gene_novel.at(index_match_novel).cuff_no_anno_difftrans_pvalue)
					global_list_gene_novel.at(index_match_novel).cuff_no_anno_difftrans_pvalue = new_cuff_res.p_value;
				global_list_gene_novel.at(index_match_novel).cuff_no_anno_isoforms.push_back(new_isoform); 
			} 
			else
			{
				outfile_error << info_transID << "\t" << info_geneID << endl;
			}			
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}


int parse_cuffdiff_splicing_novel(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open Cuffdiff splicing result (no annotation)." << endl;
		return 1;
	}

	//ofstream outfile_error( (dir_error+"/cuffdiff_splicing_no_annotation_discarded.txt").c_str() );

	string info_geneID, info_testID, info_locus, info_other, titleline, locus_chr;
	long locus_rangeLow, locus_rangeHigh;
	bool is_significant;
	double sqrtJS, p_value;

	getline(infile, titleline);

	while (infile >> info_testID)
	{
		infile >> info_geneID;
		
		// find a match in annotated genes
		infile >> info_other; // gene name
		infile >> info_locus; // locus
		infile >> info_other; // sample_1
		infile >> info_other; // sample_2
		infile >> info_other; // status

		infile >> info_other; // value_1
		infile >> info_other; // value_2
		infile >> sqrtJS;
		//if (new_cuff_res.foldchange.compare("inf") == 0 || new_cuff_res.foldchange.compare("-inf") == 0)
		//	new_cuff_res.foldchange = "";

		infile >> info_other; // test_stat
		infile >> info_other; // p_value before correction

		infile >> p_value;
		infile >> info_other; is_significant = str2bool(info_other);
		getline(infile, info_other);

		extract_locus(info_locus, locus_chr, locus_rangeLow, locus_rangeHigh); 
		long index_match = locate_gene_by_overlap(global_list_gene_anno, locus_chr, locus_rangeLow, locus_rangeHigh);
		if (index_match >= 0)
		{
			global_list_gene_anno.at(index_match).cuff_no_anno_diffsplicing = is_significant;
			global_list_gene_anno.at(index_match).cuff_no_anno_diffsplicing_pvalue = p_value;
			global_list_gene_anno.at(index_match).cuff_no_anno_diffsplicing_sqrtJS = sqrtJS;
		} 
		else
		{
			// no match in annotated list
			long index_match_novel = locate_gene_by_ID(global_list_gene_novel, info_geneID); // here match using gene ID because they are all Cufflinks ID
			if (index_match_novel >= 0)
			{
				global_list_gene_novel.at(index_match_novel).cuff_no_anno_diffsplicing = is_significant;
				global_list_gene_novel.at(index_match_novel).cuff_no_anno_diffsplicing_pvalue = p_value;
				global_list_gene_novel.at(index_match_novel).cuff_no_anno_diffsplicing_sqrtJS = sqrtJS;
			} 	
		}
	}

	infile.close();
	//outfile_error.close();

	return 0;
}


int parse_mts_expr(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open MultiSplice expression result." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/mts_expr_discarded.txt").c_str() );

	string info_geneID, info_other, titleline;

	getline(infile, titleline);

	while (infile >> info_geneID)
	{
		infile >> info_other; // position
		long index_match = locate_gene_by_ID(global_list_gene_anno, info_geneID);

		if (index_match >= 0)
		{
			// find a match in annotated genes
			infile >> global_list_gene_anno.at(index_match).mts_expr_stat; // stat_diff_expr
			infile >> global_list_gene_anno.at(index_match).mts_fold_change;
			infile >> global_list_gene_anno.at(index_match).mts_cov1; // coverage_group1
			infile >> global_list_gene_anno.at(index_match).mts_cov2; // coverage_group2
			infile >> info_other; global_list_gene_anno.at(index_match).mts_diffexpr = str2bool(info_other);
			getline(infile, info_other);
		} 
		else
		{
			// no match in annotated list
			getline(infile, info_other);
			outfile_error << info_geneID << info_other << endl;
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}

int parse_mts_trans(string filename, string dir_error)
{
	ifstream infile(filename.c_str());
	if (!infile.is_open())
	{
		cout << "Error: cannot open MultiSplice transcription result." << endl;
		return 1;
	}

	ofstream outfile_error( (dir_error+"/mts_trans_discarded.txt").c_str() );

	string info_geneID, info_other, titleline;
	
	getline(infile, titleline);

	while (infile >> info_geneID)
	{
		infile >> info_other;
		long index_match = locate_gene_by_ID(global_list_gene_anno, info_geneID);

		if (index_match >= 0)
		{
			// find a match in annotated genes
			infile >> global_list_gene_anno.at(index_match).mts_trans_stat; // stat_diff_expr
			infile >> global_list_gene_anno.at(index_match).mts_sqrtJSD;
			infile >> info_other; // coverage_group1
			infile >> info_other; // coverage_group2
			infile >> info_other; global_list_gene_anno.at(index_match).mts_difftrans = str2bool(info_other);
			getline(infile, info_other);
		} 
		else
		{
			// no match in annotated list
			getline(infile, info_other);
			outfile_error << info_geneID << info_other << endl;
		}
	}

	infile.close();
	outfile_error.close();

	return 0;
}
