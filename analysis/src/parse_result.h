
#ifndef PARSE_RESULT
#define PARSE_RESULT

#include "analysis.h"
#include "gene_entry.h"

int parse_annotationGTF(string filename);
int parse_dfs_expr(string filename, string dir_error);
int parse_dfs_trans(string filename, string dir_error);
int parse_cuffdiff_expr(string filename, string dir_error);
int parse_cuffdiff_trans(string filename, string dir_error);
int parse_cuffdiff_expr_novel(string filename, string dir_error);
int parse_cuffdiff_trans_novel(string filename, string dir_error);
int parse_cuffdiff_splicing_novel(string filename, string dir_error);
int parse_mts_expr(string filename, string dir_error);
int parse_mts_trans(string filename, string dir_error);

#endif

