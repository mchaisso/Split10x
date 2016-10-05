#include <stdint.h>
#include "vcf.h"
#include "sam.h"
#include "tbx.h"

#include <stdlib.h>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <map>

#include <sstream>
#include <numeric>
#include <iomanip>

using namespace std;
string ExtractBarcode(char *str) {
	int len = strlen(str);
	int i;
	for (i = 0; i < len; i++) {
		if (str[i] == '-' or str[i] == '_') {
			return string(str,i);
		}
	}
	return "";
}

int CountUnique(map<string,int> &l1, map<string,int> &l2, int &u1, int &u2) {
	map<string,int>::iterator lit;
	u1 = 0;
	u2 = 0;
	for (lit = l1.begin(); lit != l1.end(); ++lit) {
		if (l2.find(lit->first) == l2.end()) {
			++u1;
		}
	}
	for (lit = l2.begin(); lit != l2.end(); ++lit) {
		if (l1.find(lit->first) == l1.end()) {
			++u2;
		}
	}
}
		
	
void ParseBarcodeList(char* str, vector<string> &refBCList, vector<string> &altBCList) {
	int s, e;
	bool parseRef = true;
	for (s = 0; str[s] != '\0'; s++) {
		for (e = s; str[e] != '\0'; e++) {
			if (str[e] == '-') {
				if (e > s) {
					if (parseRef == true) {
						refBCList.push_back(string(&str[s], e-s));
					}
					else {
						altBCList.push_back(string(&str[s],e-s));
					}
				}
				break;
			}
		}
		s = e;
		for (; str[s] != '\0'; s++) {
			if (str[s] == ',') {
				parseRef = false;
				break;
			}
			if (str[s] == ';') {
 				break;
			}
		}
		if (str[s] == '\0') {
			return;
		}
	}
}
void RemoveInactiveBarcodes(map<string,int> & bc, int curDNAPos, int dist) {
	map<string, int>::iterator bcIt;
	bcIt = bc.begin();
	int before = bc.size();
	while (bcIt != bc.end()) {
		if (curDNAPos - bcIt->second > dist) {
			map<string, int>::iterator next = bcIt;
			next++;
			bc.erase(bcIt);
			bcIt = next;
		}
		else {
			++bcIt;
		}
	}
}

int main(int argc, char* argv[]) {

	if (argc < 5) {
		cout << "Usage: split10x reads.bam phased-vcf base dist" << endl;
		exit(1);
	}
			
	string bamFileName = argv[1];
	string vcfFileName = argv[2];
	string hetBase     = argv[3];
	int    dist        = atoi(argv[4]);
	bool verbose = false;
	if (argc > 5 and string(argv[5]) == "-v") {
		verbose = true;
	}
	hts_itr_t *iter=NULL;
 	hts_idx_t *bamIndex=NULL;
	samFile *readsFile = NULL;
	bam1_t *read= NULL;
	bam_hdr_t *header = NULL;


	htsFile *vcfFile = hts_open(vcfFileName.c_str(), "r");
	bcf_hdr_t * vcfHeader = bcf_hdr_read(vcfFile);
	
	bcf1_t *var = bcf_init();
	hts_idx_t *vcfIdx = bcf_index_load(vcfFileName.c_str());


	tbx_t *tabix = tbx_index_load(vcfFileName.c_str());


	string hap0Name = hetBase + ".0.bam";
	string hap1Name = hetBase + ".1.bam";

		
	
	readsFile = hts_open(bamFileName.c_str(), "r");
	if(readsFile==NULL)
		return -1;
	if ((header = sam_hdr_read(readsFile)) == 0)
		return -1;
	bamIndex = bam_index_load(bamFileName.c_str());
	if (bamIndex==NULL)
		return -1;

	htsFile *hap0 = hts_open(hap0Name.c_str(), "wb");
	htsFile *hap1 = hts_open(hap1Name.c_str(), "wb");

	sam_hdr_write(hap0, header);
	sam_hdr_write(hap1, header);
	
	read = bam_init1();	
	int targetId;
	for (targetId = 0; targetId < header->n_targets; targetId++) {
		map<string, int> activeH0Barcodes, activeH1Barcodes;
		int curVarPos = 0;
		int curDNAPos = 0;
		
		hts_itr_t *vcfIter = tbx_itr_querys(tabix, header->target_name[targetId]);


		kstring_t str={0,0,0};

		hts_itr_t* bamIter = bam_itr_querys(bamIndex, header, header->target_name[targetId]);
		int nReads =0 ;
		int nHap0=0, nHap1=0;
		int prevVCFPos = 0, prevDNAPos=0;
		int prevHap0 = 0, prevHap1 = 0;
		int nUnphased = 0;
		while (bam_itr_next(readsFile, bamIter, read) >= 0) {
			uint8_t *aux = bam_aux_get(read, "BX");
			++nReads;
			if (aux == 0 ) {
				continue;
			}
			

			char *bcStr = bam_aux2Z(aux);
			string readBarcode = ExtractBarcode(bcStr);
			curDNAPos = read->core.pos;
			//
			// First purge inactive elements.
			//
			RemoveInactiveBarcodes(activeH0Barcodes, curDNAPos, dist);
			RemoveInactiveBarcodes(activeH1Barcodes, curDNAPos, dist);			
			//
			// Advance the list of active elements
			//
			while (curVarPos - curDNAPos < dist and
						 tbx_itr_next(vcfFile, tabix, vcfIter, &str)) {
				int ndst = 0; char **dst = NULL;
				int i;
				if (str.l > 0) {
					vcf_parse(&str, vcfHeader, var);
				
					if ( bcf_get_format_string(vcfHeader, var, "BX", &dst, &ndst) > 0 ) {
						for (i=0; i<bcf_hdr_nsamples(vcfHeader); i++) {
							vector<string> refBCList, altBCList;

							ParseBarcodeList(dst[i], refBCList, altBCList);
							free(dst[i]);
							//
							// Add/update these barcodes
							//
							int b;
							int ngt, *gt_arr = NULL, ngt_arr = 0;
							int *gt=NULL;;

							ngt = bcf_get_genotypes(vcfHeader, var, &gt_arr, &ngt_arr);
							int refAllele = 0;
							int p0 = bcf_gt_is_phased(gt_arr[0]);
							int p1 = bcf_gt_is_phased(gt_arr[1]);
							if (p0 or p1) {
								int p0 = bcf_gt_allele(gt_arr[0]);
								int p1 = bcf_gt_allele(gt_arr[1]);
								string gt;
								if (p0 == 0 and p1 == 1) {
									gt= "0|1";
								}
								else if (p0== 1 and p1== 0) {
									gt = "1|0";
								}
								else {
									gt = "UNK";
								}

								if (p0 == 0) {
									for (b = 0; b < refBCList.size(); b++) {
										activeH0Barcodes[refBCList[b]] = var->pos;
									}
									for (b = 0; b < altBCList.size(); b++) {
										activeH1Barcodes[altBCList[b]] = var->pos;
									}
								}
								else {
									for (b = 0; b < altBCList.size(); b++) {
										activeH0Barcodes[altBCList[b]] = var->pos;
									}
									for (b = 0; b < refBCList.size(); b++) {
										activeH1Barcodes[refBCList[b]] = var->pos;
									}
								}
								int u1, u2;
								CountUnique(activeH0Barcodes, activeH1Barcodes, u1,u2);
								if (verbose) {
									cerr << "# Unique " << u1 << "/" << activeH0Barcodes.size() << "\t" << u2 << "/" << activeH1Barcodes.size() << endl;
								}
							}
						}
					}
					if (verbose) {
						cerr << "advanced to var " << var->pos << " dna: " << curDNAPos << " reads: " << nReads
								 << " unp,h0,h1: " << nUnphased << "," << nHap0 << "," << nHap1 << "\tdeltas:\t"
								 << var->pos - prevVCFPos << " " << curDNAPos - prevDNAPos
								 << "\t" << nHap0 - prevHap0 << "\t" << nHap1 - prevHap1 << endl;
					}
					prevVCFPos = var->pos;
					prevDNAPos = curDNAPos;
					prevHap0 = nHap0;
					prevHap1 = nHap1;
					free(str.s);
					curVarPos = var->pos;
				}
				str = {0,0,0};

			}
			if (activeH0Barcodes.find(readBarcode) != activeH0Barcodes.end() and
					activeH1Barcodes.find(readBarcode) != activeH1Barcodes.end()) {

			}
			else if (activeH0Barcodes.find(readBarcode) != activeH0Barcodes.end()) {
				bam_write1(hap0->fp.bgzf, read);
				nHap0++;
			}
			else if (activeH1Barcodes.find(readBarcode) != activeH1Barcodes.end()) {
				bam_write1(hap1->fp.bgzf, read);
				nHap1++;
			}
			else {
				++nUnphased;
			}
			if (nReads % 1000000 == 0) {
				cerr << "Processed " << nReads << endl;
			}
		}
		cerr << "Processed " << header->target_name[targetId] << " " << nReads << "\t" << nUnphased << endl;
	}



	
	bam_destroy1(read);
	bam_hdr_destroy(header);
	sam_close(readsFile);

	
	

}
