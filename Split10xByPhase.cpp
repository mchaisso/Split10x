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

#include "BarcodeTools.h"

using namespace std;
int RemoveInactiveBarcodes(map<string,int> & bc, int curDNAPos, int dist) {
	map<string, int>::iterator bcIt;
	bcIt = bc.begin();
	int before = bc.size();
	int prev = -1;
	int minPos = -1;
	while (bcIt != bc.end()) {
		if (curDNAPos - bcIt->second > dist) {
			map<string, int>::iterator next = bcIt;
			next++;
			bc.erase(bcIt);
			bcIt = next;
		}
		else {
			if (minPos == -1 or minPos > bcIt->second) {
				minPos = bcIt->second;
			}
			++bcIt;
		}
	}
	return minPos;
}

class PhaseQuery {
public:
	vector< pair<int, int> > intervals;
	vector<int> phase;
	int FindPhase(int pos) {
		int i;
		while (i < intervals.size()) {
			if (intervals[i].first <= pos and pos < intervals[i].second) {
				return phase[i];
			}
			else if (pos < intervals.size() - 1 and intervals[i].second <= pos and pos < intervals[i+1].first) {
				// Return middle
				int middle = (intervals[i].second + intervals[i+1].first)/2;
				if (pos < middle) {
					return phase[i];
				}
				else {
					return phase[i+1];
				}
			}
			i+=1;
		}
		return 3;
	}

};
	
void MakePhaseQuery(vector<int> &phase, vector<string> &chrom, vector<int> &pos, map<string, PhaseQuery> &chromQueries, int minBlockSize=20000) {
	int i =0;
	while (i < phase.size()) {
		if (i < phase.size() - 1 and phase[i] != phase[i+1]) {
			i+=1;
			continue;
		}
		if (chromQueries.find(chrom[i]) == chromQueries.end()) {
			chromQueries[chrom[i]] = PhaseQuery();
		}
		if ((pos[i+1] - pos[i]) > minBlockSize ) {
			chromQueries[chrom[i]].intervals.push_back(pair<int,int>(pos[i], pos[i+1]));
			chromQueries[chrom[i]].phase.push_back(phase[i]);
			cerr << "switching " << chrom[i] << " " << pos[i] << " " << pos[i+1] << " " << phase[i] << endl;
		}
		i+=2;
	}
	chromQueries["chrX"] = PhaseQuery();
	chromQueries["chrX"].intervals.push_back(pair<int,int>(0,100000000));
	chromQueries["chrX"].phase.push_back(2);
}

int main(int argc, char* argv[]) {

	if (argc < 5) {
		cout << "Usage: split10x reads.bam phased-vcf base dist [-v] [ --switch switchfile ]" << endl;
		cout << "switchfile is a file of phase\tchrom\tpos, one line per endpoint, " << endl
				 << "of all endpoints between switch errors." << endl;
	  exit(1);
	}
			
	string bamFileName = argv[1];
	string vcfFileName = argv[2];
	string hetBase     = argv[3];
	int    dist        = atoi(argv[4]);
	bool verbose = false;
	int argi = 5;
	string switchFileName = "";
	while (argi < argc) {
		if (string(argv[argi]) == "-v") {
			verbose = true;
		}
		else if (string(argv[argi]) == "--switch") {
			switchFileName = argv[++argi];
		}
		else {
			cout << "Invalid option " << argv[argi];
			exit(1);
		}
		argi++;
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
	vector<int> switchHap;
	vector<string> switchChrom;
	vector<int> switchPos;
	map<string, PhaseQuery> phaseQuery;
	bool checkPhase = false;
	if (switchFileName != "") {
		ifstream switchFile(switchFileName.c_str());
		string line;
		while (switchFile) {
			getline(switchFile, line);
			stringstream strm(line);
			int hap, pos;
			string chrom;
			if (!(strm >> hap >> chrom >> pos)) {
				break;
			}
			else {
				switchHap.push_back(hap);
				switchChrom.push_back(chrom);
				switchPos.push_back(pos);
			}
		}
		MakePhaseQuery(switchHap, switchChrom, switchPos, phaseQuery);
		checkPhase = true;
	}
		
	
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

	int res;
	res=sam_hdr_write(hap0, header);
	res=sam_hdr_write(hap1, header);
	
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
		int bamRetval;
		int h0BCPos = 0;
		int h1BCPos = 0;
		int nUnassigned = 0;
		int nAssigned = 0;
		int nAuto = 0;
		while ((bamRetval = bam_itr_next(readsFile, bamIter, read)) > 0) {
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
			if (curDNAPos - h0BCPos > dist) {
				h0BCPos = RemoveInactiveBarcodes(activeH0Barcodes, curDNAPos, dist);
			}
			if (curDNAPos - h1BCPos > dist) {
				h1BCPos = RemoveInactiveBarcodes(activeH1Barcodes, curDNAPos, dist);			
			}
			//
			// Advance the list of active elements
			//
			int phase = 0;
			if (checkPhase) {

				string chrom = string(header->target_name[read->core.tid]);
				if (phaseQuery.find(chrom) != phaseQuery.end()) {
					phase = phaseQuery[chrom].FindPhase(curDNAPos);
				}
			}
			
			while (curVarPos - curDNAPos < dist and
						 tbx_itr_next(vcfFile, tabix, vcfIter, &str) >=0 ) {
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
							int a0 = bcf_gt_allele(gt_arr[0]);
							int a1 = bcf_gt_allele(gt_arr[1]);
							//							cerr << p0 << " " << p1 << " " << a0 << " " << a1 << endl;
							if (p0 or p1) {
								string gt;
							
								if (a0 != a1) {
									if (a0 == 0 and a1 == 1) {
										gt= "0|1";
									}
									else if (a0== 1 and a1== 0) {
										gt = "1|0";
									}
									else {
										gt = "UNK";
									}
									if (a0 == 0 ) {
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
									if (verbose) {
										int u1, u2;
										CountUnique(activeH0Barcodes, activeH1Barcodes, u1,u2);
										cerr << "# Unique " << u1 << "/" << activeH0Barcodes.size() << "\t" << u2 << "/" << activeH1Barcodes.size() << endl;
									}
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
			if (activeH0Barcodes.find(readBarcode) == activeH0Barcodes.end() or
					activeH1Barcodes.find(readBarcode) == activeH1Barcodes.end()) {
				if (activeH0Barcodes.find(readBarcode) != activeH0Barcodes.end()) {
					if (checkPhase) {
						if (phase == 0) {
							nAssigned++;
							res=bam_write1(hap0->fp.bgzf, read);
							nHap0++;
						}
						else if (phase == 1) {
							nAssigned++;
							res=bam_write1(hap1->fp.bgzf, read);
							nHap1++;
						}
						else if (phase == 2) {
							nAssigned++;
							nAuto++;
							res=bam_write1(hap0->fp.bgzf, read);
							nHap0++;
							res=bam_write1(hap1->fp.bgzf, read);
							nHap1++;
						}							
						else {
							nUnassigned++;
						}
					}
					else {
						res=bam_write1(hap0->fp.bgzf, read);
						nAssigned++;
						nHap0++;						
					}
				}
				else if (activeH1Barcodes.find(readBarcode) != activeH1Barcodes.end()) {
					if (checkPhase) {
						if (phase == 0) {
							nAssigned++;
							res=bam_write1(hap1->fp.bgzf, read);
							nHap1++;
						}
						else if (phase == 1) {
							nAssigned++;
							res=bam_write1(hap0->fp.bgzf, read);
							nHap0++;
						}
						else if (phase == 2) {
							nAssigned++;
							nAuto++;
							res=bam_write1(hap0->fp.bgzf, read);
							nHap0++;
							res=bam_write1(hap1->fp.bgzf, read);
							nHap1++;
						}							
						
						else {
							nUnassigned++;
						}
					}
					else {
						res=bam_write1(hap1->fp.bgzf, read);
						nAssigned++;
						nHap1++;						
					}

				}
				else {
					++nUnphased;
				}
			}
			if (nReads % 1000000 == 0) {
				cerr << header->target_name[targetId] << "\t" << curDNAPos << "\t" << nReads << "\t"
						 << nAssigned << "\t" << nHap0 << "\t" << nHap1 << "\t" << nAuto << "\t"
						 << activeH0Barcodes.size() << "\t" << activeH1Barcodes.size() << endl;
			}
		}
		cerr << "Processed " << header->target_name[targetId] << " " << nReads << "\t" << nUnphased << "\t" << nAssigned << "\t" << nUnassigned << endl;
	}

	bam_destroy1(read);
	bam_hdr_destroy(header);
	sam_close(readsFile);
	sam_close(hap0);
	sam_close(hap1);
		

	
}
