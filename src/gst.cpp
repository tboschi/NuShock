#include "montecarlo/gst.h"

gst::gst(TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	Init(tree);
}

gst::~gst()
{
	delete fChain->GetCurrentFile();
}

int gst::GetEntries()
{
	return fChain->GetEntries();
}

int gst::GetEntry(long entry)
{
	// Read contents of entry.
	if (!fChain)
		return 0;
	else 
		return fChain->GetEntry(entry);
}

void gst::Init(TTree *tree)
{

	// Set branch addresses and branch pointers
	if (!tree)
		return;
	fChain = tree;

	fChain->SetBranchAddress("neu", &neu, &b_neu);
	fChain->SetBranchAddress("qel", &qel, &b_qel);
	fChain->SetBranchAddress("mec", &mec, &b_mec);
	fChain->SetBranchAddress("res", &res, &b_res);
	fChain->SetBranchAddress("dis", &dis, &b_dis);
	fChain->SetBranchAddress("coh", &coh, &b_coh);
	fChain->SetBranchAddress("cc", &cc, &b_cc);
	fChain->SetBranchAddress("nc", &nc, &b_nc);
	fChain->SetBranchAddress("wght", &wght, &b_wght);
	fChain->SetBranchAddress("Ev", &Ev, &b_Ev);
	fChain->SetBranchAddress("pxv", &pxv, &b_pxv);
	fChain->SetBranchAddress("pyv", &pyv, &b_pyv);
	fChain->SetBranchAddress("pzv", &pzv, &b_pzv);
	fChain->SetBranchAddress("En", &En, &b_En);
	fChain->SetBranchAddress("pxn", &pxn, &b_pxn);
	fChain->SetBranchAddress("pyn", &pyn, &b_pyn);
	fChain->SetBranchAddress("pzn", &pzn, &b_pzn);
	fChain->SetBranchAddress("El", &El, &b_El);
	fChain->SetBranchAddress("pxl", &pxl, &b_pxl);
	fChain->SetBranchAddress("pyl", &pyl, &b_pyl);
	fChain->SetBranchAddress("pzl", &pzl, &b_pzl);
	fChain->SetBranchAddress("nf", &nf, &b_nf);
	fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
	fChain->SetBranchAddress("Ef", Ef, &b_Ef);
	fChain->SetBranchAddress("pxf", pxf, &b_pxf);
	fChain->SetBranchAddress("pyf", pyf, &b_pyf);
	fChain->SetBranchAddress("pzf", pzf, &b_pzf);
	fChain->SetBranchAddress("pf", pf, &b_pf);
}
