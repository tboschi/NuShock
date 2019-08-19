#include "gst.h"

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

	fChain->SetBranchAddress("iev", &iev, &b_iev);
	fChain->SetBranchAddress("neu", &neu, &b_neu);
	fChain->SetBranchAddress("fspl", &fspl, &b_fspl);
	fChain->SetBranchAddress("tgt", &tgt, &b_tgt);
	fChain->SetBranchAddress("Z", &Z, &b_Z);
	fChain->SetBranchAddress("A", &A, &b_A);
	fChain->SetBranchAddress("hitnuc", &hitnuc, &b_hitnuc);
	fChain->SetBranchAddress("hitqrk", &hitqrk, &b_hitqrk);
	fChain->SetBranchAddress("resid", &resid, &b_resid);
	fChain->SetBranchAddress("sea", &sea, &b_sea);
	fChain->SetBranchAddress("qel", &qel, &b_qel);
	fChain->SetBranchAddress("mec", &mec, &b_mec);
	fChain->SetBranchAddress("res", &res, &b_res);
	fChain->SetBranchAddress("dis", &dis, &b_dis);
	fChain->SetBranchAddress("coh", &coh, &b_coh);
	fChain->SetBranchAddress("dfr", &dfr, &b_dfr);
	fChain->SetBranchAddress("imd", &imd, &b_imd);
	fChain->SetBranchAddress("imdanh", &imdanh, &b_imdanh);
	fChain->SetBranchAddress("singlek", &singlek, &b_singlek);
	fChain->SetBranchAddress("nuel", &nuel, &b_nuel);
	fChain->SetBranchAddress("em", &em, &b_em);
	fChain->SetBranchAddress("cc", &cc, &b_cc);
	fChain->SetBranchAddress("nc", &nc, &b_nc);
	fChain->SetBranchAddress("charm", &charm, &b_charm);
	fChain->SetBranchAddress("neut_code", &neut_code, &b_neut_code);
	fChain->SetBranchAddress("nuance_code", &nuance_code, &b_nuance_code);
	fChain->SetBranchAddress("wght", &wght, &b_wght);
	fChain->SetBranchAddress("xs", &xs, &b_xs);
	fChain->SetBranchAddress("ys", &ys, &b_ys);
	fChain->SetBranchAddress("ts", &ts, &b_ts);
	fChain->SetBranchAddress("Q2s", &Q2s, &b_Q2s);
	fChain->SetBranchAddress("Ws", &Ws, &b_Ws);
	fChain->SetBranchAddress("x", &x, &b_x);
	fChain->SetBranchAddress("y", &y, &b_y);
	fChain->SetBranchAddress("t", &t, &b_t);
	fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
	fChain->SetBranchAddress("W", &W, &b_W);
	fChain->SetBranchAddress("EvRF", &EvRF, &b_EvRF);
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
	fChain->SetBranchAddress("pl", &pl, &b_pl);
	fChain->SetBranchAddress("cthl", &cthl, &b_cthl);
	fChain->SetBranchAddress("nfp", &nfp, &b_nfp);
	fChain->SetBranchAddress("nfn", &nfn, &b_nfn);
	fChain->SetBranchAddress("nfpip", &nfpip, &b_nfpip);
	fChain->SetBranchAddress("nfpim", &nfpim, &b_nfpim);
	fChain->SetBranchAddress("nfpi0", &nfpi0, &b_nfpi0);
	fChain->SetBranchAddress("nfkp", &nfkp, &b_nfkp);
	fChain->SetBranchAddress("nfkm", &nfkm, &b_nfkm);
	fChain->SetBranchAddress("nfk0", &nfk0, &b_nfk0);
	fChain->SetBranchAddress("nfem", &nfem, &b_nfem);
	fChain->SetBranchAddress("nfother", &nfother, &b_nfother);
	fChain->SetBranchAddress("nip", &nip, &b_nip);
	fChain->SetBranchAddress("nin", &nin, &b_nin);
	fChain->SetBranchAddress("nipip", &nipip, &b_nipip);
	fChain->SetBranchAddress("nipim", &nipim, &b_nipim);
	fChain->SetBranchAddress("nipi0", &nipi0, &b_nipi0);
	fChain->SetBranchAddress("nikp", &nikp, &b_nikp);
	fChain->SetBranchAddress("nikm", &nikm, &b_nikm);
	fChain->SetBranchAddress("nik0", &nik0, &b_nik0);
	fChain->SetBranchAddress("niem", &niem, &b_niem);
	fChain->SetBranchAddress("niother", &niother, &b_niother);
	fChain->SetBranchAddress("ni", &ni, &b_ni);
	fChain->SetBranchAddress("pdgi", pdgi, &b_pdgi);
	fChain->SetBranchAddress("resc", resc, &b_resc);
	fChain->SetBranchAddress("Ei", Ei, &b_Ei);
	fChain->SetBranchAddress("pxi", pxi, &b_pxi);
	fChain->SetBranchAddress("pyi", pyi, &b_pyi);
	fChain->SetBranchAddress("pzi", pzi, &b_pzi);
	fChain->SetBranchAddress("nf", &nf, &b_nf);
	fChain->SetBranchAddress("pdgf", pdgf, &b_pdgf);
	fChain->SetBranchAddress("Ef", Ef, &b_Ef);
	fChain->SetBranchAddress("pxf", pxf, &b_pxf);
	fChain->SetBranchAddress("pyf", pyf, &b_pyf);
	fChain->SetBranchAddress("pzf", pzf, &b_pzf);
	fChain->SetBranchAddress("pf", pf, &b_pf);
	fChain->SetBranchAddress("cthf", cthf, &b_cthf);
	fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
	fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
	fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
	fChain->SetBranchAddress("vtxt", &vtxt, &b_vtxt);
	fChain->SetBranchAddress("sumKEf", &sumKEf, &b_sumKEf);
	fChain->SetBranchAddress("calresp0", &calresp0, &b_calresp0);
}
