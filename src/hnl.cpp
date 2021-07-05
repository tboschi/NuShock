#include "montecarlo/hnl.h"

hnl::hnl(std::string name) : fChain(NULL), own(false)
{
	New(name);
}

hnl::hnl(TTree *tree) : fChain(NULL), own(false)
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	Init(tree);
}

hnl::~hnl()
{
	if (own)
		delete fChain;
}

size_t hnl::GetEntries()
{
	return fChain->GetEntries();
}

int hnl::GetEntry(size_t entry)
{
	if (fChain)
		return fChain->GetEntry(entry);
	return 0;
}

// create internal TTree of which this class is owner
void hnl::New(std::string name)
{
	// Set branch addresses and branch pointers
	if (!fChain)
		fChain = new TTree(name.c_str(), name.c_str());

	b_True  = fChain->Branch("True", &True, "fTrue/D");
	b_Vert  = fChain->Branch("Vert", Vert, "fVert[3]/D");
	b_W     = fChain->Branch("W", &W, "fTrue/D");
	b_ChID  = fChain->Branch("ChID", &ChID, "fTrue/O");

	b_PdgA  = fChain->Branch("PdgA", &PdgA, "iPdgA/I");
	b_ChA   = fChain->Branch("ChA", &ChA, "iChA/I");
	b_E_A   = fChain->Branch("E_A", &E_A, "fE_A/D");
	b_P_A   = fChain->Branch("P_A", &P_A, "fP_A/D");
	b_T_A   = fChain->Branch("T_A", &T_A, "fT_A/D");
	b_PhiA  = fChain->Branch("PhiA", &PhiA, "fPhiA/D");
	b_TheA  = fChain->Branch("TheA", &TheA, "fTheA/D");
	b_M_A   = fChain->Branch("M_A", &M_A, "fM_A/D");
	b_In_A  = fChain->Branch("In_A", &In_A, "fIn_A/D");
	b_Out_A = fChain->Branch("Out_A", &Out_A, "fOut_A/D");

	b_PdgB  = fChain->Branch("PdgB", &PdgB, "iPdgB/I");
	b_ChB   = fChain->Branch("ChB", &ChB, "iChB/I");
	b_E_B   = fChain->Branch("E_B", &E_B, "fE_B/D");
	b_P_B   = fChain->Branch("P_B", &P_B, "fP_B/D");
	b_T_B   = fChain->Branch("T_B", &T_B, "fT_B/D");
	b_PhiB  = fChain->Branch("PhiB", &PhiB, "fPhiB/D");
	b_TheB  = fChain->Branch("TheB", &TheB, "fTheB/D");
	b_M_B   = fChain->Branch("M_B", &M_B, "fM_B/D");
	b_In_A  = fChain->Branch("In_B", &In_B, "fIn_B/D");
	b_Out_B = fChain->Branch("Out_B", &Out_B, "fOut_B/D");

	b_Angle = fChain->Branch("Angle", &Angle, "fAngle/D");
	b_E_0   = fChain->Branch("E_0", &E_0, "fE_0/D");
	b_P_0   = fChain->Branch("P_0", &P_0, "fP_0/D");
	b_T_0   = fChain->Branch("T_0", &T_0, "fT_0/D");
	b_Phi0  = fChain->Branch("Phi0", &Phi0, "fPhi0/D");
	b_The0  = fChain->Branch("The0", &The0, "fThe0/D");
	b_M_0   = fChain->Branch("M_0", &M_0, "fM_0/D");

	fChain->SetDirectory(0);
	own = true;
}

// bind addresses of external TTree of which this class is not owner
void hnl::Init(TTree *tree)
{
	std::cout << "WARNING: wrapping around TTree " << tree << "; remember to delete it externally!\n";

	// Set branch addresses and branch pointers
	if (!tree)
		return;
	fChain = tree;

	fChain->SetBranchAddress("True", &True, &b_True);
	fChain->SetBranchAddress("Vert", Vert, &b_Vert);
	fChain->SetBranchAddress("W", &W, &b_W);
	fChain->SetBranchAddress("ChID", &ChID, &b_ChID);

	fChain->SetBranchAddress("PdgA", &PdgA, &b_PdgA);
	fChain->SetBranchAddress("ChA", &ChA, &b_ChA);
	fChain->SetBranchAddress("E_A", &E_A, &b_E_A);
	fChain->SetBranchAddress("P_A", &P_A, &b_P_A);
	fChain->SetBranchAddress("T_A", &T_A, &b_T_A);
	fChain->SetBranchAddress("PhiA", &PhiA, &b_PhiA);
	fChain->SetBranchAddress("TheA", &TheA, &b_TheA);
	fChain->SetBranchAddress("M_A", &M_A, &b_M_A);
	fChain->SetBranchAddress("In_A", &In_A, &b_In_A);
	fChain->SetBranchAddress("Out_A", &Out_A, &b_Out_A);

	fChain->SetBranchAddress("PdgB", &PdgB, &b_PdgB);
	fChain->SetBranchAddress("ChB", &ChB, &b_ChB);
	fChain->SetBranchAddress("E_B", &E_B, &b_E_B);
	fChain->SetBranchAddress("P_B", &P_B, &b_P_B);
	fChain->SetBranchAddress("T_B", &T_B, &b_T_B);
	fChain->SetBranchAddress("PhiB", &PhiB, &b_PhiB);
	fChain->SetBranchAddress("TheB", &TheB, &b_TheB);
	fChain->SetBranchAddress("M_B", &M_B, &b_M_B);
	fChain->SetBranchAddress("In_B", &In_B, &b_In_B);
	fChain->SetBranchAddress("Out_B", &Out_B, &b_Out_B);

	fChain->SetBranchAddress("Angle", &Angle, &b_Angle);
	fChain->SetBranchAddress("E_0", &E_0, &b_E_0);
	fChain->SetBranchAddress("P_0", &P_0, &b_P_0);
	fChain->SetBranchAddress("T_0", &T_0, &b_T_0);
	fChain->SetBranchAddress("Phi0", &Phi0, &b_Phi0);
	fChain->SetBranchAddress("The0", &The0, &b_The0);
	fChain->SetBranchAddress("M_0", &M_0, &b_M_0);

	own = false;
}

//Load tree
void hnl::Fill(double weight, bool chid,
	       const Tracker::Event &e0, const Tracker::Event &e1, const Tracker::Event &e2)
{
	// even rate -> montecarlo weight
	W = weight;
	ChID = chid;

	True = e0.first.E();
	Vert[0]  = e0.second[0];
	Vert[1]  = e0.second[1];
	Vert[2]  = e0.second[2];

	PdgA  = e1.first.Pdg();
	ChA   = e1.first.Q();
	E_A   = e1.first.E();
	P_A   = e1.first.P();
	T_A   = e1.first.Pt();
	TheA  = e1.first.Theta();
	PhiA  = e1.first.Phi();
	M_A   = e1.first.M();
	Out_A = e1.second.Length("out");
	In_A  = e1.second.Length() - Out_A;

	PdgB  = e2.first.Pdg();
	ChB   = e2.first.Q();
	E_B   = e2.first.E();
	P_B   = e2.first.P();
	T_B   = e2.first.Pt();
	TheB  = e2.first.Theta();
	PhiB  = e2.first.Phi();
	M_B   = e2.first.M();
	Out_B = e2.second.Length("out");
	In_B  = e2.second.Length() - Out_B;

	Angle = e1.first.Vect().Angle(e2.first.Vect());

	TLorentzVector reco = e1.first + e2.first;
	E_0  = reco.E();
	P_0  = reco.P();
	T_0  = reco.Pt();
	The0 = reco.Theta();
	Phi0 = reco.Phi();
	M_0  = reco.M();

	Fill();
}

void hnl::Fill()
{
	fChain->Fill();
}

void hnl::Write()
{
	fChain->Write();
}

TTree * hnl::Chain()
{
	return fChain;
}

double hnl::Events() {
	double evts = 0;
	for (size_t i = 0; i < GetEntries(); ++i)
		evts += W;

	return evts;
}
