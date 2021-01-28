#ifndef CUTS_H
#define CUTS_H

struct Cut
{
	enum Name {
		E_A,
		P_A,
		T_A,
		TheA,
		PhiA,
		In_A,
		Out_A,
		E_B,
		P_B,
		T_B,
		TheB,
		PhiB,
		In_B,
		Out_B,
		Angle,
		E_0,
		P_0,
		T_0,
		The0,
		Phi0,
		M_0,
		//////////////////
		CosAB,
		aCosAB,
		CircAB,
		atAB0,
		T_AB,
		T_AB0,
		TTA,
		EAB,
		E0Ang
	};

	static Name fromString(const std::string &name) {
		if (name == "E_A")
			return Name::E_A;
		else if (name == "P_A")
			return Name::P_A;
		else if (name == "T_A")
			return Name::T_A;
		else if (name == "TheA")
			return Name::TheA;
		else if (name == "PhiA")
			return Name::PhiA;
		else if (name == "In_A")
			return Name::In_A;
		else if (name == "Out_A")
			return Name::Out_A;
		else if (name == "E_B")
			return Name::E_B;
		else if (name == "P_B")
			return Name::P_B;
		else if (name == "T_B")
			return Name::T_B;
		else if (name == "TheB")
			return Name::TheB;
		else if (name == "PhiB")
			return Name::PhiB;
		else if (name == "In_B")
			return Name::In_B;
		else if (name == "Out_B")
			return Name::Out_B;
		else if (name == "Angle")
			return Name::Angle;
		else if (name == "E_0")
			return Name::E_0;
		else if (name == "P_0")
			return Name::P_0;
		else if (name == "T_0")
			return Name::T_0;
		else if (name == "The0")
			return Name::The0;
		else if (name == "Phi0")
			return Name::Phi0;
		else if (name == "M_0")
			return Name::M_0;
		else if (name == "CosAB")
			return Name::CosAB;
		else if (name == "aCosAB")
			return Name::aCosAB;
		else if (name == "CircAB")
			return Name::CircAB;
		else if (name == "atAB0")
			return Name::atAB0;
		else if (name == "T_AB")
			return Name::T_AB;
		else if (name == "T_AB0")
			return Name::T_AB0;
		else if (name == "TTA")
			return Name::TTA;
		else if (name == "EAB")
			return Name::EAB;
		else if (name == "E0Ang")
			return Name::E0Ang;
		else
			return Name::undefined;
	}

	static std::string toString(Cuts::Name cut) {
		switch (cut) {
			case Name::E_A:
				return "E_A";
			case Name::P_A:
				return "P_A";
			case Name::T_A:
				return "T_A";
			case Name::TheA:
				return "TheA";
			case Name::PhiA:
				return "PhiA";
			case Name::In_A:
				return "In_A";
			case Name::Out_A:
				return "Out_A";
			case Name::E_B:
				return "E_B";
			case Name::P_B:
				return "P_B";
			case Name::T_B:
				return "T_B";
			case Name::TheB:
				return "TheB";
			case Name::PhiB:
				return "PhiB";
			case Name::In_B:
				return "In_B";
			case Name::Out_B:
				return "Out_B";
			case Name::Angle:
				return "Angle";
			case Name::E_0:
				return "E_0";
			case Name::P_0:
				return "P_0";
			case Name::T_0:
				return "T_0";
			case Name::The0:
				return "The0";
			case Name::Phi0:
				return "Phi0";
			case Name::M_0:
				return "M_0";
			case Name::CosAB:
				return "CosAB";
			case Name::aCosAB:
				return "aCosAB";
			case Name::CircAB:
				return "CircAB";
			case Name::atAB0:
				return "atAB0";
			case Name::T_AB:
				return "T_AB";
			case Name::T_AB0:
				return "T_AB0";
			case Name::TTA:
				return "TTA";
			case Name::EAB:
				return "EAB";
			case Name::E0Ang:
				return "E0Ang";
		}
	}

	static constexpr const Name for_nPI0[] = { Name::E_0,
						   Name::T_0,
						   Name::The0,
						   Name::E_A,
						   Name::TheA };

	static constexpr const Name for_nEE[] = { Name::E_A,
						  Name::E_0,
						  Name::M_0,
						  Name::T_A,
						  Name::T_0,
						  Name::TheA,
						  Name::The0,
						  Name::T_AB0,
						  Name::Angle };

	static constexpr const Name for_nMM[] = { Name::E_A,
						  Name::E_0,
						  Name::M_0,
						  Name::T_A,
						  Name::T_0,
						  Name::TheA,
						  Name::The0,
						  Name::T_AB0,
						  Name::Angle };

	static constexpr const Name for_nEM[] = { Name::E_A,
						  Name::E_B,
						  Name::E_0,
						  Name::M_0,
						  Name::T_A,
						  Name::T_B,
						  Name::T_0,
						  Name::TheA,
						  Name::TheB,
						  Name::The0,
						  Name::T_AB0,
						  Name::Angle };

	static constexpr const Name for_EPI[] = { Name::E_A,
						  Name::E_B,
						  Name::E_0,
						  Name::M_0,
						  Name::T_A,
						  Name::T_B,
						  Name::T_0,
						  Name::TheA,
						  Name::TheB,
						  Name::The0,
						  Name::Angle,
						  Name::CosAB,
						  Name::T_AB0 };

	static constexpr const Name for_MPI[] = { Name::E_A,
						  Name::E_B,
						  Name::E_0,
						  Name::M_0,
						  Name::T_A,
						  Name::T_B,
						  Name::T_0,
						  Name::TheA,
						  Name::TheB,
						  Name::The0,
						  Name::Angle,
						  Name::CosAB,
						  Name::T_AB0 };
};

#endif
