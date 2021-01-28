#include "detector/Track.h"

//manual ctor

Track::Track(double x, double y, double z) :
	TVector3(x, y, z),
	kShower(false)
{
}

Track::Track(const TVector3 &v) :
	TVector3(v),
	kShower(false)
{
}

std::ostream & operator<<(std::ostream &os, const Track &t) {
	return os << "<start: (" << t.X() << ", " << t.Y() << ", " << t.Z() << ")"
	          << ", length: " << t.Length()
		  << ", shower: " << std::boolalpha << t.kShower << ">";
}

double Track::Length() const
{
	return std::accumulate(_tracks.begin(), _tracks.end(), 0.,
			[](double l, const std::pair<std::string, double> &tt) {
				return l + tt.second; } );
}

double Track::Length(const std::string &mod) const
{
	if (!_tracks.count(mod))
		return 0.;

	return _tracks.at(mod);
}

double Track::EnergyDeposited() const
{
	return std::accumulate(_calors.begin(), _calors.end(), 0.,
			[](double l, const std::pair<std::string, double> &tt) {
				return l + tt.second; } );
}

double Track::EnergyDeposited(const std::string &mod) const
{
	if (!_calors.count(mod))
		return 0.;

	return _calors.at(mod);
}

double Track::Importance() const
{
	return Length() * EnergyDeposited();
}

double Track::Importance(const std::string &mod) const
{
	return Length(mod) * EnergyDeposited(mod);
}

bool Track::IsShower() const
{
	return kShower;
}

//////// non const functions
//
void Track::SetLength(const std::string &mod, double track)
{
	_tracks[mod] = track;
}

void Track::SetEnergyDeposited(const std::string &mod, double energy)
{
	_calors[mod] = energy;
}

void Track::SetShower(bool shower)
{
	kShower = shower;
}


void Track::SetPosition(double x, double y, double z)
{
	TVector3::SetX(x);
	TVector3::SetY(y);
	TVector3::SetZ(z);
}

void Track::SetPosition(const TVector3 &v)
{
	TVector3::SetX(v.X());
	TVector3::SetY(v.Y());
	TVector3::SetZ(v.Z());
}
