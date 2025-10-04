#include "Beamspot.h"

Beamspot::Beamspot() = default;

Beamspot::~Beamspot() {
	delete profile;
}

Beamspot::Beamspot(Beamspot&& other) noexcept {
	profile = other.profile;
	xOffset = other.xOffset;
	yOffset = other.yOffset;
	other.profile = nullptr;
}

Beamspot& Beamspot::operator=(Beamspot&& other) noexcept{
	if(this != &other){
		delete profile;
		profile = other.profile;
		xOffset = other.xOffset;
		yOffset = other.yOffset;
		other.profile = nullptr;
	}

	return *this;
}

void Beamspot::SetProfile(BeamProfile* p){
	delete profile;
	profile = p;
}

void Beamspot::SetBeamAxisOffset(double xoff, double yoff){
	xOffset = xoff;
	yOffset = yoff;
}

Vec3 Beamspot::GeneratePoint(double z) const {
	if(!profile){
		return Vec3(xOffset, yOffset, z);//fallback of beam axis if no profile
	}

	auto[dx,dy] = profile->Sample();
	return Vec3(dx+xOffset, dy+yOffset, z);
}