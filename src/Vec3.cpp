/*
	Class to represent a 3-space vector in both cartesian and spherical coordinates. Can perform vector
	addition, subtraction, and dot product. 

	--GWM Dec 2020
*/
#include "Vec3.h"

Vec3::Vec3() {
	m_data[0] = 0.;
	m_data[1] = 0.;
	m_data[2] = 0.;
}

Vec3::Vec3(double x, double y, double z) {
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

Vec3::~Vec3() {}

void Vec3::SetVectorCartesian(double x, double y, double z) {
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

void Vec3::SetVectorSpherical(double r, double theta, double phi) {
	m_data[0] = r*std::cos(phi)*std::sin(theta);
	m_data[1] = r*std::sin(phi)*std::sin(theta);
	m_data[2] = r*std::cos(theta);
}

double Vec3::Dot(const Vec3& rhs) const {
	return GetX()*rhs.GetX() + GetY()*rhs.GetY() + GetZ()*rhs.GetZ();
}

Vec3 Vec3::Cross(const Vec3& rhs) const {
	double x = GetY()*rhs.GetZ() - GetZ()*rhs.GetY();
	double y = GetZ()*rhs.GetX() - GetX()*rhs.GetZ();
	double z = GetX()*rhs.GetY() - GetY()*rhs.GetX();
	return Vec3(x,y,z);
}

inline Vec3 operator*(const int& lhs, const Vec3& rhs) { return Vec3(rhs.GetX()*lhs, rhs.GetY()*lhs, rhs.GetZ()*lhs); };