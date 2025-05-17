/*
	Class to represent a 3-space vector in both cartesian and spherical coordinates. Can perform vector
	addition, subtraction, and dot product. 

	--GWM Dec 2020
*/
#ifndef VEC3_H
#define VEC3_H

#include <cmath>

class Vec3 {
public:
	Vec3();
	Vec3(double x, double y, double z);
	~Vec3();

	void SetVectorCartesian(double x, double y, double z);
	void SetVectorSpherical(double r, double theta, double phi);
	inline double GetX() const { return m_data[0]; };
	inline double GetY() const { return m_data[1]; };
	inline double GetZ() const { return m_data[2]; };
	inline double GetRho() const { return std::sqrt(std::pow(m_data[0], 2.0) + std::pow(m_data[1], 2.0)); };
	inline double GetR() const { return std::sqrt(std::pow(m_data[0], 2.0) + std::pow(m_data[1], 2.0) + std::pow(m_data[2], 2.0)); }
	inline double GetTheta() const { return Atan2(GetRho(), GetZ()); };
	inline double GetPhi() const {
		double phi = Atan2(GetY(), GetX());
		if(phi < 0) phi += M_PI*2.0;
		return phi;
	};

	inline const double operator[](int index) const { return index>2 || index<0 ? 0.0 : m_data[index]; };
	inline Vec3& operator=(const Vec3& rhs) { SetVectorCartesian(rhs.GetX(), rhs.GetY(), rhs.GetZ()); return *this; };
	inline Vec3 operator+(const Vec3& rhs) const { return Vec3(this->GetX()+rhs.GetX(), this->GetY()+rhs.GetY(), this->GetZ()+rhs.GetZ()); };
	inline Vec3 operator-(const Vec3& rhs) const { return Vec3(this->GetX()-rhs.GetX(), this->GetY()-rhs.GetY(), this->GetZ()-rhs.GetZ()); };
	inline Vec3 operator*(const int& rhs) const { return Vec3(GetX()*rhs,GetY()*rhs,GetZ()*rhs);};


	double Dot(const Vec3& rhs) const;
	Vec3 Cross(const Vec3& rhs) const;

	double Mag2() const { return GetX()*GetX() + GetY()*GetY() + GetZ()*GetZ(); };
	double Mag() const { return sqrt(Mag2()); };



private:

	//Use instead of std::atan2. Better control over values close to x=0
	inline double Atan2(double y, double x) const {
		if(x != 0.0) return std::atan2(y, x);
		else if(y > 0.0) return M_PI/2.0;
		else if(y < 0.0) return -M_PI/2.0;
		else return 0.0;
	}

	double m_data[3];

};

inline Vec3 operator*(const int& lhs, const Vec3& rhs);

#endif