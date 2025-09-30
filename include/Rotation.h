/*
	
	This code taken from "SabreRecon" provided to me by Rachel Shaffer
	This code was originally written by kgh at FSU and rewritten by gwm at FSU
	This code is now being rewritten by zmp at LSU in 5/2025

*/

#ifndef ROTATION_H
#define ROTATION_H

#include "Vec3.h"

class XRotation {
public:
	XRotation();
	XRotation(double ang);
	~XRotation();
	Vec3 Rotate(const Vec3& vector);
	inline void SetAngle(double ang) { m_angle = ang; GenerateMatrix(); };
	inline XRotation GetInverse() { return XRotation(-m_angle); };
	inline Vec3 operator*(const Vec3& vector){
		double x = m_matrix[0][0]*vector[0] + m_matrix[0][1]*vector[1] + m_matrix[0][2]*vector[2];
		double y = m_matrix[1][0]*vector[0] + m_matrix[1][1]*vector[1] + m_matrix[1][2]*vector[2];
		double z = m_matrix[2][0]*vector[0] + m_matrix[2][1]*vector[1] + m_matrix[2][2]*vector[2];
		return Vec3(x, y, z);
	};

private:
	void GenerateMatrix();
	double m_angle;
	double m_matrix[3][3];
};

class YRotation {
public:
	YRotation();
	YRotation(double ang);
	~YRotation();
	Vec3 Rotate(const Vec3& vector);
	inline void SetAngle(double ang) { m_angle = ang; GenerateMatrix(); };
	inline YRotation GetInverse() { return YRotation(-m_angle); };
	inline Vec3 operator*(const Vec3& vector){
		double x = m_matrix[0][0]*vector[0] + m_matrix[0][1]*vector[1] + m_matrix[0][2]*vector[2];
		double y = m_matrix[1][0]*vector[0] + m_matrix[1][1]*vector[1] + m_matrix[1][2]*vector[2];
		double z = m_matrix[2][0]*vector[0] + m_matrix[2][1]*vector[1] + m_matrix[2][2]*vector[2];
		return Vec3(x, y, z);
	}

private:
	void GenerateMatrix();
	double m_angle;
	double m_matrix[3][3];
};

class ZRotation {
public:
	ZRotation();
	ZRotation(double ang);
	~ZRotation();
	Vec3 Rotate(const Vec3& vector);
	inline void SetAngle(double ang) { m_angle = ang; GenerateMatrix();};
	inline ZRotation GetInverse() { return ZRotation(-m_angle); };
	inline Vec3 operator*(const Vec3& vector) {
		double x = m_matrix[0][0]*vector[0] + m_matrix[0][1]*vector[1] + m_matrix[0][2]*vector[2];
		double y = m_matrix[1][0]*vector[0] + m_matrix[1][1]*vector[1] + m_matrix[1][2]*vector[2];
		double z = m_matrix[2][0]*vector[0] + m_matrix[2][1]*vector[1] + m_matrix[2][2]*vector[2];
		return Vec3(x, y, z);
	};

private:
	void GenerateMatrix();
	double m_angle;
	double m_matrix[3][3];
};

#endif