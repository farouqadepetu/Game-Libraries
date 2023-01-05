#include "FADirect3DMath.h"

#if defined(_DEBUG)
#include <iostream>
#endif


#define EPSILON 1e-6
#define PI 3.14159265

namespace FAD3DMath
{
	bool compareFloats(float x, float y, float epsilon)
	{
		float diff = fabs(x - y);
		//exact epsilon
		if (diff < epsilon)
		{
			return true;
		}

		//adapative epsilon
		return diff <= epsilon * std::max(fabs(x), fabs(y));
	}

	bool compareDoubles(double x, double y, double epsilon)
	{
		double diff = fabs(x - y);
		//exact epsilon
		if (diff < epsilon)
		{
			return true;
		}

		//adapative epsilon
		return  diff <= epsilon * std::max(fabs(x), fabs(y));
	}

	//--------------------------------------------------------------------------------------
	//Constructors

	Vector2D::Vector2D() : m_x{ 0.0f }, m_y{ 0.0f }
	{}

	Vector2D::Vector2D(float x, float y) : m_x{ x }, m_y{ y }
	{}

	//--------------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------------
	//Getters and Setters

	float Vector2D::x() const
	{
		return m_x;
	}

	float Vector2D::y() const
	{
		return m_y;
	}

	void Vector2D::setX(float x)
	{
		m_x = x;
	}

	void Vector2D::setY(float y)
	{
		m_y = y;
	}

	//--------------------------------------------------------------------------------------


	//--------------------------------------------------------------------------------------
	//Memeber functions

	Vector2D& Vector2D::operator+=(const Vector2D& b)
	{
		this->m_x += (double)b.m_x;
		this->m_y += (double)b.m_y;

		return *this;
	}

	Vector2D& Vector2D::operator-=(const Vector2D& b)
	{
		this->m_x -= (double)b.m_x;
		this->m_y -= (double)b.m_y;

		return *this;
	}

	Vector2D& Vector2D::operator*=(const float& k)
	{
		this->m_x *= (double)k;
		this->m_y *= (double)k;

		return *this;
	}

	Vector2D& Vector2D::operator/=(const float& k)
	{
		if (compareFloats(k, 0.0f, EPSILON))
		{
			return *this;
		}

		this->m_x /= (double)k;
		this->m_y /= (double)k;

		return *this;
	}

	//-------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------
	//Non-member functions

	bool zeroVector(const Vector2D& a)
	{
		if (compareFloats(a.x(), 0.0f, EPSILON) && compareFloats(a.y(), 0.0f, EPSILON))
		{
			return true;
		}

		return false;
	}

	Vector2D operator+(const Vector2D& a, const Vector2D& b)
	{
		return Vector2D((double)a.x() + b.x(), (double)a.y() + b.y());
	}

	Vector2D operator-(const Vector2D& v)
	{
		return Vector2D(-v.x(), -v.y());
	}

	Vector2D operator-(const Vector2D& a, const Vector2D& b)
	{
		return Vector2D((double)a.x() - b.x(), (double)a.y() - b.y());
	}

	Vector2D operator*(const Vector2D& a, const float& k)
	{
		return Vector2D((double)a.x() * k, (double)a.y() * k);
	}

	Vector2D operator*(const float& k, const Vector2D& a)
	{
		return Vector2D((double)k * a.x(), (double)k * a.y());
	}

	Vector2D operator/(const Vector2D& a, const float& k)
	{
		if (compareFloats(k, 0.0f, EPSILON))
		{
			return Vector2D();
		}

		return Vector2D(a.x() / (double)k, a.y() / (double)k);
	}

	//a dot b = axbx + ayby
	float dotProduct(const Vector2D& a, const Vector2D& b)
	{
		return (double)a.x() * b.x() + (double)a.y() * b.y();
	}

	//length(v) = sqrt(vx^2 + vy^2)
	float length(const Vector2D& v)
	{
		return sqrt((double)v.x() * v.x() + (double)v.y() * v.y());
	}

	//norm(v) = v / length(v) == (vx / length(v), vy / length(v))
	Vector2D norm(const Vector2D& v)
	{
		//v is the zero vector
		if (zeroVector(v))
		{
			return v;
		}

		double mag = length(v);

		return Vector2D(v.x() / mag, v.y() / mag);
	}

	//v = (r, theta)
	//x = rcos((theta)
	//y = rsin(theta)
	Vector2D PolarToCartesian(const Vector2D& v)
	{
		float angle = v.y() * PI / 180.0f;

		return Vector2D(v.x() * cos(angle), v.x() * sin(angle));
	}

	//v = (x, y)
	//r = sqrt(vx^2 + vy^2)
	//theta = arctan(y / x)
	Vector2D CartesianToPolar(const Vector2D& v)
	{
		if (compareFloats(v.x(), 0.0f, EPSILON))
		{
			return v;
		}

		double theta{ atan2(v.y(), v.x()) * 180.0 / PI };
		return Vector2D(length(v), theta);
	}

	//Projb(a) = (a dot b)b
	//normalize b before projecting
	Vector2D Projection(const Vector2D& a, const Vector2D& b)
	{
		Vector2D normB(norm(b));
		return Vector2D(dotProduct(a, normB) * normB);
	}


#if defined(_DEBUG)
	void print(const Vector2D& v)
	{
		std::cout << "(" << v.x() << ", " << v.y() << ")";
	}
#endif


	//--------------------------------------------------------------------------------------
	//Constructors

	Vector3D::Vector3D() : m_x{ 0.0f }, m_y{ 0.0f }, m_z{ 0.0f }
	{}

	Vector3D::Vector3D(float x, float y, float z) : m_x{ x }, m_y{ y }, m_z{ z }
	{}

	//--------------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------------
	//Getters and Setters

	float Vector3D::x() const
	{
		return m_x;
	}

	float Vector3D::y() const
	{
		return m_y;
	}

	float Vector3D::z() const
	{
		return m_z;
	}

	void Vector3D::setX(float x)
	{
		m_x = x;
	}

	void Vector3D::setY(float y)
	{
		m_y = y;
	}

	void Vector3D::setZ(float z)
	{
		m_z = z;
	}

	//--------------------------------------------------------------------------------------


	//--------------------------------------------------------------------------------------
	//Memeber functions

	Vector3D& Vector3D::operator+=(const Vector3D& b)
	{
		this->m_x += (double)b.m_x;
		this->m_y += (double)b.m_y;
		this->m_z += (double)b.m_z;

		return *this;
	}

	Vector3D& Vector3D::operator-=(const Vector3D& b)
	{
		this->m_x -= (double)b.m_x;
		this->m_y -= (double)b.m_y;
		this->m_z -= (double)b.m_z;

		return *this;
	}

	Vector3D& Vector3D::operator*=(const float& k)
	{
		this->m_x *= (double)k;
		this->m_y *= (double)k;
		this->m_z *= (double)k;

		return *this;
	}

	Vector3D& Vector3D::operator/=(const float& k)
	{
		if (compareFloats(k, 0.0f, EPSILON))
		{
			return *this;
		}

		this->m_x /= (double)k;
		this->m_y /= (double)k;
		this->m_z /= (double)k;

		return *this;
	}

	//-------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------
	//Non-member functions

	bool zeroVector(const Vector3D& a)
	{
		if (compareFloats(a.x(), 0.0f, EPSILON) && compareFloats(a.y(), 0.0f, EPSILON) &&
			compareFloats(a.z(), 0.0f, EPSILON))
		{
			return true;
		}

		return false;
	}

	Vector3D operator+(const Vector3D& a, const Vector3D& b)
	{
		return Vector3D((double)a.x() + b.x(), (double)a.y() + b.y(), (double)a.z() + b.z());
	}

	Vector3D operator-(const Vector3D& v)
	{
		return Vector3D(-v.x(), -v.y(), -v.z());
	}

	Vector3D operator-(const Vector3D& a, const Vector3D& b)
	{
		return Vector3D(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
	}

	Vector3D operator*(const Vector3D& a, const float& k)
	{
		return Vector3D(a.x() * (double)k, a.y() * (double)k, a.z() * (double)k);
	}


	Vector3D operator*(const float& k, const Vector3D& a)
	{
		return Vector3D((double)k * a.x(), (double)k * a.y(), (double)k * a.z());
	}

	Vector3D operator/(const Vector3D& a, const float& k)
	{
		if (compareFloats(k, 0.0f, EPSILON))
		{
			return Vector3D();
		}

		return Vector3D(a.x() / (double)k, a.y() / (double)k, a.z() / (double)k);
	}

	//a dot b = axbx + ayby + azbz
	float dotProduct(const Vector3D& a, const Vector3D& b)
	{
		return (double)a.x() * b.x() + (double)a.y() * b.y() + (double)a.z() * b.z();
	}

	//a x b = (aybz - azby, azbx - axbz, axby - aybx)
	Vector3D crossProduct(const Vector3D& a, const Vector3D& b)
	{
		return Vector3D((double)a.y() * b.z() - (double)a.z() * b.y(),
			(double)a.z() * b.x() - (double)a.x() * b.z(),
			(double)a.x() * b.y() - (double)a.y() * b.x());
	}

	//length(v) = sqrt(vx^2 + vy^2 + vz^2)
	float length(const Vector3D& v)
	{
		return sqrt((double)v.x() * v.x() + (double)v.y() * v.y() + (double)v.z() * v.z());
	}

	//norm(v) = v / length(v) == (vx / length(v), vy / length(v))
	Vector3D norm(const Vector3D& v)
	{
		//v is the zero vector
		if (zeroVector(v))
		{
			return v;
		}

		double mag = length(v);

		return Vector3D(v.x() / mag, v.y() / mag, v.z() / mag);
	}

	//v = (r, theta, z)
	//x = rcos(theta)
	//y = rsin(theta)
	//z = z
	Vector3D CylindricalToCartesian(const Vector3D& v)
	{
		double angle = v.y() * PI / 180.0;

		return Vector3D(v.x() * cos(angle), v.x() * sin(angle), v.z());
	}

	//v = (x, y ,z)
	//r = sqrt(vx^2 + vy^2 + vz^2)
	//theta = arctan(y / x)
	//z = z
	Vector3D CartesianToCylindrical(const Vector3D& v)
	{
		if (compareFloats(v.x(), 0.0f, EPSILON))
		{
			return v;
		}

		double theta{ atan2(v.y(), v.x()) * 180.0 / PI };
		return Vector3D(length(v), theta, v.z());
	}

	// v = (pho, phi, theta)
	//x = pho * sin(phi) * cos(theta)
	//y = pho * sin(phi) * sin(theta)
	//z = pho * cos(theta);
	Vector3D SphericalToCartesian(const Vector3D& v)
	{
		double phi{ v.y() * PI / 180.0 };
		double theta{ v.z() * PI / 180.0 };

		return Vector3D(v.x() * sin(phi) * cos(theta), v.x() * sin(phi) * sin(theta), v.x() * cos(theta));
	}

	//v = (x, y ,z)
	//pho = sqrt(vx^2 + vy^2 + vz^2)
	//phi = arcos(z / pho)
	//theta = arctan(y / x)
	Vector3D CartesianToSpherical(const Vector3D& v)
	{
		if (compareFloats(v.x(), 0.0f, EPSILON) || zeroVector(v))
		{
			return v;
		}

		double pho{ length(v) };
		double phi{ acos(v.z() / pho) * 180.0 / PI };
		double theta{ atan2(v.y(), v.x()) * 180.0 / PI };

		return Vector3D(pho, phi, theta);
	}


	//Projb(a) = (a dot b)b
	//normalize b before projecting
	Vector3D Projection(const Vector3D& a, const Vector3D& b)
	{
		Vector3D normB(norm(b));
		return Vector3D(dotProduct(a, normB) * normB);
	}


#if defined(_DEBUG)
	void print(const Vector3D& v)
	{
		std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
	}
#endif







	//--------------------------------------------------------------------------------------
	//Constructors

	Vector4D::Vector4D() : m_x{ 0.0f }, m_y{ 0.0f }, m_z{ 0.0f }, m_w{ 0.0f }
	{}

	Vector4D::Vector4D(float x, float y, float z, float w) : m_x{ x }, m_y{ y }, m_z{ z }, m_w{ w }
	{}

	//--------------------------------------------------------------------------------------

	//--------------------------------------------------------------------------------------
	//Getters and Setters

	float Vector4D::x() const
	{
		return m_x;
	}

	float Vector4D::y() const
	{
		return m_y;
	}

	float Vector4D::z() const
	{
		return m_z;
	}

	float Vector4D::w() const
	{
		return m_w;
	}

	void Vector4D::setX(float x)
	{
		m_x = x;
	}

	void Vector4D::setY(float y)
	{
		m_y = y;
	}

	void Vector4D::setZ(float z)
	{
		m_z = z;
	}

	void Vector4D::setW(float w)
	{
		m_w = w;
	}

	//--------------------------------------------------------------------------------------


	//--------------------------------------------------------------------------------------
	//Memeber functions

	Vector4D& Vector4D::operator+=(const Vector4D& b)
	{
		this->m_x += (double)b.m_x;
		this->m_y += (double)b.m_y;
		this->m_z += (double)b.m_z;
		this->m_w += (double)b.m_w;

		return *this;
	}

	Vector4D& Vector4D::operator-=(const Vector4D& b)
	{
		this->m_x -= (double)b.m_x;
		this->m_y -= (double)b.m_y;
		this->m_z -= (double)b.m_z;
		this->m_w -= (double)b.m_w;

		return *this;
	}

	Vector4D& Vector4D::operator*=(const float& k)
	{
		this->m_x *= (double)k;
		this->m_y *= (double)k;
		this->m_z *= (double)k;
		this->m_w *= (double)k;

		return *this;
	}

	Vector4D& Vector4D::operator/=(const float& k)
	{
		if (compareFloats(k, 0.0f, EPSILON))
		{
			return *this;
		}

		this->m_x /= (double)k;
		this->m_y /= (double)k;
		this->m_z /= (double)k;
		this->m_w /= (double)k;

		return *this;
	}

	//-------------------------------------------------------------------------------------

	//-------------------------------------------------------------------------------------
	//Non-member functions

	bool zeroVector(const Vector4D& a)
	{
		if (compareFloats(a.x(), 0.0f, EPSILON) && compareFloats(a.y(), 0.0f, EPSILON) &&
			compareFloats(a.z(), 0.0f, EPSILON) && compareFloats(a.w(), 0.0f, EPSILON))
		{
			return true;
		}

		return false;
	}

	Vector4D operator+(const Vector4D& a, const Vector4D& b)
	{
		return Vector4D((double)a.x() + b.x(), (double)a.y() + b.y(), (double)a.z() + b.z(), (double)a.w() + b.w());
	}

	Vector4D operator-(const Vector4D& v)
	{
		return Vector4D(-v.x(), -v.y(), -v.z(), -v.w());
	}

	Vector4D operator-(const Vector4D& a, const Vector4D& b)
	{
		return Vector4D((double)a.x() - b.x(), (double)a.y() - b.y(), (double)a.z() - b.z(), (double)a.w() - b.w());
	}

	Vector4D operator*(const Vector4D& a, const float& k)
	{
		return Vector4D(a.x() * (double)k, a.y() * (double)k, a.z() * (double)k, a.w() * (double)k);
	}


	Vector4D operator*(const float& k, const Vector4D& a)
	{
		return Vector4D((double)k * a.x(), (double)k * a.y(), (double)k * a.z(), (double)k * a.w());
	}

	Vector4D operator/(const Vector4D& a, const float& k)
	{
		if (compareFloats(k, 0.0f, EPSILON))
		{
			return Vector4D();
		}

		return Vector4D(a.x() / (double)k, a.y() / (double)k, a.z() / (double)k, a.w() / (double)k);
	}

	//a dot b = axbx + ayby + azbz + awbw
	float dotProduct(const Vector4D& a, const Vector4D& b)
	{
		return (double)a.x() * b.x() + (double)a.y() * b.y() + (double)a.z() * b.z() + (double)a.w() * b.w();
	}

	//length(v) = sqrt(vx^2 + vy^2 + vz^2 + vw^2)
	float length(const Vector4D& v)
	{
		return sqrt((double)v.x() * v.x() + (double)v.y() * v.y() + (double)v.z() * v.z() + (double)v.w() * v.w());
	}

	//norm(v) = v / length(v) == (vx / length(v), vy / length(v))
	Vector4D norm(const Vector4D& v)
	{
		//v is the zero vector
		if (zeroVector(v))
		{
			return v;
		}

		double mag = length(v);

		return Vector4D(v.x() / mag, v.y() / mag, v.z() / mag, v.w() / mag);
	}

	//Projb(a) = (a dot b)b
	//normalize b before projecting
	Vector4D Projection(const Vector4D& a, const Vector4D& b)
	{
		Vector4D normB(norm(b));
		return Vector4D(dotProduct(a, normB) * normB);
	}


#if defined(_DEBUG)
	void print(const Vector4D& v)
	{
		std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.w() << ")";
	}
#endif








	Matrix4x4::Matrix4x4()
	{
		//1st row
		m_mat[0][0] = 1.0f;
		m_mat[0][1] = 0.0f;
		m_mat[0][2] = 0.0f;
		m_mat[0][3] = 0.0f;

		//2nd
		m_mat[1][0] = 0.0f;
		m_mat[1][1] = 1.0f;
		m_mat[1][2] = 0.0f;
		m_mat[1][3] = 0.0f;

		//3rd row
		m_mat[2][0] = 0.0f;
		m_mat[2][1] = 0.0f;
		m_mat[2][2] = 1.0f;
		m_mat[2][3] = 0.0f;

		//4th row
		m_mat[3][0] = 0.0f;
		m_mat[3][1] = 0.0f;
		m_mat[3][2] = 0.0f;
		m_mat[3][3] = 1.0f;
	}



	Matrix4x4::Matrix4x4(float a[][4])
	{
		//1st row
		m_mat[0][0] = a[0][0];
		m_mat[0][1] = a[0][1];
		m_mat[0][2] = a[0][2];
		m_mat[0][3] = a[0][3];

		//2nd
		m_mat[1][0] = a[1][0];
		m_mat[1][1] = a[1][1];
		m_mat[1][2] = a[1][2];
		m_mat[1][3] = a[1][3];

		//3rd row
		m_mat[2][0] = a[2][0];
		m_mat[2][1] = a[2][1];
		m_mat[2][2] = a[2][2];
		m_mat[2][3] = a[2][3];

		//4th row
		m_mat[3][0] = a[3][0];
		m_mat[3][1] = a[3][1];
		m_mat[3][2] = a[3][2];
		m_mat[3][3] = a[3][3];
	}

	float Matrix4x4::getElement(unsigned int row, unsigned int col) const
	{
		if (row > 3 || col > 3)
		{
			return m_mat[0][0];
		}
		else
		{
			return m_mat[row][col];
		}
	}


	const float& Matrix4x4::operator()(unsigned int row, unsigned int col) const
	{
		if (row > 3 || col > 3)
		{
			return m_mat[0][0];
		}
		else
		{
			return m_mat[row][col];
		}
	}

	float& Matrix4x4::operator()(unsigned int row, unsigned int col)
	{
		if (row > 3 || col > 3)
		{
			return m_mat[0][0];
		}
		else
		{
			return m_mat[row][col];
		}
	}

	void Matrix4x4::setElement(unsigned int row, unsigned int col, float val)
	{
		if (row > 3 || col > 3)
		{
			m_mat[0][0] = val;
		}
		else
		{
			m_mat[row][col] = val;
		}
	}

	//set to identity matrix by setting the diagonals to 1.0f and all other elements to 0.0f
	void Matrix4x4::setToIdentity()
	{
		//1st row
		m_mat[0][0] = 1.0f;
		m_mat[0][1] = 0.0f;
		m_mat[0][2] = 0.0f;
		m_mat[0][3] = 0.0f;

		//2nd
		m_mat[1][0] = 0.0f;
		m_mat[1][1] = 1.0f;
		m_mat[1][2] = 0.0f;
		m_mat[1][3] = 0.0f;

		//3rd row
		m_mat[2][0] = 0.0f;
		m_mat[2][1] = 0.0f;
		m_mat[2][2] = 1.0f;
		m_mat[2][3] = 0.0f;

		//4th row
		m_mat[3][0] = 0.0f;
		m_mat[3][1] = 0.0f;
		m_mat[3][2] = 0.0f;
		m_mat[3][3] = 1.0f;
	}

	//Is the identity matrix if the diagonals are equal to 1.0f and all other elements equals to 0.0f
	bool Matrix4x4::isIdentity()
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				if (i == j)
				{
					if (!compareFloats(m_mat[i][j], 1.0f, EPSILON))
					{
						return false;
					}
				}
				else
				{
					if (!compareFloats(m_mat[i][j], 0.0f, EPSILON))
					{
						return false;
					}
				}
			}
		}
	}

	Matrix4x4& Matrix4x4::operator+=(const Matrix4x4& m)
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				this->m_mat[i][j] += (double)m.m_mat[i][j];
			}
		}

		return *this;
	}

	Matrix4x4& Matrix4x4::operator-=(const Matrix4x4& m)
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				this->m_mat[i][j] -= (double)m.m_mat[i][j];
			}
		}

		return *this;
	}

	Matrix4x4& Matrix4x4::operator*=(const float& k)
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				this->m_mat[i][j] *= (double)k;
			}
		}

		return *this;
	}

	Matrix4x4& Matrix4x4::operator*=(const Matrix4x4& m)
	{
		Matrix4x4 res;

		for (int i = 0; i < 4; ++i)
		{
			res.m_mat[i][0] = ((double)m_mat[i][0] * m.m_mat[0][0]) +
				((double)m_mat[i][1] * m.m_mat[1][0]) +
				((double)m_mat[i][2] * m.m_mat[2][0]) +
				((double)m_mat[i][3] * m.m_mat[3][0]);

			res.m_mat[i][1] = ((double)m_mat[i][0] * m.m_mat[0][1]) +
				((double)m_mat[i][1] * m.m_mat[1][1]) +
				((double)m_mat[i][2] * m.m_mat[2][1]) +
				((double)m_mat[i][3] * m.m_mat[3][1]);

			res.m_mat[i][2] = ((double)m_mat[i][0] * m.m_mat[0][2]) +
				((double)m_mat[i][1] * m.m_mat[1][2]) +
				((double)m_mat[i][2] * m.m_mat[2][2]) +
				((double)m_mat[i][3] * m.m_mat[3][2]);

			res.m_mat[i][3] = ((double)m_mat[i][0] * m.m_mat[0][3]) +
				((double)m_mat[i][1] * m.m_mat[1][3]) +
				((double)m_mat[i][2] * m.m_mat[2][3]) +
				((double)m_mat[i][3] * m.m_mat[3][3]);
		}

		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				m_mat[i][j] = res.m_mat[i][j];
			}
		}

		return *this;
	}

	Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		Matrix4x4 res;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				res.setElement(i, j, (double)m1.getElement(i, j) + m2.getElement(i, j));
			}
		}

		return res;
	}

	Matrix4x4 operator-(const Matrix4x4& m)
	{
		Matrix4x4 res;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				res(i, j) = -m(i, j);
			}
		}

		return res;
	}


	Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		Matrix4x4 res;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				res.setElement(i, j, (double)m1.getElement(i, j) - m2.getElement(i, j));
			}
		}

		return res;
	}

	Matrix4x4 operator*(const Matrix4x4& m, const float& k)
	{
		Matrix4x4 res;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				res.setElement(i, j, (double)m.getElement(i, j) * k);
			}
		}

		return res;
	}

	Matrix4x4 operator*(const float& k, const Matrix4x4& m)
	{
		Matrix4x4 res;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				res.setElement(i, j, (double)k * m.getElement(i, j));
			}
		}

		return res;
	}

	Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		Matrix4x4 res;

		for (int i = 0; i < 4; ++i)
		{
			res(i, 0) = ((double)m1(i, 0) * m2(0, 0)) +
				((double)m1(i, 1) * m2(1, 0)) +
				((double)m1(i, 2) * m2(2, 0)) +
				((double)m1(i, 3) * m2(3, 0));

			res(i, 1) = ((double)m1(i, 0) * m2(0, 1)) +
				((double)m1(i, 1) * m2(1, 1)) +
				((double)m1(i, 2) * m2(2, 1)) +
				((double)m1(i, 3) * m2(3, 1));

			res(i, 2) = ((double)m1(i, 0) * m2(0, 2)) +
				((double)m1(i, 1) * m2(1, 2)) +
				((double)m1(i, 2) * m2(2, 2)) +
				((double)m1(i, 3) * m2(3, 2));

			res(i, 3) = ((double)m1(i, 0) * m2(0, 3)) +
				((double)m1(i, 1) * m2(1, 3)) +
				((double)m1(i, 2) * m2(2, 3)) +
				((double)m1(i, 3) * m2(3, 3));
		}

		return res;
	}

	Vector4D operator*(const Matrix4x4& m, const Vector4D& v)
	{
		Vector4D res;

		res.setX((double)m(0, 0) * v.x() + (double)m(0, 1) * v.y() + (double)m(0, 2) * v.z() + (double)m(0, 3) * v.w());
		res.setY((double)m(1, 0) * v.x() + (double)m(1, 1) * v.y() + (double)m(1, 2) * v.z() + (double)m(1, 3) * v.w());
		res.setZ((double)m(2, 0) * v.x() + (double)m(2, 1) * v.y() + (double)m(2, 2) * v.z() + (double)m(2, 3) * v.w());
		res.setW((double)m(3, 0) * v.x() + (double)m(3, 1) * v.y() + (double)m(3, 2) * v.z() + (double)m(3, 3) * v.w());

		return res;
	}

	Vector4D operator*(const Vector4D& v, const Matrix4x4& m)
	{
		Vector4D res;

		res.setX((double)v.x() * m(0, 0) + (double)v.y() * m(1, 0) + (double)v.z() * m(2, 0) + (double)v.w() * m(3, 0));
		res.setY((double)v.x() * m(0, 1) + (double)v.y() * m(1, 1) + (double)v.z() * m(2, 1) + (double)v.w() * m(3, 1));
		res.setZ((double)v.x() * m(0, 2) + (double)v.y() * m(1, 2) + (double)v.z() * m(2, 2) + (double)v.w() * m(3, 2));
		res.setW((double)v.x() * m(0, 3) + (double)v.y() * m(1, 3) + (double)v.z() * m(2, 3) + (double)v.w() * m(3, 3));

		return res;
	}

	//make the rows into cols
	Matrix4x4 transpose(const Matrix4x4& m)
	{
		Matrix4x4 res;

		//1st col = 1st row
		res(0, 0) = m(0, 0);
		res(1, 0) = m(0, 1);
		res(2, 0) = m(0, 2);
		res(3, 0) = m(0, 3);

		//2nd col = 2nd row
		res(0, 1) = m(1, 0);
		res(1, 1) = m(1, 1);
		res(2, 1) = m(1, 2);
		res(3, 1) = m(1, 3);

		//3rd col = 3rd row
		res(0, 2) = m(2, 0);
		res(1, 2) = m(2, 1);
		res(2, 2) = m(2, 2);
		res(3, 2) = m(2, 3);

		//4th col = 4th row
		res(0, 3) = m(3, 0);
		res(1, 3) = m(3, 1);
		res(2, 3) = m(3, 2);
		res(3, 3) = m(3, 3);

		return res;
	}

	//1 0 0 0
	//0 1 0 0
	//0 0 1 0
	//x y z 1
	Matrix4x4 translate(const Matrix4x4& cm, float x, float y, float z)
	{
		Matrix4x4 t;
		t(3, 0) = x;
		t(3, 1) = y;
		t(3, 2) = z;

		return cm * t;
	}

	//x 0 0 0
	//0 y 0 0
	//0 0 z 0
	//0 0 0 1
	Matrix4x4 scale(const Matrix4x4& cm, float x, float y, float z)
	{
		Matrix4x4 s;
		s(0, 0) = x;
		s(1, 1) = y;
		s(2, 2) = z;

		return cm * s;
	}

	//c + (1 - c)x^2	(1 - c)xy + sz	(1 - c)xz - sy	0
	//(1 - c)xy - sz	c + (1 - c)y^2	(1 - c)yz + sx	0
	//(1 - c)xz + sy	(1 - c)yz - sx	c + (1 - c)z^2	0
	//0					0				0				1
	//c = cos(angle)
	//s = sin(angle)
	Matrix4x4 rotate(const Matrix4x4& cm, float angle, float x, float y, float z)
	{
		double c = cos(angle * PI / 180.0);
		double s = sin(angle * PI / 180.0);


		Matrix4x4 r;

		//1st row
		r(0, 0) = c + (1.0 - c) * ((double)x * x);
		r(0, 1) = (1.0 - c) * ((double)x * y) + (s * z);
		r(0, 2) = (1.0 - c) * ((double)x * z) - (s * y);

		//2nd row
		r(1, 0) = (1.0 - c) * ((double)x * y) - (s * z);
		r(1, 1) = c + (1.0 - c) * ((double)y * y);
		r(1, 2) = (1.0 - c) * ((double)y * z) + (s * x);

		//3rd row
		r(2, 0) = (1.0 - c) * ((double)x * z) + (s * y);
		r(2, 1) = (1.0 - c) * ((double)y * z) - (s * x);
		r(2, 2) = c + (1.0 - c) * ((double)z * z);

		return cm * r;
	}

	double det(const Matrix4x4& m)
	{
		//m00m11(m22m33 - m23m32)
		double c1 = (double)m(0, 0) * m(1, 1) * m(2, 2) * m(3, 3) - (double)m(0, 0) * m(1, 1) * m(2, 3) * m(3, 2);

		//m00m12(m23m31 - m21m33)
		double c2 = (double)m(0, 0) * m(1, 2) * m(2, 3) * m(3, 1) - (double)m(0, 0) * m(1, 2) * m(2, 1) * m(3, 3);

		//m00m13(m21m32 - m22m31)
		double c3 = (double)m(0, 0) * m(1, 3) * m(2, 1) * m(3, 2) - (double)m(0, 0) * m(1, 3) * m(2, 2) * m(3, 1);

		//m01m10(m22m33 - m23m32)
		double c4 = (double)m(0, 1) * m(1, 0) * m(2, 2) * m(3, 3) - (double)m(0, 1) * m(1, 0) * m(2, 3) * m(3, 2);

		//m01m12(m23m30 - m20m33)
		double c5 = (double)m(0, 1) * m(1, 2) * m(2, 3) * m(3, 0) - (double)m(0, 1) * m(1, 2) * m(2, 0) * m(3, 3);

		//m01m13(m20m32 - m22m30)
		double c6 = (double)m(0, 1) * m(1, 3) * m(2, 0) * m(3, 2) - (double)m(0, 1) * m(1, 3) * m(2, 2) * m(3, 0);

		//m02m10(m21m33 - m23m31)
		double c7 = (double)m(0, 2) * m(1, 0) * m(2, 1) * m(3, 3) - (double)m(0, 2) * m(1, 0) * m(2, 3) * m(3, 1);

		//m02m11(m23m30 - m20m33)
		double c8 = (double)m(0, 2) * m(1, 1) * m(2, 3) * m(3, 0) - (double)m(0, 2) * m(1, 1) * m(2, 0) * m(3, 3);

		//m02m13(m20m31 - m21m30)
		double c9 = (double)m(0, 2) * m(1, 3) * m(2, 0) * m(3, 1) - (double)m(0, 2) * m(1, 3) * m(2, 1) * m(3, 0);

		//m03m10(m21m32 - m22m21)
		double c10 = (double)m(0, 3) * m(1, 0) * m(2, 1) * m(3, 2) - (double)m(0, 3) * m(1, 0) * m(2, 2) * m(3, 1);

		//m03m11(m22m30 - m20m32)
		double c11 = (double)m(0, 3) * m(1, 1) * m(2, 2) * m(3, 0) - (double)m(0, 3) * m(1, 1) * m(2, 0) * m(3, 2);

		//m03m12(m20m31 - m21m30)
		double c12 = (double)m(0, 3) * m(1, 2) * m(2, 0) * m(3, 1) - (double)m(0, 3) * m(1, 2) * m(2, 1) * m(3, 0);

		return (c1 + c2 + c3) - (c4 + c5 + c6) + (c7 + c8 + c9) - (c10 + c11 + c12);
	}

	//cij = (-1)^i + j * det of minor(i, j);
	double cofactor(const Matrix4x4& m, unsigned int row, unsigned int col)
	{
		double tempMat[3][3];
		int tr{ 0 };
		int tc{ 0 };

		//minor(i, j)
		for (int i = 0; i < 4; ++i)
		{
			if (i == row)
				continue;

			for (int j = 0; j < 4; ++j)
			{
				if (j == col)
					continue;

				tempMat[tr][tc] = m(i, j);
				++tc;

			}
			tc = 0;
			++tr;
		}

		//determinant of minor(i, j)
		double det3x3 = (tempMat[0][0] * tempMat[1][1] * tempMat[2][2]) + (tempMat[0][1] * tempMat[1][2] * tempMat[2][0]) +
			(tempMat[0][2] * tempMat[1][0] * tempMat[2][1]) - (tempMat[0][2] * tempMat[1][1] * tempMat[2][0]) -
			(tempMat[0][1] * tempMat[1][0] * tempMat[2][2]) - (tempMat[0][0] * tempMat[1][2] * tempMat[2][1]);

		return pow(-1, row + col) * det3x3;
	}

	//Cofactor of each ijth position put into matrix cA.
	//Adjoint is the tranposed matrix of cA.
	Matrix4x4 adjoint(const Matrix4x4& m)
	{
		Matrix4x4 cA;
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				cA(i, j) = cofactor(m, i, j);
			}
		}

		return transpose(cA);
	}

	//Inverse of m = adjoint of m / det of m
	Matrix4x4 inverse(const Matrix4x4& m)
	{
		double determinant = det(m);
		if (compareDoubles(determinant, 0.0, EPSILON))
			return Matrix4x4();

		return adjoint(m) * (1.0 / determinant);
	}


#if defined(_DEBUG)
	void print(const Matrix4x4& m)
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				std::cout << m(i, j) << " ";
			}

			std::cout << std::endl;
		}
	}
#endif




	Quaternion::Quaternion() : m_scalar{ 1.0f }, m_x{ 0.0f }, m_y{ 0.0f }, m_z{ 0.0f }
	{
	}

	Quaternion::Quaternion(float scalar, float x, float y, float z) : m_scalar{ scalar }, m_x{ x }, m_y{ y }, m_z{ z }
	{
	}

	Quaternion::Quaternion(float scalar, const Vector3D& v) : m_scalar{ scalar }, m_x{ v.x() }, m_y{ v.y() }, m_z{ v.z() }
	{
	}

	Quaternion::Quaternion(const Vector4D& v) : m_scalar{ v.x() }, m_x{ v.y() }, m_y{ v.z() }, m_z{ v.w() }
	{
	}

	float& Quaternion::scalar()
	{
		return m_scalar;
	}

	const float& Quaternion::scalar() const
	{
		return m_scalar;
	}

	float& Quaternion::x()
	{
		return m_x;
	}

	const float& Quaternion::x() const
	{
		return m_x;
	}

	float& Quaternion::y()
	{
		return m_y;
	}

	const float& Quaternion::y() const
	{
		return m_y;
	}

	float& Quaternion::z()
	{
		return m_z;
	}

	const float& Quaternion::z() const
	{
		return m_z;
	}

	Vector3D Quaternion::vector()
	{
		return Vector3D(m_x, m_y, m_z);
	}

	Quaternion& Quaternion::operator+=(const Quaternion& q)
	{
		this->m_scalar += (double)q.m_scalar;
		this->m_x += (double)q.m_x;
		this->m_y += (double)q.m_y;
		this->m_z += (double)q.m_z;

		return *this;
	}

	Quaternion& Quaternion::operator-=(const Quaternion& q)
	{
		this->m_scalar -= (double)q.m_scalar;
		this->m_x -= (double)q.m_x;
		this->m_y -= (double)q.m_y;
		this->m_z -= (double)q.m_z;

		return *this;
	}


	Quaternion& Quaternion::operator*=(float k)
	{
		this->m_scalar *= (double)k;
		this->m_x *= (double)k;
		this->m_y *= (double)k;
		this->m_z *= (double)k;

		return *this;
	}

	Quaternion& Quaternion::operator*=(const Quaternion& q)
	{
		Vector3D thisVector(this->m_x, this->m_y, this->m_z);
		Vector3D qVector(q.m_x, q.m_y, q.m_z);

		double s{ (double)this->m_scalar * q.m_scalar };
		double dP{ dotProduct(thisVector, qVector) };
		double resultScalar{ s - dP };

		Vector3D a(this->m_scalar * qVector);
		Vector3D b(q.m_scalar * thisVector);
		Vector3D cP(crossProduct(thisVector, qVector));
		Vector3D resultVector(a + b + cP);

		this->m_scalar = resultScalar;
		this->m_x = resultVector.x();
		this->m_y = resultVector.y();
		this->m_z = resultVector.z();

		return *this;
	}


	Quaternion operator+(const Quaternion& q1, const Quaternion& q2)
	{
		return Quaternion((double)q1.scalar() + q2.scalar(), (double)q1.x() + q2.x(), (double)q1.y() + q2.y(), (double)q1.z() + q2.z());
	}

	Quaternion operator-(const Quaternion& q)
	{
		return Quaternion(-q.scalar(), -q.x(), -q.y(), -q.z());
	}

	Quaternion operator-(const Quaternion& q1, const Quaternion& q2)
	{
		return Quaternion((double)q1.scalar() - q2.scalar(), (double)q1.x() - q2.x(), (double)q1.y() - q2.y(), (double)q1.z() - q2.z());
	}

	Quaternion operator*(float k, const Quaternion& q)
	{
		return Quaternion((double)k * q.scalar(), (double)k * q.x(), (double)k * q.y(), (double)k * q.z());
	}

	Quaternion operator*(const Quaternion& q, float k)
	{
		return Quaternion(q.scalar() * (double)k, q.x() * (double)k, q.y() * (double)k, q.z() * (double)k);
	}

	//scalar part = q1scalar * q2scalar - q1Vector dot q2Vector
	//vector part = q1Scalar * q2Vector + q2Scalar * q1Vector + q1Vector cross q2Vector
	Quaternion operator*(const Quaternion& q1, const Quaternion& q2)
	{
		Vector3D q1Vector(q1.x(), q1.y(), q1.z());
		Vector3D q2Vector(q2.x(), q2.y(), q2.z());

		double s{ (double)q1.scalar() * q2.scalar() };
		double dP{ dotProduct(q1Vector, q2Vector) };
		double resultScalar{ s - dP };

		Vector3D a(q1.scalar() * q2Vector);
		Vector3D b(q2.scalar() * q1Vector);
		Vector3D cP(crossProduct(q1Vector, q2Vector));
		Vector3D resultVector(a + b + cP);

		return Quaternion(resultScalar, resultVector);
	}

	//zero quaternion = (0, 0, 0, 0)
	bool isZeroQuaternion(const Quaternion& q)
	{
		return compareFloats(q.scalar(), 0.0f, EPSILON) && compareFloats(q.x(), 0.0f, EPSILON) &&
			compareFloats(q.y(), 0.0f, EPSILON) && compareFloats(q.z(), 0.0f, EPSILON);
	}

	//identity quaternion = (1, 0, 0, 0)
	bool isIdentity(const Quaternion& q)
	{
		return compareFloats(q.scalar(), 1.0f, EPSILON) && compareFloats(q.x(), 0.0f, EPSILON) &&
			compareFloats(q.y(), 0.0f, EPSILON) && compareFloats(q.z(), 0.0f, EPSILON);
	}

	//conjugate of a quaternion is the quaternion with its vector part negated
	Quaternion conjugate(const Quaternion& q)
	{
		return Quaternion(q.scalar(), -q.x(), -q.y(), -q.z());
	}

	//length of a quaternion = sqrt(scalar^2 + x^2 + y^2 + z^2)
	float length(const Quaternion& q)
	{
		return sqrt((double)q.scalar() * q.scalar() + (double)q.x() * q.x() + (double)q.y() * q.y() + (double)q.z() * q.z());
	}

	//to normalize a quaternion you do q / |q|
	Quaternion normalize(const Quaternion& q)
	{
		if (isZeroQuaternion(q))
			return q;

		double d{ length(q) };

		return Quaternion(q.scalar() / d, q.x() / d, q.y() / d, q.z() / d);
	}

	//inverse = conjugate of q / |q|^2
	Quaternion inverse(const Quaternion& q)
	{
		if (isZeroQuaternion(q))
			return q;

		Quaternion conjugateOfQ(conjugate(q));

		double d{ length(q) };
		d *= d;

		return Quaternion(conjugateOfQ.scalar() / d, conjugateOfQ.x() / d, conjugateOfQ.y() / d, conjugateOfQ.z() / d);
	}

	//A roatation quaternion is a quaternion where the
	//scalar part = cos(theta / 2)
	//vector part = sin(theta / 2) * axis
	//the axis needs to be normalized
	Quaternion rotationQuaternion(float angle, float x, float y, float z)
	{
		double ang{ angle / 2.0 };
		double c{ cos(ang * PI / 180.0) };
		double s{ sin(ang * PI / 180.0) };

		Vector3D axis(x, y, z);
		axis = norm(axis);

		return Quaternion(c, s * axis.x(), s * axis.y(), s * axis.z());
	}

	//A roatation quaternion is a quaternion where the
	//scalar part = cos(theta / 2)
	//vector part = sin(theta / 2) * axis
	//the axis needs to be normalized
	Quaternion rotationQuaternion(float angle, const Vector3D& axis)
	{
		double ang{ angle / 2.0 };
		double c{ cos(ang * PI / 180.0) };
		double s{ sin(ang * PI / 180.0) };

		Vector3D axisN(norm(axis));

		return Quaternion(c, s * axisN.x(), s * axisN.y(), s * axisN.z());
	}

	//A roatation quaternion is a quaternion where the
	//scalar part = cos(theta / 2)
	//vector part = sin(theta / 2) * axis
	//the axis needs to be normalized
	Quaternion rotationQuaternion(const Vector4D& angAxis)
	{
		double angle{ angAxis.x() / 2.0 };
		double c{ cos(angle * PI / 180.0) };
		double s{ sin(angle * PI / 180.0) };

		Vector3D axis(angAxis.y(), angAxis.z(), angAxis.w());
		axis = norm(axis);

		return Quaternion(c, s * axis.x(), s * axis.y(), s * axis.z());
	}

	//1 - 2q3^2 - 2q4^2		2q2q3 - 2q1q4		2q2q4 + 2q1q3		0
	//2q2q3 + 2q1q4			1 - 2q2^2 - 2q4^2	2q3q4 - 2q1q2		0
	//2q2q4 - 2q1q3			2q3q4 + 2q1q2		1 - 2q2^2 - 2q3^2	0
	//0						0					0					1
	//q1 = scalar
	//q2 = x
	//q3 = y
	//q4 = z
	Matrix4x4 quaternionRotationMatrixCol(const Quaternion& q)
	{
		float colMat[4][4] = {};

		colMat[0][0] = 1.0 - 2.0 * q.y() * q.y() - 2.0 * q.z() * q.z();
		colMat[0][1] = 2.0 * q.x() * q.y() - 2.0 * q.scalar() * q.z();
		colMat[0][2] = 2.0 * q.x() * q.z() + 2.0 * q.scalar() * q.y();
		colMat[0][3] = 0.0f;

		colMat[1][0] = 2.0 * q.x() * q.y() + 2.0 * q.scalar() * q.z();
		colMat[1][1] = 1.0 - 2.0 * q.x() * q.x() - 2.0 * q.z() * q.z();
		colMat[1][2] = 2.0 * q.y() * q.z() - 2.0 * q.scalar() * q.x();
		colMat[1][3] = 0.0f;

		colMat[2][0] = 2.0 * q.x() * q.z() - 2.0 * q.scalar() * q.y();
		colMat[2][1] = 2.0 * q.y() * q.z() + 2.0 * q.scalar() * q.x();
		colMat[2][2] = 1.0 - 2.0 * q.x() * q.x() - 2.0 * q.y() * q.y();
		colMat[2][3] = 0.0f;

		colMat[3][0] = 0.0f;
		colMat[3][1] = 0.0f;
		colMat[3][2] = 0.0f;
		colMat[3][3] = 1.0f;

		return Matrix4x4(colMat);
	}


	//1 - 2q3^2 - 2q4^2		2q2q3 + 2q1q4		2q2q4 - 2q1q3		0
	//2q2q3 - 2q1q4			1 - 2q2^2 - 2q4^2	2q3q4 + 2q1q2		0
	//2q2q4 + 2q1q3			2q3q4 - 2q1q2		1 - 2q2^2 - 2q3^2	0
	//0						0					0					1
	//q1 = scalar
	//q2 = x
	//q3 = y
	//q4 = z
	Matrix4x4 quaternionRotationMatrixRow(const Quaternion& q)
	{
		float rowMat[4][4] = {};

		rowMat[0][0] = 1.0 - 2.0 * q.y() * q.y() - 2.0 * q.z() * q.z();
		rowMat[0][1] = 2.0 * q.x() * q.y() + 2.0 * q.scalar() * q.z();
		rowMat[0][2] = 2.0 * q.x() * q.z() - 2.0 * q.scalar() * q.y();
		rowMat[0][3] = 0.0f;

		rowMat[1][0] = 2.0 * q.x() * q.y() - 2.0 * q.scalar() * q.z();
		rowMat[1][1] = 1.0 - 2.0 * q.x() * q.x() - 2.0 * q.z() * q.z();
		rowMat[1][2] = 2.0 * q.y() * q.z() + 2.0 * q.scalar() * q.x();
		rowMat[1][3] = 0.0f;

		rowMat[2][0] = 2.0 * q.x() * q.z() + 2.0 * q.scalar() * q.y();
		rowMat[2][1] = 2.0 * q.y() * q.z() - 2.0 * q.scalar() * q.x();
		rowMat[2][2] = 1.0 - 2.0 * q.x() * q.x() - 2.0 * q.y() * q.y();
		rowMat[2][3] = 0.0f;

		rowMat[3][0] = 0.0f;
		rowMat[3][1] = 0.0f;
		rowMat[3][2] = 0.0f;
		rowMat[3][3] = 1.0f;

		return Matrix4x4(rowMat);
	}

#if defined(_DEBUG)
	void print(const Quaternion& q)
	{
		std::cout << "(" << q.scalar() << ", " << q.x() << ", " << q.y() << ", " << q.z();
	}
#endif

}
