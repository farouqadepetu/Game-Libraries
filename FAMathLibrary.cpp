#include <iostream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <stdlib.h>
#include "FAMathLibrary.h"
#include "FAUtil.h"

#define PI 3.14159265
#define EPSILON 1e-7

namespace FAMath
{

//VECTOR2
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Constructors for Vector2 class

	//Creates a new 2D vector/point with the components initialized to 0.0
	Vector2::Vector2() : m_x{ 0.0 }, m_y{ 0.0 }
	{}

	//Creates a new 2D vector/point with the components initialized to the arguments
	Vector2::Vector2(double x, double y) : m_x{ x }, m_y{ y }
	{}

	//Creates a new 2D vector/point using v's x and y coordinate
	Vector2::Vector2(const Vector3& v) : m_x{ v.x() }, m_y{ v.y() }
	{}

	//Creates a new 2D vector/point using v's x and y coordinate
	Vector2::Vector2(const Vector4& v) : m_x{ v.x() }, m_y{ v.y() }
	{}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Getters for Vector2 class
	
	//Returns the value of the x-coordinate
	double Vector2::x() const
	{
		return m_x;
	}

	//Returns the value of the y-coordinate
	double Vector2::y() const
	{
		return m_y;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Setters for Vector2 class

	//Sets the x-coordinate of the 2D vector/point
	void Vector2::setX(double x)
	{
		m_x = x;
	}

	//Sets the y-coordinate of the 2D vector/point
	void Vector2::setY(double y)
	{
		m_y = y;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Member Functions

	//Returns true if all the components of the 2D vector equals to zero, false otherwise.
	bool Vector2::isZeroVector() const
	{
		return FAUtil::compareDoubles(m_x, 0.0, EPSILON) && FAUtil::compareFloats(m_y, 0.0, EPSILON);
	}

	//Adds two 2D vectors
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector2 object + Vector2 object b
	Vector2& Vector2::operator+=(const Vector2& b)
	{
		this->m_x += b.m_x;
		this->m_y += b.m_y;

		return *this;
	}

	//Subtracts two 2D vectors
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector2 object - Vector2 object b
	Vector2& Vector2::operator-=(const Vector2& b)
	{
		this->m_x -= b.m_x;
		this->m_y -= b.m_y;

		return *this;
	}

	//Multiply each component in the 2D vector by scalar
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector2 object * scalar
	Vector2& Vector2::operator*=(const double& scalar)
	{
		this->m_x *= scalar;
		this->m_y *= scalar;

		return *this;
	}

	//Divide each component in the 2D vector by scalar
	//Stores result in current vector object
	//Throws a std::invlaid_argument if scalar is zero
	//Returns a reference to the current vector object
	//That has the result of the current Vector2 object / scalar
	Vector2& Vector2::operator/=(const double& scalar)
	{
		//if scalar is 0.0f throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		this->m_x /= scalar;
		this->m_y /= scalar;

		return *this;
	}

	//Stores the x and y values of Vector3s v into the x and y values of this Vector2
	void Vector2::operator=(const Vector3& v)
	{
		this->m_x = v.x();
		this->m_y = v.y();
	}

	//Stores the x and y values of Vector4s v into the x and y values of this Vector2
	void Vector2::operator=(const Vector4& v)
	{
		this->m_x = v.x();
		this->m_y = v.y();
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Friend Functions for Vector2 class

	//Adds two 2D vectors
	//Vector addition for 2D is (ax + bx, ay + by)
	//Returns a new Vector2 object that has the result of a + b
	Vector2 operator+(const Vector2& a, const Vector2& b)
	{
		return Vector2(a.m_x + b.m_x, a.m_y + b.m_y);
	}

	//Subtracts two 2D vectors
	//Vector addition for 2D is (ax - bx, ay - by)
	//Returns a new Vector2 object that has the result of a - b
	Vector2 operator-(const Vector2& a, const Vector2& b)
	{
		return Vector2(a.m_x - b.m_x, a.m_y - b.m_y);
	}

	//Multiplies each component in the 2D vector by scalar
	//Vector multiplicaton by a scalar for 2D is (vx * scalar , vy * scalar)
	//Returns a new Vector2 object that has the result of v * scalar
	Vector2 operator*(const Vector2& v, const double& scalar)
	{
		return Vector2(v.m_x * scalar, v.m_y * scalar);
	}

	//Multiplies each component in the 2D vector by scalar
	//Vector multiplicaton by a scalar for 2D is (vx * scalar , vy * scalar)
	//Returns a new Vector2 object that has the result of v * scalar
	Vector2 operator*(const double& scalar, const Vector2& v)
	{
		return Vector2(v.m_x * scalar, v.m_y * scalar);
	}

	//Divides each componet in the 2D vector by a scalar
	//Vector division by a scalar for 2D is (vx / scalar, vy / scalar)
	//Throws a std::invalid_argument if scalar is zero
	//Returns a new Vector2 object with the result of v / scalar
	Vector2 operator/(const Vector2& v, const double& scalar)
	{
		//if scalar is 0.0f throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		return Vector2(v.m_x / scalar, v.m_y / scalar);
	}

	//Returns true if a is equal to b, false otherwise
	bool operator==(const Vector2& a, const Vector2& b)
	{
		return FAUtil::compareDoubles(a.x(), b.x(), EPSILON) && FAUtil::compareDoubles(a.y(), b.y(), EPSILON);
	}

	//Returns false if a is equal to b, true otherwise
	bool operator!=(const Vector2& a, const Vector2& b)
	{
		return !(a == b);
	}

	//Returns the magnitude/length of a 2D vector
	//Length of a 2D vector is sqrt(vx * vx + vy * vy)
	double length(const Vector2& v)
	{
		return sqrt(v.m_x * v.m_x + v.m_y * v.m_y);
	}

	//Normalizes a vector
	//Returns a new Vector2 object that is a 2D unit vector
	//A unit vector is a vector that has a length of 1
	Vector2 normalize(const Vector2& v)
	{
		if (v.isZeroVector())
		{
			return v;
		}

		return v / length(v);
	}

	//Returns the distance between two 2D points
	//To get the distance between two 2D points:
	//Subtract the two 2D points and take the length of the resultant 2D point
	//length(a-b)
	//Or sqrt((ax - bx)^2 + (ay - by)^2)
	double distance(const Vector2& a, const Vector2& b)
	{
		Vector2 c{ a - b };

		return length(c);
	}

	//Does the dot product between two 2D vectors and returns the value
	//Dot Product between two 2D vectors: ax * bx + ay * by
	double dotProduct(const Vector2& a, const Vector2& b)
	{
		return a.m_x * b.m_x + a.m_y * b.m_y;
	}

	//Returns the angle between two vectors
	//Angle between two vectors = acos(a dot b)
	//a and b should be unit vectors before calling this function
	double angle(const Vector2& a, const Vector2& b)
	{
		double dP = dotProduct(a, b);
		if (dP < -1.0)
		{
			dP = -1.0;
		}
		else if (dP > 1.0)
		{
			dP = 1.0;
		}

		//radians to degrees = 180 degrees / PI rad
		return acos(dP) * 180.0f / PI;
	}


	void print(Vector2 v)
	{
		std::cout << "(" << v.m_x << ", " << v.m_y << ")" << std::endl;
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------







//VECTOR3
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Constructors for Vector3 class

	//Creates a new 3D vector/point with the components initialized to 0.0
	Vector3::Vector3() : m_x{ 0.0 }, m_y{ 0.0 }, m_z{ 0.0 }
	{}

	//Creates a new 3D vector/point with the components initialized to the arguments
	Vector3::Vector3(double x, double y, double z) : m_x{ x }, m_y{ y }, m_z{ z }
	{}

	//Creates a new 3D vector/point using v's x and y coordinate and sets z to 0.0f
	Vector3::Vector3(const Vector2& v) : m_x{ v.x() }, m_y{ v.y() }, m_z{ 0.0 }
	{}

	//Creates a new 3D vector/point using v's x and y coordinates and sets this Vector3 z component to the given z
	Vector3::Vector3(const Vector2& v, const double& z) : m_x{ v.x() }, m_y{ v.y() }, m_z{ z }
	{}

	//Creates a new 3D vector/point using v's x, y  and z coordinate
	Vector3::Vector3(const Vector4& v): m_x{ v.x() }, m_y{ v.y() }, m_z{ v.z() }
	{}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Getters for Vector3 class

	//Returns the value of the x-coordinate
	double Vector3::x() const
	{
		return m_x;
	}

	//Returns the value of the y-coordinate
	double Vector3::y() const
	{
		return m_y;
	}

	//Returns the value of the z-coordinate
	double Vector3::z() const
	{
		return m_z;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Setters for Vector3 class

	//Sets the x-coordinate of the 3D vector/point
	void Vector3::setX(double x)
	{
		m_x = x;
	}

	//Sets the y-coordinate of the 3D vector/point
	void Vector3::setY(double y)
	{
		m_y = y;
	}

	//Sets the z-coordinate of the 3D vector/point
	void Vector3::setZ(double z)
	{
		m_z = z;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Member Functions

	//Returns true if all the components of the 3D vector equals to zero, false otherwise.
	bool Vector3::isZeroVector() const
	{
		return FAUtil::compareDoubles(m_x, 0.0, EPSILON) && FAUtil::compareDoubles(m_y, 0.0, EPSILON) && FAUtil::compareDoubles(m_z, 0.0, EPSILON);
	}

	//Adds two 3D vectors
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector3 object + Vector3 object b
	Vector3& Vector3::operator+=(const Vector3& b)
	{
		this->m_x += b.m_x;
		this->m_y += b.m_y;
		this->m_z += b.m_z;

		return *this;
	}

	//Subtracts two 3D vectors
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector3 object - Vector3 object b
	//Time Complexity: O(1)
	//Space Complexity: O(1)
	Vector3& Vector3::operator-=(const Vector3& b)
	{
		this->m_x -= b.m_x;
		this->m_y -= b.m_y;
		this->m_z -= b.m_z;

		return *this;
	}

	//Multiply each component in the 3D vector by scalar
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector3 object * scalar
	//Time Complexity: O(1)
	//Space Complexity: O(1)
	Vector3& Vector3::operator*=(const double& scalar)
	{
		this->m_x *= scalar;
		this->m_y *= scalar;
		this->m_z *= scalar;

		return *this;
	}

	//Divide each component in the 3D vector by scalar
	//Stores result in current vector object
	//Throws a std::invlaid_argument if scalar is zero
	//Returns a reference to the current vector object
	//That has the result of the current Vector3 object / scalar
	Vector3& Vector3::operator/=(const double& scalar)
	{
		//if scalar is 0.0 throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		this->m_x /= scalar;
		this->m_y /= scalar;
		this->m_z /= scalar;

		return *this;
	}

	//Stores the xand y values of Vector2s v into the xand y values of this Vector3and sets z = 0.0f
	void Vector3::operator=(const Vector2& v)
	{
		this->m_x = v.x();
		this->m_y = v.y();
		this->m_z = 0.0f;
	}

	//Stores the x, y and z values of Vector4s v into the x, y and z values of this Vector3
	void Vector3::operator=(const Vector4& v)
	{
		this->m_x = v.x();
		this->m_y = v.y();
		this->m_z = v.z();
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Friend Functions for Vector3 class

	//Adds two 3D vectors
	//Vector addition for 3D is (ax + bx, ay + by, az + bz)
	//Returns a new Vector3 object that has the result of a + b
	Vector3 operator+(const Vector3& a, const Vector3& b)
	{
		return Vector3(a.m_x + b.m_x, a.m_y + b.m_y, a.m_z + b.m_z);
	}

	//Subtracts two 3D vectors
	//Vector addition for 3D is (ax - bx, ay - by, az - bz)
	//Returns a new Vector3 object that has the result of a - b
	Vector3 operator-(const Vector3& a, const Vector3& b)
	{
		return Vector3(a.m_x - b.m_x, a.m_y - b.m_y, a.m_z - b.m_z);
	}

	//Multiplies each component in the 3D vector by scalar
	//Vector multiplicaton by a scalar for 3D is (vx * scalar , vy * scalar, vz * scalar)
	//Returns a new Vector3 object that has the result of v * scalar
	Vector3 operator*(const Vector3& v, const double& scalar)
	{
		return Vector3(v.m_x * scalar, v.m_y * scalar, v.m_z * scalar);
	}

	//Multiplies each component in the 3D vector by scalar
	//Vector multiplicaton by a scalar for 3D is (vx * scalar , vy * scalar, vz * scalar)
	//Returns a new Vector3 object that has the result of v * scalar
	Vector3 operator*(const double& scalar, const Vector3& v)
	{
		return Vector3(v.m_x * scalar, v.m_y * scalar, v.m_z * scalar);
	}

	//Divides each componet in the 3D vector by a scalar
	//Vector division by a scalar for 3D is (vx / scalar, vy / scalar, vz / scalar)
	//Throws a std::invalid_argument if scalar is zero
	//Returns a new Vector3 object with the result of v / scalar
	Vector3 operator/(const Vector3& v, const double& scalar)
	{
		//if scalar is 0.0 throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		return Vector3(v.m_x / scalar, v.m_y / scalar, v.m_z / scalar);
	}

	//Returns true if a is equal to b, false otherwise
	bool operator==(const Vector3& a, const Vector3& b)
	{
		return FAUtil::compareDoubles(a.x(), b.x(), EPSILON) && FAUtil::compareDoubles(a.y(), b.y(), EPSILON) && FAUtil::compareDoubles(a.z(), b.z(), EPSILON);
	}

	//Returns false if a is equal to b, true otherwise
	bool operator!=(const Vector3& a, const Vector3& b)
	{
		return !(a == b);
	}

	//Returns the magnitude/length of a 3D vector
	//Length of a 3D vector is sqrt(vx * vx + vy * vy + vz * vz)
	//Time Complexity: O(1)
	//Space Complexity: O(1)
	double length(const Vector3& v)
	{
		return sqrt(v.m_x * v.m_x + v.m_y * v.m_y + v.m_z * v.m_z);
	}

	//Normalizes a vector
	//Returns a new Vector3 object that is a 3D unit vector
	//A unit vector is a vector that has a length of 1
	Vector3 normalize(const Vector3& v)
	{
		if (v.isZeroVector())
		{
			return v;
		}

		return v / length(v);
	}

	//Returns the distance between two 3D points
	//To get the distance between two 3D points:
	//Subtract the two 3D points and take the length of the resultant 3D point
	//length(a-b)
	//Or sqrt((ax - bx)^2 + (ay - by)^2 + (az - bz)^2)
	double distance(const Vector3& a, const Vector3& b)
	{
		Vector3 c{ a - b };

		return length(c);
	}

	//Does the dot product between two 3D vectors and returns the value
	//Dot Product between two 3D vectors: ax * bx + ay * by + az * bz
	double dotProduct(const Vector3& a, const Vector3& b)
	{
		return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
	}

	//Returns the angle between two vectors
	//Angle between two vectors = acos(a dot b)
	//a and b should be unit vectors before calling this function
	double angle(const Vector3& a, const Vector3& b)
	{
		double dP = dotProduct(a, b);
		if (dP < -1.0)
		{
			dP = -1.0;
		}
		else if (dP > 1.0)
		{
			dP = 1.0;
		}

		//radians to degrees = 180 degrees / PI rad
		return acos(dP) * 180.0 / PI;
	}

	//Does the cross product between two 3D vector and returns the vector produced
	//Cross Product between two 3D vectors is:
	//x: ay * bz - az * by
	//y: az * bz - ax * bz
	//z: ax * by - ay * bx
	Vector3 crossProduct(const Vector3& a, const Vector3& b)
	{
		return Vector3(a.y() * b.z() - a.z() * b.y(), a.z() * b.x() - a.x() * b.z(), a.x() * b.y() - a.y() * b.x());
	}

	void print(Vector3 v)
	{
		std::cout << std::setprecision(15) << "(" << v.m_x << ", " << v.m_y << ", " << v.m_z << ")" << std::endl;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------












//Vector4
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Constructors for Vector4 class

	//Creates a new 4D vector/point with the components initialized to 0.0.
	Vector4::Vector4() : m_x{ 0.0 }, m_y{ 0.0 }, m_z{ 0.0 }, m_w{ 0.0 }
	{}

	//Creates a new 4D vector/point with the components initialized to the arguments

	Vector4::Vector4(double x, double y, double z, double w) : m_x{ x }, m_y{ y }, m_z{ z }, m_w{ w }
	{}

	//Creates a new 4D vector/point using v's x and y coordinates and sets z and w to 0.0f

	Vector4::Vector4(const Vector2& v) : m_x{ v.x() }, m_y{ v.y() }, m_z{ 0.0 }, m_w{ 0.0 }
	{}

	//Creates a new 4D vector / point using v's x and y coordinates and the given z. Sets w to 0.0f
	Vector4::Vector4(const Vector2& v, const double& z) : m_x{ v.x() }, m_y{ v.y() }, m_z{ z }, m_w{ 0.0 }
	{}

	//Creates a new 4D vector / point using v's x and y coordinates and the given z and w values
	Vector4::Vector4(const Vector2& v, const double& z, const double& w): m_x { v.x() }, m_y{ v.y() }, m_z{ z }, m_w{ w }
	{}

	//Creates a new 4D vector/point using v's x, y and z coordinates and sets w to 0.0
	Vector4::Vector4(const Vector3& v) : m_x{ v.x() }, m_y{ v.y() }, m_z{ v.z() }, m_w{ 0.0 }
	{}

	//Creates a new 4D vector/point using v's x, y and z coordinates and the given w value
	Vector4::Vector4(const Vector3& v, const double& w) : m_x{ v.x() }, m_y{ v.y() }, m_z{ v.z() }, m_w{ w }
	{}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Getters for Vector4 class

	//Returns the value of the x-coordinate
	double Vector4::x() const
	{
		return m_x;
	}

	//Returns the value of the y-coordinate
	double Vector4::y() const
	{
		return m_y;
	}

	//Returns the value of the z-coordinate
	double Vector4::z() const
	{
		return m_z;
	}

	//Returns the value of the w-coordinate
	double Vector4::w() const
	{
		return m_w;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Setters for Vector4 class

	//Sets the x-coordinate of the 4D vector/point
	void Vector4::setX(double x)
	{
		m_x = x;
	}

	//Sets the y-coordinate of the 4D vector/point
	void Vector4::setY(double y)
	{
		m_y = y;
	}

	//Sets the z-coordinate of the 4D vector/point
	void Vector4::setZ(double z)
	{
		m_z = z;
	}

	//Sets the w-coordinate of the 4D vector/point)
	void Vector4::setW(double w)
	{
		m_w = w;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Member Functions
	//Returns true if all the components of the 4D vector equals to zero, false otherwise.
	bool Vector4::isZeroVector() const
	{
		return FAUtil::compareDoubles(m_x, 0.0, EPSILON) && FAUtil::compareDoubles(m_y, 0.0, EPSILON) && FAUtil::compareDoubles(m_z, 0.0, EPSILON)
			&& FAUtil::compareDoubles(m_w, 0.0, EPSILON);
	}

	//Adds two 4D vectors
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector4 object + Vector4 object b
	Vector4& Vector4::operator+=(const Vector4& b)
	{
		this->m_x += b.m_x;
		this->m_y += b.m_y;
		this->m_z += b.m_z;
		this->m_w += b.m_w;

		return *this;
	}

	//Subtracts two 4D vectors
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector4 object - Vector4 object b
	Vector4& Vector4::operator-=(const Vector4& b)
	{
		this->m_x -= b.m_x;
		this->m_y -= b.m_y;
		this->m_z -= b.m_z;
		this->m_w -= b.m_w;

		return *this;
	}

	//Multiply each component in the 4D vector by scalar
	//Stores result in current vector object
	//Returns a reference to the current vector object
	//That has the result of the current Vector4 object * scalar
	Vector4& Vector4::operator*=(const double& scalar)
	{
		this->m_x *= scalar;
		this->m_y *= scalar;
		this->m_z *= scalar;
		this->m_w *= scalar;

		return *this;
	}

	//Divide each component in the 4D vector by scalar
	//Stores result in current vector object
	//Throws a std::invlaid_argument if scalar is zero
	//Returns a reference to the current vector object
	//That has the result of the current Vector4 object / scalar
	Vector4& Vector4::operator/=(const double& scalar)
	{
		//if scalar is 0.0 throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		this->m_x /= scalar;
		this->m_y /= scalar;
		this->m_z /= scalar;
		this->m_w /= scalar;

		return *this;
	}

	//Stores the x and y values of Vector2s v into the x and y values of this Vector4 and sets z and w to 0.0f
	void Vector4::operator=(const Vector2& v)
	{
		this->m_x = v.x();
		this->m_y = v.y();
		this->m_z = 0.0;
		this->m_w = 0.0;
	}

	//Stores the x, y and z values of Vector3s v into the x, y and z values of this Vector4 and sets w to 0.0f
	void Vector4::operator=(const Vector3& v)
	{
		this->m_x = v.x();
		this->m_y = v.y();
		this->m_z = v.z();
		this->m_w = 0.0;
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Friend Functions for Vector4 class

	//Adds two 4D vectors
	//Vector addition for 4D is (ax + bx, ay + by, az + bz, aw + bw)
	//Returns a new Vector4 object that has the result of a + b
	Vector4 operator+(const Vector4& a, const Vector4& b)
	{
		return Vector4(a.m_x + b.m_x, a.m_y + b.m_y, a.m_z + b.m_z, a.m_w + b.m_w);
	}

	//Subtracts two 4D vectors
	//Vector addition for 4D is (ax - bx, ay - by, az - bz, aw - bw)
	//Returns a new Vector4 object that has the result of a - b
	Vector4 operator-(const Vector4& a, const Vector4& b)
	{
		return Vector4(a.m_x - b.m_x, a.m_y - b.m_y, a.m_z - b.m_z, a.m_w - b.m_w);
	}

	//Multiplies each component in the 4D vector by scalar
	//Vector multiplicaton by a scalar for 4D is (vx * scalar , vy * scalar, vz * scalar, vw * scalar)
	//Returns a new Vector4 object that has the result of v * scalar
	Vector4 operator*(const Vector4& v, const double& scalar)
	{
		return Vector4(v.m_x * scalar, v.m_y * scalar, v.m_z * scalar, v.m_w * scalar);
	}

	//Multiplies each component in the 4D vector by scalar
	//Vector multiplicaton by a scalar for 4D is (vx * scalar , vy * scalar, vz * scalar, vw * scalar)
	//Returns a new Vector4 object that has the result of v * scalar
	Vector4 operator*(const double& scalar, const Vector4& v)
	{
		return Vector4(v.m_x * scalar, v.m_y * scalar, v.m_z * scalar, v.m_w * scalar);
	}

	//Divides each componet in the 4D vector by a scalar
	//Vector division by a scalar for 4D is (vx / scalar, vy / scalar, vz / scalar, vw / scalar)
	//Throws a std::invalid_argument if scalar is zero
	//Returns a new Vector4 object with the result of v / scalar
	Vector4 operator/(const Vector4& v, const double& scalar)
	{
		//if scalar is 0.0 throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		return Vector4(v.m_x / scalar, v.m_y / scalar, v.m_z / scalar, v.m_w / scalar);
	}

	//Returns true if a is equal to b, false otherwise.
	bool operator==(const Vector4& a, const Vector4& b)
	{
		return FAUtil::compareDoubles(a.x(), b.x(), EPSILON) && FAUtil::compareDoubles(a.y(), b.y(), EPSILON) && FAUtil::compareDoubles(a.z(), b.z(), EPSILON)
			&& FAUtil::compareDoubles(a.w(), b.w(), EPSILON);
	}

	//Returns false if a is equal to b, true otherwise.
	bool operator!=(const Vector4& a, const Vector4& b)
	{
		return !(a == b);
	}

	//Returns the magnitude/length of a 4D vector
	//Length of a 4D vector is sqrt(vx * vx + vy * vy + vz * vz + vw * vw)
	double length(const Vector4& v)
	{
		return sqrt(v.m_x * v.m_x + v.m_y * v.m_y + v.m_z * v.m_z + v.m_w * v.m_w);
	}

	//Normalizes a vector
	//Returns a new Vector4 object that is a 4D unit vector
	//A unit vector is a vector that has a length of 1
	Vector4 normalize(const Vector4& v)
	{
		if (v.isZeroVector())
		{
			return v;
		}

		return v / length(v);
	}

	//Returns the distance between two 4D points
	//To get the distance between two 4D points:
	//Subtract the two 4D points and take the length of the resultant 4D point
	//length(a-b)
	//Or sqrt((ax - bx)^2 + (ay - by)^2 + (az - bz)^2 + (aw - bw)^2)
	double distance(const Vector4& a, const Vector4& b)
	{
		Vector4 c{ a - b };

		return length(c);
	}

	//Does the dot product between two 4D vectors and returns the value
	//Dot Product between two 4D vectors: ax * bx + ay * by + az * bz + aw * bw
	double dotProduct(const Vector4& a, const Vector4& b)
	{
		return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z + a.m_w * b.m_w;
	}

	//Returns the angle between two vectors
	//Angle between two vectors = acos(a dot b)
	//a and b should be unit vectors before calling this function
	double angle(const Vector4& a, const Vector4& b)
	{
		double dP = dotProduct(a, b);
		if (dP < -1.0)
		{
			dP = -1.0;
		}
		else if (dP > 1.0)
		{
			dP = 1.0;
		}

		//radians to degrees = 180 degrees / PI rad
		return acos(dP) * 180.0f / PI;
	}


	void print(Vector4 v)
	{
		std::cout << std::setprecision(15) << "(" << v.m_x << ", " << v.m_y << ", " << v.m_z << ", " << v.m_w << ")" << std::endl;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------










//Matrix4x4
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

	//Constructors for Matrix4x4 class
	
	//Default Constructor
	//Constructos an identity matrix
	//the diagonal elements are at indcies 0, 5, 10, 15
	Matrix4x4::Matrix4x4() 
		: m_matrix{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 }
	{}

	//Overloaded Constructor
	//Constructors a 4x4 matrix from the given std::array
	Matrix4x4::Matrix4x4(const double* values)
		: m_matrix{ values[0], values[4], values[8], values[12],
					values[1], values[5], values[9], values[13],
					values[2], values[6], values[10], values[14],
					values[3], values[7], values[11], values[15] }
	{}

	//Overloaded Constructor
	//Constructs a  4x4 matrix from the specified 16 elements
	//The elements are specified in row-major order
	//Elements are stored in column major order)
	Matrix4x4::Matrix4x4(double m11, double m12, double m13, double m14,
		double m21, double m22, double m23, double m24,
		double m31, double m32, double m33, double m34,
		double m41, double m42, double m43, double m44)
		: m_matrix{m11, m21, m31, m41, m12, m22, m32, m42, m13, m23, m33, m43, m14, m24, m34, m44}
	{}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Getter for Matrix4x4 class

	//Returns the element of column index as a Vector4
	//Column index starts at 0
	Vector4 Matrix4x4::column(const unsigned int& index) const
	{
		if (index > 3)
		{
			throw std::out_of_range("Not a valid column index\n");
		}

		//index * 4 to go the start of column index
		unsigned int i = index * 4;

		//go to each successive element at column index and store in a Vector4 object
		return Vector4(m_matrix[i], m_matrix[i + 1], m_matrix[i + 2], m_matrix[i + 3]);
	}

	//Returns the element of row index as a Vector4
	//Column index starts at 0
	Vector4 Matrix4x4::row(const unsigned int& index) const
	{	
		if (index > 3)
		{
			throw std::out_of_range("Not a valid row index\n");
		}

		//go to each successive element at row index and store in a Vector4 object
		return Vector4(m_matrix[index], m_matrix[index + 4], m_matrix[index + 8], m_matrix[index + 12]);
	}


	//Returns a pointer to a double that has the raw data of m_matrix
	double* Matrix4x4::data()
	{
		return m_matrix;
	}

	//Returns a pointer to a const double that has the raw data of m_matrix
	const double* Matrix4x4::data() const
	{
		return m_matrix;
	}

	//Returns a pointer to a const double that has the raw data of m_matrix
	const double* Matrix4x4::constData() const
	{
		return m_matrix;
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Setter for Matrix4x4 class

	//Sets the element at [row][col] to the specified value
	void Matrix4x4::set(const unsigned int& row, const unsigned int& col, const double& value)
	{
		if (row > 3)
		{
			throw std::out_of_range("Not a valid row index\n");
		}

		if (col > 3)
		{
			throw std::out_of_range("Not a valid column index\n");
		}

		//turn [row][col] to column major index
		//(col * # of rows per col) + row
		m_matrix[col * 4 + row] = value;
	}

	//Sets the matrix to the identity matrix
	//Sets the diagonal elements to 1.0f and the nondiagonal elements to 0.0f
	void Matrix4x4::setToIdentity()
	{
		m_matrix[0] = 1.0;
		m_matrix[1] = 0.0;
		m_matrix[2] = 0.0;
		m_matrix[3] = 0.0;

		m_matrix[4] = 0.0;
		m_matrix[5] = 1.0;
		m_matrix[6] = 0.0;
		m_matrix[7] = 0.0;

		m_matrix[8] = 0.0;
		m_matrix[9] = 0.0;
		m_matrix[10] = 1.0;
		m_matrix[11] = 0.0;

		m_matrix[12] = 0.0;
		m_matrix[13] = 0.0;
		m_matrix[14] = 0.0;
		m_matrix[15] = 1.0;
	}

	//Sets all elements in the matrix to value
	//Time Complexity : O(1)
	//Space Complexity : O(1)
	void Matrix4x4::fill(const double& value)
	{
		for (double& i : m_matrix)
		{
			i = value;
		}
	}

	//Sets the elements of the column index to the components of the Vector4 object value
	//Time Complexity : O(1)
	//Space Complexity : O(1)
	void Matrix4x4::setColumn(const unsigned int& index, const Vector4& value)
	{
		if (index > 3)
		{
			throw std::out_of_range("Not a valid row index\n");
		}

		//calculate the index that has the first value of column index
		const unsigned int i = index * 4;

		m_matrix[i] = value.x();
		m_matrix[i + 1] = value.y();
		m_matrix[i + 2] = value.z();
		m_matrix[i + 3] = value.w();
	}

	//Sets the elements of the row index to the components of the Vector4 object value
	//Time Complexity : O(1)
	//Space Complexity : O(1)
	void Matrix4x4::setRow(const unsigned int& index, const Vector4& value)
	{
		if (index > 3)
		{
			throw std::out_of_range("Not a valid row index\n");
		}

		m_matrix[index] = value.x();
		m_matrix[index + 4] = value.y();
		m_matrix[index + 8] = value.z();
		m_matrix[index + 12] = value.w();
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Member Functions

	//Returns true if the Matrix4x4 object is the identity matrix, false otherwise
	bool Matrix4x4::isIdentity() const
	{
		return	FAUtil::compareDoubles(m_matrix[0], 1.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[1], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[2], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[3], 0.0, EPSILON) &&

				FAUtil::compareDoubles(m_matrix[4], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[5], 1.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[6], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[7], 0.0, EPSILON) &&

				FAUtil::compareDoubles(m_matrix[8], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[9], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[10], 1.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[11], 0.0, EPSILON) &&

				FAUtil::compareDoubles(m_matrix[12], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[13], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[14], 0.0, EPSILON) &&
				FAUtil::compareDoubles(m_matrix[15], 1.0, EPSILON);
	}


	//Returns this matrix transposed about its diagonal
	Matrix4x4 Matrix4x4::transposed() const
	{
		Matrix4x4 t;

		//turn rows into columns
		for (unsigned int i = 0; i < 4; ++i)
		{
			t.m_matrix[i * 4] = this->m_matrix[i];
			t.m_matrix[i * 4 + 1] = this->m_matrix[i + 4];
			t.m_matrix[i * 4 + 2] = this->m_matrix[i + 8];
			t.m_matrix[i* 4 + 3] = this->m_matrix[i + 12];
		}

		return t;
	}

	//Returns a constant reference to the element at [row][col] in this matrix
	const double& Matrix4x4::operator()(const unsigned int& row, const unsigned int& col) const
	{
		if (row > 3)
		{
			throw std::out_of_range("Not a valid row index\n");
		}

		if (col > 3)
		{
			throw std::out_of_range("Not a valid column index\n");
		}

		//turn [row][col] to column major index
		//(col * # of rows per col) + row
		return m_matrix[col * 4 + row];
	}

	//Returns a reference to the element at [row][col] in this matrix
	double& Matrix4x4::operator()(const unsigned int& row, const unsigned int& col)
	{
		if (row > 3)
		{
			throw std::out_of_range("Not a valid row index\n");
		}

		if (col > 3)
		{
			throw std::out_of_range("Not a valid column index\n");
		}

		//turn [row][col] to column major index
		//(col * # of rows per col) + row
		return m_matrix[col * 4 + row];
	}

	//Adds the corresponding elements in the current Matrix4x4 object and m
	//and places the results in the current Matrix4x4 object
	//Returns a reference to the current Matrix4x4 object with the result of the
	//current Matrix4x4 object + m
	Matrix4x4& Matrix4x4::operator+=(const Matrix4x4& m)
	{
		for(unsigned int i = 0; i < 16; ++i)
		{
			this->m_matrix[i] += m.m_matrix[i];
		}

		return *this;
	}

	//Subtracts the corresponding elements in the current Matrix4x4 object and m
	//and places the results in the current Matrix4x4 object
	//Returns a reference to the current Matrix4x4 object with the result of the
	//current Matrix4x4 object - m
	Matrix4x4& Matrix4x4::operator-=(const Matrix4x4& m)
	{
		for (unsigned int i = 0; i < 16; ++i)
		{
			this->m_matrix[i] -= m.m_matrix[i];
		}

		return *this;
	}

	//Multiplies each element in the current Matrix4x4 object by scalar
	//Returns a reference to the current Matrix4x4 object with the result of the current Matrix4x4 object * scalar
	Matrix4x4& Matrix4x4::operator*=(const double& scalar)
	{
		//multiply each element in the matrix by scalar
		for (double& i : m_matrix)
		{
			i *= scalar;
		}

		return *this;
	}

	//Multiplies two matrices
	//Returns a reference to the current Matrix4x4 object with the result of the current Matrix4x4 object *  m
	Matrix4x4& Matrix4x4::operator*=(const Matrix4x4& m)
	{
		//store results in t
		Matrix4x4 t;

		for (unsigned int i = 0; i < 4; ++i)
		{
			//i row dot product with 1st column
			//store result in i row 1st column of t
			t.m_matrix[i] =	(this->m_matrix[i] * m.m_matrix[0]) +
							(this->m_matrix[i + 4] * m.m_matrix[1]) +
							(this->m_matrix[i + 8] * m.m_matrix[2]) +
							(this->m_matrix[i + 12] * m.m_matrix[3]);

			//i row dot product with 2nd column
			//store result in i row 2nd column of t
			t.m_matrix[i + 4] = (this->m_matrix[i] * m.m_matrix[4]) +
								(this->m_matrix[i + 4] * m.m_matrix[5]) +
								(this->m_matrix[i + 8] * m.m_matrix[6]) +
								(this->m_matrix[i + 12] * m.m_matrix[7]);

			//i row dot product with 3rd column
			//store result in i row 3rd column of t
			t.m_matrix[i + 8] =	(this->m_matrix[i] * m.m_matrix[8]) +
								(this->m_matrix[i + 4] * m.m_matrix[9]) +
								(this->m_matrix[i + 8] * m.m_matrix[10]) +
								(this->m_matrix[i + 12] * m.m_matrix[11]);

			//i row dot product with 4th column
			//store result in i row 4th column of t
			t.m_matrix[i + 12] = (this->m_matrix[i] * m.m_matrix[12]) +
								 (this->m_matrix[i + 4] * m.m_matrix[13]) +
								 (this->m_matrix[i + 8] * m.m_matrix[14]) +
								 (this->m_matrix[i + 12] * m.m_matrix[15]);

		}

		//copy t's matrix into this matrix
		memcpy(this->m_matrix, t.m_matrix, sizeof(Matrix4x4));

		return *this;
	}

	//Divides each element in matrix the current Matrix4x4 object by scalar
	//Returns a reference to the current Matrix4x4 object with the result of the current Matrix4x4 object / scalar
	Matrix4x4& Matrix4x4::operator/=(const double& scalar)
	{
		//if scalar is 0.0f throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}
		
		//divde each element in the matrix by scalar
		for (double& i : m_matrix)
		{
			i /= scalar;
		}
		return *this;
	}

	//Constructs a 4x4 matrix that rotates angle degrees about vector v and is multipled by this matrix
	//c = cos(angle)
	//s = sin(angle)
	//vx^2 * (1 - c) + c						vx * vy * (1 - c) - vz * s			vx * vz * (1 - c) + vy * s		0
	//vx * vy * (1 - c) + vz * s				vy^2 * (1 - c) + c					vy * vz * (1 - c) - vx * s		0
	//vx * vz * (1 - c) - vy * s				vy * vz * (1 - c) + vx * s			vz^2 * (1 - c) + c				0
	//		0											0									0						1
	void Matrix4x4::rotate(double angle, const Vector3& v)
	{
		Matrix4x4 rot;
		//convert from degrees to radians
		double rad = angle * PI / 180.0;

		//cos(angle)
		double c = cos(rad);

		//sin(angle)
		double s = sin(rad);

		rot.m_matrix[0] = v.x() * v.x() * (1 - c) + c;
		rot.m_matrix[1] = v.x() * v.y() * (1 - c) + v.z() * s;
		rot.m_matrix[2] = v.x() * v.z() * (1 - c) - v.y() * s;

		rot.m_matrix[4] = v.x() * v.y() * (1 - c) - v.z() * s;
		rot.m_matrix[5] = v.y() * v.y() * (1 - c) + c;
		rot.m_matrix[6] = v.y() * v.z() * (1 - c) + v.x() * s;

		rot.m_matrix[8] = v.x() * v.z() * (1 - c) + v.y() * s;
		rot.m_matrix[9] = v.y() * v.z() * (1 - c) - v.x() * s;
		rot.m_matrix[10] = v.z() * v.z() * (1 - c) + c;

		//multiply this matrix by the rot matirx and store result in this matrix
		*this *= rot;
	}

	//Constructs a 4x4 matrix that rotates angle degrees about vector (x, y, z) and is multipled by this matrix
	//c = cos(angle)
	//s = sin(angle)
	//x^2 * (1 - c) + c						x * y * (1 - c) - z * s			x * z * (1 - c) + y * s		0
	//x * y * (1 - c) + z * s				y^2 * (1 - c) + c				y * z * (1 - c) - x * s		0
	//x * z * (1 - c) - y * s				y * z * (1 - c) + x * s			z^2 * (1 - c) + c			0
	//		0										0								0					1
	void Matrix4x4::rotate(const double& angle, const double& x, const double& y, const double& z)
	{
		Matrix4x4 rot;

		//convert from degrees to radians
		double rad = angle * PI / 180.0;

		//cos(angle)
		double c = cos(rad);

		//sin(angle)
		double s = sin(rad);

		rot.m_matrix[0] = x * x * (1 - c)  + c;
		rot.m_matrix[1] = x * y * (1 - c) + z * s;
		rot.m_matrix[2] = x * z * (1 - c) - y * s;

		rot.m_matrix[4] = x * y * (1 - c) - z * s;
		rot.m_matrix[5] = y * y * (1 - c) + c;
		rot.m_matrix[6] = y * z * (1 - c) + x * s;

		rot.m_matrix[8] = x * z * (1 - c) + y * s;
		rot.m_matrix[9] = y * z * (1 - c) - x * s;
		rot.m_matrix[10] = z * z * (1 - c) + c;


		//multiply this matrix by the rot matirx and store result in this matrix
		*this *= rot;
	}

	//Multiplies this matrix by another that rotates the coordinates using the quaternion matrix.
	//The x, y and z values in the Vector4 is the axis you want to rotate around normalized
	//The w value is the angle in degrees
	//x = xsin(w / 2)
	//y = ysin(w / 2)
	//z = zsin(w / 2)
	//w = cos(w / 2)
	//1 - 2y^2 - 2z^2	2xy - 2wz			2xz + 2wy			0
	//2xy + 2wz			1 - 2x^2 - 2z^2		2yz - 2wx			0
	//2xz - 2wy			2yz + 2wx			1- 2x^2 - 2y^2		0
	//0						0					0				1
	void Matrix4x4::rotateUsingQuaternion(const Vector4& v)
	{
		Matrix4x4 r;
		
		double x = v.x() * sin((v.w() / 2) * PI / 180.0);
		double y = v.y() * sin((v.w() / 2) * PI / 180.0);
		double z = v.z() * sin((v.w() / 2) * PI / 180.0);
		double w = cos((v.w() / 2) * PI / 180.0);

		r.m_matrix[0] = 1 - 2 * y * y - 2 * z * z;
		r.m_matrix[1] = 2 * x * y + 2 * w * z;
		r.m_matrix[2] = 2 * x * z - 2 * w * y;

		r.m_matrix[4] = 2 * x * y - 2 * w * z;
		r.m_matrix[5] = 1 - 2 * x * x - 2 * z * z;
		r.m_matrix[6] = 2 * y * z + 2 * w * x;

		r.m_matrix[8] = 2 * x * z + 2 * w * y;
		r.m_matrix[9] = 2 * y * z - 2 * w * x;
		r.m_matrix[10] = 1 - 2 * x * x - 2 * y * y;

		//multiply this matrix by the rot matirx and store result in this matrix
		*this *= r;
	}

	//Multiplies this matrix by another that rotates the coordinates using the quaternion matrix
	//v is the axis to rotate around normalized
	//The angle should be given in degrees
	//x = xsin(angle / 2)
	//y = ysin(angle / 2)
	//z = zsin(angle / 2)
	//w = cos(angle / 2)
	//1 - 2y^2 - 2z^2	2xy - 2wz			2xz + 2wy			0
	//2xy + 2wz			1 - 2x^2 - 2z^2		2yz - 2wx			0
	//2xz - 2wy			2yz + 2wx			1- 2x^2 - 2y^2		0
	//0						0					0				1
	void Matrix4x4::rotateUsingQuaternion(const double& angle, const Vector3& v)
	{
		Matrix4x4 r;

		double x = v.x() * sin((angle / 2) * PI / 180.0);
		double y = v.y() * sin((angle / 2) * PI / 180.0);
		double z = v.z() * sin((angle / 2) * PI / 180.0);
		double w = cos((angle / 2) * PI / 180.0);

		r.m_matrix[0] = 1 - 2 * y * y - 2 * z * z;
		r.m_matrix[1] = 2 * x * y + 2 * w * z;
		r.m_matrix[2] = 2 * x * z - 2 * w * y;

		r.m_matrix[4] = 2 * x * y - 2 * w * z;
		r.m_matrix[5] = 1 - 2 * x * x - 2 * z * z;
		r.m_matrix[6] = 2 * y * z + 2 * w * x;

		r.m_matrix[8] = 2 * x * z + 2 * w * y;
		r.m_matrix[9] = 2 * y * z - 2 * w * x;
		r.m_matrix[10] = 1 - 2 * x * x - 2 * y * y;

		//multiply this matrix by the rot matirx and store result in this matrix
		*this *= r;
	}

	//Multiplies this matrix by another that rotates the coordinates using the quaternion matrix
	//(x, y, z) is the axis to rotate around normalized.
	//angle is in degrees
	//x = xsin(angle / 2)
	//y = ysin(angle / 2)
	//z = zsin(angle / 2)
	//w = cos(angle / 2)
	//1 - 2y^2 - 2z^2	2xy - 2wz			2xz + 2wy			0
	//2xy + 2wz			1 - 2x^2 - 2z^2		2yz - 2wx			0
	//2xz - 2wy			2yz + 2wx			1- 2x^2 - 2y^2		0
	//0						0					0				1
	void Matrix4x4::rotateUsingQuaternion(const double& angle, const double& x, const double& y, const double& z)
	{
		Matrix4x4 r;

		double newX = x * sin((angle / 2) * PI / 180.0);
		double newY = y * sin((angle / 2) * PI / 180.0);
		double newZ = z * sin((angle / 2) * PI / 180.0);
		double w = cos((angle / 2) * PI / 180.0);

		r.m_matrix[0] = 1 - 2 * newY * newY - 2 * newZ * newZ;
		r.m_matrix[1] = 2 * newX * newY + 2 * w * newZ;
		r.m_matrix[2] = 2 * newX * newZ - 2 * w * newY;

		r.m_matrix[4] = 2 * newX * newY - 2 * w * newZ;
		r.m_matrix[5] = 1 - 2 * newX * newX - 2 * newZ * newZ;
		r.m_matrix[6] = 2 * newY * newZ + 2 * w * newX;

		r.m_matrix[8] = 2 * newX * newZ + 2 * w * newY;
		r.m_matrix[9] = 2 * newY * newZ - 2 * w * newX;
		r.m_matrix[10] = 1 - 2 * newX * newX - 2 * newY * newY;

		//multiply this matrix by the rot matirx and store result in this matrix
		*this *= r;
	}

	//Constructs a matrix thats scales the coordinates by the components of vector v and is multiplied by this matrix)
	//vx	 0		0		0
	//0		 vy     0		0
	//0		 0		vz		0
	//0		 0		0		1
	void Matrix4x4::scale(const Vector3& v)
	{
		Matrix4x4 s;

		//set the diagonals of the matrix to the respective components of vector v
		s.m_matrix[0] = v.x();
		s.m_matrix[5] = v.y();
		s.m_matrix[10] = v.z();

		//multiply this matrix by the s matirx and store result in this matrix
		*this *= s;
	}

	//Constructs a matrix thats scales the coordinates by the components of x and y and is multiplied by this matrix
	//x		 0		0		0
	//0		 y		0		0
	//0		 0		1		0
	//0		 0		0		1
	void Matrix4x4::scale(const double& x, const double& y)
	{
		Matrix4x4 s;

		//set the diagonals of the matrix to the respective components of x and y
		s.m_matrix[0] = x;
		s.m_matrix[5] = y;

		//multiply this matrix by the s matirx and store result in this matrix
		*this *= s;
	}

	//Constructs a matrix thats scales the coordinates by the components of x, y and z and is multiplied by this matrix
	//x		 0		0		0
	//0		 y		0		0
	//0		 0		z		0
	//0		 0		0		1
	void Matrix4x4::scale(const double& x, const double& y, const double& z)
	{
		Matrix4x4 s;

		//set the diagonals of the matrix to the respective components of x, y and z
		s.m_matrix[0] = x;
		s.m_matrix[5] = y;
		s.m_matrix[10] = z;

		//multiply this matrix by the s matirx and store result in this matrix
		*this *= s;
	}

	//Constructs a matrix thats scales the coordinates by the specified factor and is multiplied by this matrix
	//factor		  0				   0			0
	//0				factor			   0			0
	//0				  0				factor			0
	//0				  0				   0			1
	void Matrix4x4::scale(const double& factor)
	{
		Matrix4x4 s;

		//set the diagonals of the matrix to factor
		s.m_matrix[0] = factor;
		s.m_matrix[5] = factor;
		s.m_matrix[10] = factor;

		//multiply this matrix by the s matirx and store result in this matrix
		*this *= s;
	}

	//Constructs a matrix thats scales the coordinates by the specified factor along vector v and is multiplied by this matrix
	//1 + (factor - 1) * vx * vx		(factor - 1) * vx * vy			(factor - 1) * vx * vz			0
	//(factor - 1) * v.x * v.y			1 + (factor - 1) * vy * vy		factor - 1) * vy * vz			0
	//(factor - 1) * vx * vz			(factor - 1) * vy * vz			1 + (factor - 1) * vz * vz		0
	//		0										0								0					1
	void Matrix4x4::scale(const Vector3& v, const double& factor)
	{
		Matrix4x4 s;

		s.m_matrix[0] = 1 + (factor - 1) * v.x() * v.x();
		s.m_matrix[1] = (factor - 1) * v.x() * v.y();
		s.m_matrix[2] = (factor - 1) * v.x() * v.z();

		s.m_matrix[4] = (factor - 1) * v.x() * v.y();
		s.m_matrix[5] = 1 + (factor - 1) * v.y() * v.y();
		s.m_matrix[6] = (factor - 1) * v.y() * v.z();

		s.m_matrix[8] = (factor - 1) * v.x() * v.z();
		s.m_matrix[9] = (factor - 1) * v.y() * v.z();
		s.m_matrix[10] = 1 + (factor - 1) * v.z() * v.z();

		//multiply this matrix by the s matirx and store result in this matrix
		*this *= s;
	}

	//Multiplies this matrix by another that translates coordinates by the components of v
	//1			0				0			vx
	//0			1				0			vy
	//0			0				1			vz
	//0			0				0			1
	void Matrix4x4::translate(const Vector3& v)
	{
		Matrix4x4 t;

		//set the last column to (vx, vy, vz, 1.0)
		t.m_matrix[12] = v.x();
		t.m_matrix[13] = v.y();
		t.m_matrix[14] = v.z();

		//multiply this matrix by the t matirx and store result in this matrix
		*this *= t;
	}

	//Multiplies this matrix by another that translates coordinates by the components of x and y
	//1		 0				0			x
	//0		 1				0			y
	//0		 0				1			0
	//0		 0				0			1
	void Matrix4x4::translate(const double& x, const double& y)
	{
		Matrix4x4 t;

		//set the last column to (vx, vy, vz, 1.0)
		t.m_matrix[12] = x;
		t.m_matrix[13] = y;


		//multiply this matrix by the t matirx and store result in this matrix
		*this *= t;
	}

	//Multiplies this matrix by another that translates coordinates by the components of x, y and z
	//1			0				0			x
	//0			1				0			y
	//0			0				1			z
	//0			0				0			1
	void Matrix4x4::translate(const double& x, const double& y, const double& z)
	{
		Matrix4x4 t;

		//set the last column to (vx, vy, vz, 1.0)
		t.m_matrix[12] = x;
		t.m_matrix[13] = y;
		t.m_matrix[14] = z;

		//multiply this matrix by the t matirx and store result in this matrix
		*this *= t;
	}

	//Constructs a matrix that does othrographic projection for a window with lower left corner at (left, bottom), 
	//upper right corner at (right, top) and the specified near and far clipping planes
	//Then multiplies this matrix by the constructed matrix
	//2 / (right - left)						0							0						-(right + left) / (right - left)
	//		0							2 / (top - bottom)					0						-(top + bottom) / (top - bottom)
	//		0									0					  -2 / (far - near)				-(far + near)/(far - near)	
	//		0									0							0									1
	void Matrix4x4::ortho(const double& left, const double& right, const double& bottom, const double& top, const double& near, const double& far)
	{
		Matrix4x4 t;

		//set the diagonals to the scaling factor
		t.m_matrix[0] = 2 / (right - left);
		t.m_matrix[5] = 2 / (top - bottom);
		t.m_matrix[10] = -2 / (far - near);

		//set the last column to the translation values
		t.m_matrix[12] = -(right + left) / (right - left);
		t.m_matrix[13] = -(top + bottom) / (top - bottom);
		t.m_matrix[14] = -(far + near) / (far - near);
		
		//multiply this matrix by t and store result in this matrix
		*this *= t;
	}

	//Multiplies this matrix by another that applies a persepective projection.
	//The fov is the vertical angle in degrees
	//Aspect ratio is the aspect ratio of your window
	//Near and far are the distances from the viewer to the corresponding planes
	//(h/w)(1/tan(fov / 2))							0							0							0 
	//		0								(1/tan(fov / 2))					0							0
	//		0										0					  far/(far - near)		 -(far * near)/(far - near)	
	//		0										0							1							0
	void Matrix4x4::perspective(const double& fov, const double& aspectRatio, const double& near, const double& far)
	{
		Matrix4x4 t;

		double f = fov / 2;
		t.m_matrix[0] = (1 / aspectRatio) * (1 / tan(f * PI / 180.0f));
		t.m_matrix[5] = 1 / tan(f * PI / 180.0f);
		t.m_matrix[10] = far / (far - near);
		t.m_matrix[11] = 1.0f;
		t.m_matrix[14] = -(far * near) / (far - near);
		t.m_matrix[15] = 0.0f;

		*this *= t;
	}

	//Returns the determinant of this matrix
	//detM = m_matrix[i][j] * (-1)^(i + j) * M(i, j)
	//i = row, j = col, M = minor
	//Go through each element in the first row and multiply it by its cofactor((-1)^(0 + j) * M(0, j))
	double Matrix4x4::determinant() const
	{
		//copy all elements in m_matrix into m
		//want the values to be double for larger range of values
		double m[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			m[i] = m_matrix[i];
		}

		//m11[m22(m33m44 - m34m43) + (m23(m34m42 - m32m44) + m24(m32m43 - m33m42)]
		//{0}[{5}({10}{15} - {14}{11}) + {9}({14}{7} - {6}{15}) + {13}({6}{11} - {10}{7})}
		double cofactor1{ m[5] * ((m[10] * m[15]) - (m[14] * m[11])) };
		double cofactor2{ m[9] * ((m[14] * m[7]) - (m[6] * m[15])) };
		double cofactor3{ m[13] * ((m[6] * m[11]) - (m[10] * m[7])) };

		//-m12[m21(m33m44 - m34m43) + m23(m34m41 - m31m44) + m24(m31m43 - m33m41)]
		//{4}[{1}({10}{15} - {14}{11}) + {9}({14}{3} - {2}{15}) + {13}({2}{11} - {10}{3})]
		double cofactor4{ m[1] * ((m[10] * m[15]) - (m[14] * m[11])) };
		double cofactor5{ m[9] * ((m[14] * m[3]) - (m[2] * m[15])) };
		double cofactor6{ m[13] * ((m[2] * m[11]) - (m[10] * m[3])) };

		//m13[m21(m32m44 - m34m42) + m22(m34m41 - m31m44) + m24(m31m42 - m32m41)]
		//{8}[{1}({6}{15} - {14}{7}) + {5}({14}{3} - {2}{15}) + {13}({2}{7} - {6}{3})]
		double cofactor7{ m[1] * ((m[6] * m[15]) - (m[14] * m[7])) };
		double cofactor8{ m[5] * ((m[14] * m[3]) - (m[2] * m[15])) };
		double cofactor9{ m[13] * ((m[2] * m[7]) - (m[6] * m[3])) };

		//-m14[m21(m32m43 - m33m42) + m22(m33m41 - m31m43) + m23(m31m42 - m32m41)]
		//{12}[{1}({6}{11} - {10}{7}) + {5}({10}{3} - {2]{11}) + {9}({2}{7} - {6}{3})]
		double cofactor10{ m[1] * ((m[6] * m[11]) - (m[10] * m[7])) };
		double cofactor11{ m[5] * ((m[10] * m[3]) - (m[2] * m[11])) };
		double cofactor12{ m[9] * ((m[2] * m[7]) - (m[6] * m[3])) };


		return (m[0] * (cofactor1 + cofactor2 + cofactor3)) - (m[4] * (cofactor4 + cofactor5 + cofactor6)) +
			   (m[8] * (cofactor7 + cofactor8 + cofactor9)) - (m[12] * (cofactor10 + cofactor11 + cofactor12));
	}

	//Minor of a given row and column of a matrix
	//Get rid of the row and col numbers and take the determinant of remaining 3x3 matrix
	double Matrix4x4::minor(const unsigned int& row, const unsigned int& col) const
	{
		double mat[9]{};
		int index{ 0 };
		//copy all of the numbers not in the given row or col into mat
		for (int i = 0; i < 16; ++i)
		{
			if (i == row || i == row + 4 || i == row + 8 || i == row + 12 || i == col * 4 || i == col * 4 + 1 || i == col * 4 + 2 || i == col * 4 + 3)
				continue;

			mat[index] = m_matrix[i];
			++index;
		}

		//calculate the the determiant of a 3x3 matrix
		//m11m22m33 + m12m23m31 + m12m21m32 - m13m22m31 - m12m21m33 - m11m23m32
		double m = (mat[0] * mat[4] * mat[8]) + (mat[3] * mat[7] * mat[2]) + (mat[6] * mat[1] * mat[5])
			- (mat[6] * mat[4] * mat[2]) - (mat[3] * mat[1] * mat[8]) - (mat[0] * mat[7] * mat[5]);

		return m;
	}

	//Returns the adjoint of this matrix
	//adjM is a 4x4 matrix of the cofactors of this matrix
	Matrix4x4 Matrix4x4::adjoint() const
	{
		Matrix4x4 adj;
		int index{ 0 };
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j, ++index)
			{
				double cofactor = pow(-1, i + j)  * minor(i, j);
				adj.m_matrix[index] = cofactor;
			}
		}

		return adj;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Friend non-member functions

	//Adds two 4x4 matrices
	//Returns a Matrix4x4 object with the result of m1 + m2
	Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		Matrix4x4 t;

		for (unsigned int i = 0; i < 16; ++i)
		{
			t.m_matrix[i] = m1.m_matrix[i] + m2.m_matrix[i];
		}

		return t;
	}

	//Subtracts two 4x4 matrices
	//Returns a Matrix4x4 object with the result of m1 - m2
	Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		Matrix4x4 t;

		for (unsigned int i = 0; i < 16; ++i)
		{
			t.m_matrix[i] = m1.m_matrix[i] - m2.m_matrix[i];
		}

		return t;
	}

	//Negates the 4x4 matrix
	//Returns a Matrix4x4 object with the result of -m
	Matrix4x4 operator-(Matrix4x4& m)
	{
		Matrix4x4 t;
		//negate each element in the matrix
		for(unsigned int i = 0; i < 16; ++i)
		{
			t.m_matrix[i] = -m.m_matrix[i];
		}

		return t;
	}

	//Multiplies each element in matrix m1 by scalar
	//Returns a Matrix4x4 object with the result of m1 * scalar
	Matrix4x4 operator*(const Matrix4x4& m1, const double& scalar)
	{
		Matrix4x4 t;

		for (unsigned int i = 0; i < 16; ++i)
		{
			t.m_matrix[i] = m1.m_matrix[i] * scalar;
		}

		return t;
	}

	//Multiplies each element in matrix m1 by scalar
	//Returns a Matrix4x4 object with the result of scalar * m1
	Matrix4x4 operator*(const double& scalar, const Matrix4x4& m1)
	{
		Matrix4x4 t;

		for (unsigned int i = 0; i < 16; ++i)
		{
			t.m_matrix[i] = scalar * m1.m_matrix[i];
		}

		return t;
	}

	//Multiplies two matrices
	//Returns a Matrix4x4 object with the result of m1 * m2
	Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		//store results in t
		Matrix4x4 t;

		for (unsigned int i = 0; i < 4; ++i)
		{
			//i row dot product with 1st column
			//store result in i row 1st column of t
			t.m_matrix[i] = (m1.m_matrix[i] * m2.m_matrix[0]) +
							(m1.m_matrix[i + 4] * m2.m_matrix[1]) +
							(m1.m_matrix[i + 8] * m2.m_matrix[2]) +
							(m1.m_matrix[i + 12] * m2.m_matrix[3]);

			//i row dot product with 2nd column
			//store result in i row 2nd column of t
			t.m_matrix[i + 4] = (m1.m_matrix[i] * m2.m_matrix[4]) +
								(m1.m_matrix[i + 4] * m2.m_matrix[5]) +
								(m1.m_matrix[i + 8] * m2.m_matrix[6]) +
								(m1.m_matrix[i + 12] * m2.m_matrix[7]);

			//i row dot product with 3rd column
			//store result in i row 3rd column of t
			t.m_matrix[i + 8] = (m1.m_matrix[i] * m2.m_matrix[8]) +
								(m1.m_matrix[i + 4] * m2.m_matrix[9]) +
								(m1.m_matrix[i + 8] * m2.m_matrix[10]) +
								(m1.m_matrix[i + 12] * m2.m_matrix[11]);

			//i row dot product with 4th column
			//store result in i row 4th column of t
			t.m_matrix[i + 12] = (m1.m_matrix[i] * m2.m_matrix[12]) +
								 (m1.m_matrix[i + 4] * m2.m_matrix[13]) +
								 (m1.m_matrix[i + 8] * m2.m_matrix[14]) +
								 (m1.m_matrix[i + 12] * m2.m_matrix[15]);
		}

		return t;
	}

	//Multiplies a 4x4 matrix with a column vector(4x1)
	//Returns a Vector4 object with the result of m * vec
	Vector4 operator*(const Matrix4x4& m, const Vector4& vec)
	{
		//Dot Product each row in m matrix with vec
		double x = m.m_matrix[0] * vec.x() + m.m_matrix[4] * vec.y() + m.m_matrix[8] * vec.z() + m.m_matrix[12] * vec.w();
		double y = m.m_matrix[1] * vec.x() + m.m_matrix[5] * vec.y() + m.m_matrix[9] * vec.z() + m.m_matrix[13] * vec.w();
		double z = m.m_matrix[2] * vec.x() + m.m_matrix[6] * vec.y() + m.m_matrix[10] * vec.z() + m.m_matrix[14] * vec.w();
		double w = m.m_matrix[3] * vec.x() + m.m_matrix[7] * vec.y() + m.m_matrix[11] * vec.z() + m.m_matrix[15] * vec.w();

		return Vector4(x, y, z, w);
	}

	//Multiplies a 4x4 matrix with a row vector(1x4)
	//Returns a Vector4 object with the result of vec * m
	Vector4 operator*(const Vector4& vec, const Matrix4x4& m)
	{
		//Dot Product each column in m matrix with vec
		double x = m.m_matrix[0] * vec.x() + m.m_matrix[1] * vec.y() + m.m_matrix[2] * vec.z() + m.m_matrix[3] * vec.w();
		double y = m.m_matrix[4] * vec.x() + m.m_matrix[5] * vec.y() + m.m_matrix[6] * vec.z() + m.m_matrix[7] * vec.w();
		double z = m.m_matrix[8] * vec.x() + m.m_matrix[9] * vec.y() + m.m_matrix[10] * vec.z() + m.m_matrix[11] * vec.w();
		double w = m.m_matrix[12] * vec.x() + m.m_matrix[13] * vec.y() + m.m_matrix[14] * vec.z() + m.m_matrix[15] * vec.w();

		return Vector4(x, y, z, w);
	}

	//Divides each element in matrix m1 by scalar
	//Returns a Matrix4x4 object with the result of m1 / scalar
	Matrix4x4 operator/(const Matrix4x4& m1, const double& scalar)
	{
		//if scalar is 0.0 throw
		if (FAUtil::compareDoubles(scalar, 0.0, EPSILON))
		{
			throw std::invalid_argument("Scalar argument is zero\n");
		}

		Matrix4x4 t;

		for (unsigned int i = 0; i < 16; ++i)
		{
			t.m_matrix[i] = m1.m_matrix[i] / scalar;
		}

		return t;
	}

	///Returns true if m1 is identical to m2, false otherwise.
	bool operator==(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		for (int i = 0; i < 16; ++i)
		{
			if (!FAUtil::compareDoubles(m1.m_matrix[i], m2.m_matrix[i], EPSILON))
				return false;
		}

		return true;
	}

	//Returns false if m1 is identical to m2, true otherwise.
	bool operator!=(const Matrix4x4& m1, const Matrix4x4& m2)
	{
		return !(m1 == m2);
	}

	//Returns the inverse of the matrix m.
	//If the matrix can't be inverted then the identity matrix is returned.
	//Inverse(M) = adjM / det of M;
	Matrix4x4 inverse(const Matrix4x4& m)
	{
		double det = m.determinant();
		Matrix4x4 inv;
		if (FAUtil::compareDoubles(det, 0.0, EPSILON))
		{
			return inv;
		}

		Matrix4x4 adj = m.adjoint();
		inv = adj / det;

		return inv;
	}

	void print(const Matrix4x4& m)
	{
		for (unsigned int i = 0; i < 4; ++i)
		{
			std::cout << std::setprecision(20);
			std::cout << m.m_matrix[i * 4] << " ";
			std::cout << m.m_matrix[i * 4 + 1] << " ";
			std::cout << m.m_matrix[i * 4 + 2] << " ";
			std::cout << m.m_matrix[i * 4 + 3] << std::endl;
		}
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------


//QUATERNION
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Constructors for Quaternion class

	//Creates a new Quaternion with scalar value = 1 and the 3D vector = (0, 0, 0).
	Quaternion::Quaternion() : m_w{ 1 }, m_v{ 0.0, 0.0, 0.0 }
	{}

	//Creates a new Quaternion with scalar value equal to w and
	//the 3D vector equal to v.
	Quaternion::Quaternion(const double& w, const Vector3& v) : m_w{ w }, m_v{ v }
	{}

	//Creates a new Quaternion with scalar value equal to w and
	//the 3D vector equal to (x, y ,z).
	Quaternion::Quaternion(const double& w, const double& x, const double& y, const double& z) : m_w{ w }, m_v{ x, y, z }
	{}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Getters for Quaternion Class
	
	//Returns the scalar component of the quaternion.
	double Quaternion::scalar() const
	{
		return m_w;
	}

	//Returns the 3D vector component of the quaternion.
	Vector3 Quaternion::vector() const
	{
		return m_v;
	}

	//Returns the x component of the quaternion's 3D vector.
	double Quaternion::x() const
	{
		return m_v.x();
	}

	//Returns the y component of the quaternion's 3D vector.
	double Quaternion::y() const
	{
		return m_v.y();
	}
	
	//Returns the z component of the quaternion's 3D vector.
	double Quaternion::z() const
	{
		return m_v.z();
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

	//Setters for Quaternion class

	//Sets the quaternion values.
	//The scalar value equals to w and the 3D vector equal to(x, y, z).
	void Quaternion::setQuaternion(const double& w, const Vector3& v)
	{
		m_w = w;
		m_v = v;
	}

	//Sets the quaternion values.
	//The scalar value equals to w and the 3D vector equals to (x, y, z).
	void Quaternion::setQuaternion(const double& w, double& x, double& y, double& z)
	{
		m_w = w;
		m_v.setX(x);
		m_v.setY(y);
		m_v.setZ(z);
	}

	//Sets the scalar value in the quaternion to w.
	void Quaternion::setScalar(const double& w)
	{
		m_w = w;
	}

	//Sets the 3D vector in the quaternion to v.
	void Quaternion::setVector(const Vector3& v)
	{
		m_v = v;
	}

	//Sets the 3D vector in the quaternion to(x, y, z).
	void Quaternion::setVector(const double& x, const double& y, const double& z)
	{
		m_v = Vector3(x, y, z);
	}

	//Sets the x component of the 3D vector in the quaternion to the given x.
	void Quaternion::setX(const double& x)
	{
		m_v.setX(x);
	}

	//Sets the y component of the 3D vector in the quaternion to the given y.
	void Quaternion::setY(const double& y)
	{
		m_v.setY(y);
	}

	//Sets the z component of the 3D vector in the quaternion to the given z.
	void Quaternion::setZ(const double& z)
	{
		m_v.setZ(z);
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Member Functions for Quternion class

	//Returns true if the quaternion scalar value is equal to 0 and all the components of the 3D vector is equal to zero, false otherwise.
	bool Quaternion::isZeroQuaternion() const
	{
		return FAUtil::compareDoubles(m_w, 0.0, EPSILON) && FAUtil::compareDoubles(m_v.x(), 0.0, EPSILON) &&
			FAUtil::compareDoubles(m_v.y(), 0.0, EPSILON) && FAUtil::compareDoubles(m_v.z(), 0.0, EPSILON);
	}

	//Multiplies this quaternion by a scalar and stores the result in this quaternion.
	Quaternion& Quaternion::operator*=(const double& k)
	{
		m_w *= k;
		m_v *= k;

		return *this;
	}

	//Multiplies this quaternion and q and stores result in this quaternion.
	//q1q2 = [w1 (x1 y1 z1)][w2 (x2 y2 z2)] = [w1 v1][w2 v2] = [w1w2 - v1 dot v2 (w1v2 + w2v1 + v1 x v2)]
	Quaternion& Quaternion::operator*=(const Quaternion& q)
	{
		double scalar = m_w * q.m_w - dotProduct(m_v, q.m_v);
		Vector3 v = m_w * q.m_v + q.m_w * m_v + crossProduct(m_v, q.m_v);

		m_w = scalar;
		m_v = v;

		return *this;
	}

	//Adds this quaternion to the given quaternion q and returns a reference to this quaternion.
	Quaternion& Quaternion::operator+=(const Quaternion& q)
	{
		m_w += q.m_w;
		m_v += q.m_v;

		return *this;
	}

	//Subtracts this quaternion from the given quaternion q and returns a reference to this quaternion.
	Quaternion& Quaternion::operator-=(const Quaternion& q)
	{
		m_w -= q.m_w;
		m_v -= q.m_v;

		return *this;
	}

	//Creates a rotation matrix from this quaternion.
	//Normalize the quaternion before using this function.
	//1 - 2y^2 - 2z^2	2xy - 2wz			2xz + 2wy			0
	//2xy + 2wz			1 - 2x^2 - 2z^2		2yz - 2wx			0
	//2xz - 2wy			2yz + 2wx			1- 2x^2 - 2y^2		0
	//0						0					0				1
	Matrix4x4 Quaternion::toRotationMatrix()
	{
		Matrix4x4 rot;

		double* r = rot.data();
		r[0] = 1 - 2 * m_v.y() * m_v.y() - 2 * m_v.z() * m_v.z();
		r[1] = 2 * m_v.x() * m_v.y() + 2 * m_w * m_v.z();
		r[2] = 2 * m_v.x() * m_v.z() - 2 * m_w * m_v.y();

		r[4] = 2 * m_v.x() * m_v.y() - 2 * m_w * m_v.z();
		r[5] = 1 - 2 * m_v.x() * m_v.x() - 2 * m_v.z() * m_v.z();
		r[6] = 2 * m_v.y() * m_v.z() + 2 * m_w * m_v.x();

		r[8] = 2 * m_v.x() * m_v.z() + 2 * m_w * m_v.y();
		r[9] = 2 * m_v.y() * m_v.z() - 2 * m_w * m_v.x();
		r[10] = 1 - 2 * m_v.x() * m_v.x() - 2 * m_v.y() * m_v.y();

		return rot;
	}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
	//Friend Functions for Quaternion class

	//Negates the scalar value of q and each componenet in the 3D vector of q.
	Quaternion operator-(const Quaternion& q)
	{
		return Quaternion(-q.m_w, -q.x(), -q.y(), -q.z());
	}

	//Returns the magnitude of the quaternion.
	//The magnitude of a quaternion is sqrt(w^2 + x^2 + y^2 + z^2)
	double length(const Quaternion& q)
	{
		return sqrt(q.m_w * q.m_w + q.x() * q.x() + q.y() * q.y() + q.z() * q.z());
	}

	//Normalizes the quaternion q.
	//returns a Quaternion object that has the result q / |q|.
	Quaternion normalize(const Quaternion& q)
	{
		if (q.isZeroQuaternion())
		{
			return q;
		}

		double mag = length(q);
		return Quaternion(q.m_w / mag, q.x() / mag, q.y() / mag, q.z() / mag);
	}

	//Returns the conjugate of quaternion q.
	Quaternion conjugate(const Quaternion& q)
	{
		return Quaternion(q.m_w, -q.x(), -q.y(), -q.z());
	}

	//Returns the inverse of quaternion q.
	Quaternion inverse(const Quaternion& q)
	{
		if (q.isZeroQuaternion())
		{
			return q;
		}

		Quaternion conjugateOfQ = conjugate(q);
		double mag = length(q);

		return Quaternion(conjugateOfQ.m_w / mag, conjugateOfQ.x() / mag, conjugateOfQ.y() / mag, conjugateOfQ.z() / mag);
	}

	//Returns the product of q1 and q2 using quaternion multiplication.
	//q1q2 = [w1 (x1 y1 z1)][w2 (x2 y2 z2)] = [w1 v1][w2 v2] = [w1w2 - v1 dot v2 (w1v2 + w2v1 + v1 x v2)]
	Quaternion operator*(const Quaternion& q1, const Quaternion& q2)
	{
		double scalar = q1.m_w * q2.m_w - dotProduct(q1.m_v, q2.m_v);
		Vector3 v = q1.m_w * q2.m_v + q2.m_w * q1.m_v + crossProduct(q1.m_v, q2.m_v);

		return Quaternion(scalar, v);
	}

	//Returns a Quaternion object that has the result of q * k.
	Quaternion operator*(const Quaternion& q, const double& k)
	{
		return Quaternion(q.m_w * k, q.m_v * k);
	}

	//Returns a Quaternion object that has the result of k * q.
	Quaternion operator*(const double& k, const Quaternion& q)
	{
		return Quaternion(k * q.m_w, k * q.m_v);
	}

	//Rotates 3D vector v by quaternion q to produce a new 3D vector.
	Vector3 operator*(const Quaternion& q, const Vector3& v)
	{
		//turn v to a quaternion p = [0, v]
		Quaternion p(0, v);

		//to rotate the vector to a new vector you do qpq^-1
		Quaternion res = q * p * inverse(q);

		Vector3 result(res.m_v);

		return result;
	}

	//Returns a Quaternion object that is the sum of q1 and q2.
	Quaternion operator+(const Quaternion& q1, const Quaternion& q2)
	{
		return Quaternion(q1.m_w + q2.m_w, q1.m_v + q2.m_v);
	}

	//Returns a Quaternion object that has the result of q1 - q2;
	Quaternion operator-(const Quaternion& q1, const Quaternion& q2)
	{
		return Quaternion(q1.m_w - q2.m_w, q1.m_v - q2.m_v);
	}

	//Return true if q1 and q2 are equal, false otherwise.
	bool operator==(const Quaternion& q1, const Quaternion& q2)
	{
		return FAUtil::compareDoubles(q1.m_w, q2.m_w, EPSILON) && FAUtil::compareDoubles(q1.x(), q2.x(), EPSILON) &&
			FAUtil::compareDoubles(q1.y(), q2.y(), EPSILON) && FAUtil::compareDoubles(q1.z(), q2.z(), EPSILON);
	}

	//Return true if q1 and q2 aren't equal, false otherwise.
	bool operator!=(const Quaternion& q1, const Quaternion& q2)
	{
		return !(q1 == q2);
	}

	//Returns the dot product between q1 and q2.
	double dotProduct(const Quaternion& q1, const Quaternion& q2)
	{
		return q1.m_w * q2.m_w + dotProduct(q1.m_v, q2.m_v);
	}

	//Spherically Interpolates between rotations q1 and q2.
	//t should be between 0 and 1. If t < 0 q1 will be returned. If t > 1 q2 will be returned.
	Quaternion slerp(const Quaternion& q1, const Quaternion& q2, const double& t)
	{
		if (t < 0.0)
			return q1;
		if (t > 1.0)
			return q2;

		//dot product of two quaternions is the cosine of the angle btw the two
		double c = dotProduct(q1, q2);

		//the dot product is negative so negate q1 to take the shorter arc
		if (c < 0.0f)
		{
			-q1;
		}

		//if the angle is very small meaning q1 and q2 are close then linearly interpolate
		double k0 = 0.0f;
		double k1 = 0.0f;
		if (c > 0.999f)
		{
			k0 = 1.0f - t;
			k1 = t;
		}
		else //the angle is not too small between q1 and q2
		{
			//sin^2(theta) = 1 - cos^2(theta)
			double s = sqrt(1.0f - c * c);

			//retrieve the angle
			double theta = atan2(s, c);

			double oneS = 1.0 / s;

			//Interploation parameters
			k0 = sin((1.0 - t) * theta) * oneS;
			k1 = sin(t * theta) * oneS;
		}

		//Slerp or lerp
		//Slerp - k0 = sin(1-t)theta / sin(theta), k1 = sin(t * theta) / sin(theta)
		//lerp - k0 = 1 - t, k1 = t
		//newQ = k0 * q1 + k1 * q2;
		Quaternion newQ = k0 * q1 + k1 * q2;

		return newQ;
	}

	void print(const Quaternion& q)
	{
		std::cout << "Scalar = " << std::setprecision(15) << q.m_w << std::endl;
		std::cout << "Vector = "; print(q.m_v);
		std::cout << std::endl;
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------

}
