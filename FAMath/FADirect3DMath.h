#pragma once

/** @file FADirect3DMAth.h
*	@brief File that has namespace FAMD3Dath. Withn the namespace are the classes Vector2D, Vector3D, Vector4D, Matrix4x4 and Quaternion.
*/

/**
*										FADIRECT3DMATH_H FILE
*
*/

/** @namespace FAD3DMath
*	@brief Has Vector2D, Vector3D, Vector4D, Matrix4x4, and Quaternion classes
*/

namespace FAD3DMath
{
	bool compareFloats(float x, float y, float epsilon);
	bool compareDoubles(double x, double y, double epsilon);

	/** @class Vector2D ""
	*	@brief A vector class used for 2D vectors/points and their manipulations.
	*
	*	The datatype for the components is float
	*/
	class Vector2D
	{
	public:


		/**@brief Default Constructor.
		*
		*	Creates a new 2D vector/point with the components initialized to 0.0.
		*/
		Vector2D();

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 2D vector/point with the components initialized to the arguments.
		*/
		Vector2D(float x, float y);

		/**@brief Returns the x component.
		*/
		float x() const;

		/**@brief Returns the y component.
		*/
		float y() const;

		/**@brief Sets the x component to the provided float.
		*/
		void setX(float x);

		/**@brief Sets the y component to the provided float.
		*/
		void setY(float y);


		/**@brief 2D vector addition through overloading operator +=.
		*/
		Vector2D& operator+=(const Vector2D& b);

		/**@brief 2D vector subtraction through overloading operator -=.
		*/
		Vector2D& operator-=(const Vector2D& b);

		/**@brief 2D vector scalar multiplication through overloading operator *=.
		*/
		Vector2D& operator*=(const float& k);

		/**@brief 2D vector scalar division through overloading operator /=.
		*
		* If k is zero, the vector is unchanged.
		*/
		Vector2D& operator/=(const float& k);

	private:
		float m_x;
		float m_y;
	};

	/**@brief Returns true if a is the zero vector.
	*/
	bool zeroVector(const Vector2D& a);

	/**@brief 2D vector addition.
	*/
	Vector2D operator+(const Vector2D& a, const Vector2D& b);

	/**@brief 2D vector negation.
	*/
	Vector2D operator-(const Vector2D& v);

	/**@brief 2D vector subtraction.
	*/
	Vector2D operator-(const Vector2D& a, const Vector2D& b);

	/**@brief 2D vector scalar multiplication.
	* Returns a * k, where a is a vector and k is a scalar(float)
	*/
	Vector2D operator*(const Vector2D& a, const float& k);

	/**@brief 2D vector scalar multiplication.
	* Returns k * a,  where a is a vector and k is a scalar(float)
	*/
	Vector2D operator*(const float& k, const Vector2D& a);

	/**@brief 2D vector scalar division.
	* Returns a / k,  where a is a vector and k is a scalar(float)
	* If k = 0 the returned vector is the zero vector.
	*/
	Vector2D operator/(const Vector2D& a, const float& k);

	/**@brief Returns the dot product between two 2D vectors.
	*
	*/
	float dotProduct(const Vector2D& a, const Vector2D& b);

	/**@brief Returns the length(magnitude) of the 2D vector v.
	*/
	float length(const Vector2D& v);

	/**@brief Normalizes the 2D vector v.
	*
	* If the 2D vector is the zero vector v is returned.
	*/
	Vector2D norm(const Vector2D& v);

	/**@brief Converts the 2D vector v from polar coordinates to cartesian coordinates.
	* v should = (r, theta(degrees))
	* The returned 2D vector = (x, y)
	*/
	Vector2D PolarToCartesian(const Vector2D& v);

	/**@brief Converts the 2D vector v from cartesian coordinates to polar coordinates.
	* v should = (x, y, z)
	* If vx is zero then no conversion happens and v is returned.\n
	* The returned 2D vector = (r, theta(degrees)).
	*/
	Vector2D CartesianToPolar(const Vector2D& v);

	/**@brief Returns a 2D vector that is the projection of a onto b.
	* If b is the zero vector a is returned.
	*/
	Vector2D Projection(const Vector2D& a, const Vector2D& b);

#if defined(_DEBUG)
	void print(const Vector2D& v);
#endif








	/** @class Vector3D ""
	*	@brief A vector class used for 3D vectors/points and their manipulations.
	*
	*	The datatype for the components is float
	*/
	class Vector3D
	{
	public:


		/**@brief Default Constructor.
		*
		*	Creates a new 3D vector/point with the components initialized to 0.0.
		*/
		Vector3D();

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 3D vector/point with the components initialized to the arguments.
		*/
		Vector3D(float x, float y, float z);

		/**@brief Returns the x component.
		*/
		float x() const;

		/**@brief Returns the y component.
		*/
		float y() const;

		/**@brief Returns the z component.
		*/
		float z() const;

		/**@brief Sets the x component to the provided float.
		*/
		void setX(float x);

		/**@brief Sets the y component to the provided float.
		*/
		void setY(float y);

		/**@brief Sets the z component to the provided float.
		*/
		void setZ(float z);


		/**@brief 3D vector addition through overloading operator +=.
		*/
		Vector3D& operator+=(const Vector3D& b);

		/**@brief 3D vector subtraction through overloading operator -=.
		*/
		Vector3D& operator-=(const Vector3D& b);

		/**@brief 3D vector scalar multiplication through overloading operator *=.
		*/
		Vector3D& operator*=(const float& k);

		/**@brief 3D vector scalar division through overloading operator /=.
		*
		* If k is zero, the vector is unchanged.
		*/
		Vector3D& operator/=(const float& k);

	private:
		float m_x;
		float m_y;
		float m_z;
	};

	/**@brief Returns true if a is the zero vector.
	*/
	bool zeroVector(const Vector3D& a);

	/**@brief 3D vector addition.
	*/
	Vector3D operator+(const Vector3D& a, const Vector3D& b);

	/**@brief 3D vector negeation.
	*/
	Vector3D operator-(const Vector3D& v);

	/**@brief 3D vector subtraction.
	*/
	Vector3D operator-(const Vector3D& a, const Vector3D& b);

	/**@brief 3D vector scalar multiplication.
	* Returns a * k, where a is a vector and k is a scalar(float)
	*/
	Vector3D operator*(const Vector3D& a, const float& k);

	/**@brief 3D vector scalar multiplication.
	* Returns k * a,  where a is a vector and k is a scalar(float)
	*/
	Vector3D operator*(const float& k, const Vector3D& a);

	/**@brief 3D vector scalar division.
	* Returns a / k,  where a is a vector and k is a scalar(float)
	* If k = 0 the returned vector is the zero vector.
	*/
	Vector3D operator/(const Vector3D& a, const float& k);

	/**@brief Returns the dot product between two 3D vectors.
	*/
	float dotProduct(const Vector3D& a, const Vector3D& b);

	/**@brief Returns the cross product between two 3D vectors.
	*/
	Vector3D crossProduct(const Vector3D& a, const Vector3D& b);

	/**@brief Returns the length(magnitude) of the 3D vector v.
	*/
	float length(const Vector3D& v);

	/**@brief Normalizes the 3D vector v.
	*
	* If the 3D vector is the zero vector v is returned.
	*/
	Vector3D norm(const Vector3D& v);

	/**@brief Converts the 3D vector v from cylindrical coordinates to cartesian coordinates.
	* v should = (r, theta(degrees), z).\n
	* The returned 3D vector = (x, y ,z).
	*/
	Vector3D CylindricalToCartesian(const Vector3D& v);

	/**@brief Converts the 3D vector v from cartesian coordinates to cylindrical coordinates.
	* v should = (x, y, z).\n
	* If vx is zero then no conversion happens and v is returned.\n
	* The returned 3D vector = (r, theta(degrees), z).
	*/
	Vector3D CartesianToCylindrical(const Vector3D& v);

	/**@brief Converts the 3D vector v from spherical coordinates to cartesian coordinates.
	* v should = (pho, phi(degrees), theta(degrees)).\n
	* The returned 3D vector = (x, y, z)
	*/
	Vector3D SphericalToCartesian(const Vector3D& v);

	/**@brief Converts the 3D vector v from cartesian coordinates to spherical coordinates.
	* If v is the zero vector or if vx is zero then no conversion happens and v is returned.\n
	* The returned 3D vector = (r, phi(degrees), theta(degrees)).
	*/
	Vector3D CartesianToSpherical(const Vector3D& v);

	/**@brief Returns a 3D vector that is the projection of a onto b.
	* If b is the zero vector a is returned.
	*/
	Vector3D Projection(const Vector3D& a, const Vector3D& b);


#if defined(_DEBUG)
	void print(const Vector3D& v);
#endif










	/** @class Vector4D ""
	*	@brief A vector class used for 4D vectors/points and their manipulations.
	*
	*	The datatype for the components is float
	*/
	class Vector4D
	{
	public:


		/**@brief Default Constructor.
		*
		*	Creates a new 4D vector/point with the components initialized to 0.0.
		*/
		Vector4D();

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector/point with the components initialized to the arguments.
		*/
		Vector4D(float x, float y, float z, float w);

		/**@brief Returns the x component.
		*/
		float x() const;

		/**@brief Returns the y component.
		*/
		float y() const;

		/**@brief Returns the z component.
		*/
		float z() const;

		/**@brief Returns the w component.
		*/
		float w() const;

		/**@brief Sets the x component to the provided float.
		*/
		void setX(float x);

		/**@brief Sets the y component to the provided float.
		*/
		void setY(float y);

		/**@brief Sets the z component to the provided float.
		*/
		void setZ(float z);

		/**@brief Sets the w component to the provided float.
		*/
		void setW(float w);


		/**@brief 4D vector addition through overloading operator +=.
		*/
		Vector4D& operator+=(const Vector4D& b);

		/**@brief 4D vector subtraction through overloading operator -=.
		*/
		Vector4D& operator-=(const Vector4D& b);

		/**@brief 4D vector scalar multiplication through overloading operator *=.
		*/
		Vector4D& operator*=(const float& k);

		/**@brief 4D vector scalar division through overloading operator /=.
		*
		* If k is zero, the vector is unchanged.
		*/
		Vector4D& operator/=(const float& k);

	private:
		float m_x;
		float m_y;
		float m_z;
		float m_w;
	};

	/**@brief Returns true if a is the zero vector.
	*/
	bool zeroVector(const Vector4D& a);

	/**@brief 4D vector addition.
	*/
	Vector4D operator+(const Vector4D& a, const Vector4D& b);

	/**@brief 4D vector negation.
	*/
	Vector4D operator-(const Vector4D& v);

	/**@brief 4D vector subtraction.
	*/
	Vector4D operator-(const Vector4D& a, const Vector4D& b);

	/**@brief 4D vector scalar multiplication.
	* Returns a * k, where a is a vector and k is a scalar(float)
	*/
	Vector4D operator*(const Vector4D& a, const float& k);

	/**@brief 4D vector scalar multiplication.
	* Returns k * a,  where a is a vector and k is a scalar(float)
	*/
	Vector4D operator*(const float& k, const Vector4D& a);

	/**@brief 4D vector scalar division.
	* Returns a / k,  where a is a vector and k is a scalar(float)
	* If k = 0 the returned vector is the zero vector.
	*/
	Vector4D operator/(const Vector4D& a, const float& k);

	/**@brief Returns the dot product between two 4D vectors.
	*/
	float dotProduct(const Vector4D& a, const Vector4D& b);

	/**@brief Returns the length(magnitude) of the 4D vector v.
	*/
	float length(const Vector4D& v);

	/**@brief Normalizes the 4D vector v.
	*
	* If the 4D vector is the zero vector v is returned.
	*/
	Vector4D norm(const Vector4D& v);

	/**@brief Returns a 4D vector that is the projection of a onto b.
	* If b is the zero vector a is returned.
	*/
	Vector4D Projection(const Vector4D& a, const Vector4D& b);


#if defined(_DEBUG)
	void print(const Vector4D& v);
#endif



	/** @class Matrix4x4 ""
	*	@brief A matrix class used for 4x4 matrices and their manipulations.
	*
	*	The datatype for the components is float.\n
	*	The 4x4 matrix is treated as a row-major matrix.
	*
	*/
	class Matrix4x4
	{
	public:

		/**@brief Default Constructor.
		*
		*	Creates a new 4x4 identity matrix.
		*/
		Matrix4x4();

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4x4 matrix with elements initialized to the given 2D array.\n
		*	If the passed in 2D array isn't a 4x4 matrix, the behavior is undefined.
		*/
		Matrix4x4(float a[][4]);

		/**@brief Returns the element at the given (row, col).
		* The row and col values should be between [0,3]. If any of them are out of that range, the first element will be returned.
		*/
		float getElement(unsigned int row, unsigned int col) const;

		/**@brief Returns a constant reference to the element at the given (row, col).
		* The row and col values should be between [0,3]. If any of them are out of that range, the first element will be returned.
		*/
		const float& operator()(unsigned int row, unsigned int col) const;

		/**@brief Returns a reference to the element at the given (row, col).
		* The row and col values should be between [0,3]. If any of them are out of that range, the first element will be returned.
		*/
		float& operator()(unsigned int row, unsigned int col);

		/**@brief Set the element at the given (row, col) to the given val.
		* The row and col values should be between [0,3]. If any of them are out of that range, the first element will be set to the given val.
		*/
		void setElement(unsigned int row, unsigned int col, float val);


		/**@brief Sets the matrix to the identity matrix.
		*/
		void setToIdentity();

		/**@brief Returns true if the matrix is the identity matirx, false otherwise.
		*/
		bool isIdentity();

		/**@brief Adds this 4x4 matrix with given matrix m and stores the result in this 4x4 matrix.
		*/
		Matrix4x4& operator+=(const Matrix4x4& m);

		/**@brief Subtracts this 4x4 matrix with given matrix m and stores the result in this 4x4 matrix.
		*/
		Matrix4x4& operator-=(const Matrix4x4& m);

		/**@brief Multiplies this 4x4 matrix with given scalar k and stores the result in this 4x4 matrix.
		*/
		Matrix4x4& operator*=(const float& k);


		/**@brief Multiplies this 4x4 matrix with given matrix m and stores the result in this 4x4 matrix.
		*/
		Matrix4x4& operator*=(const Matrix4x4& m);



	private:

		float m_mat[4][4];
	};

	/**@brief Adds the two given 4x4 matrices and returns a Matrix4x4 object with the result.
	*/
	Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2);

	/**@brief Negates the 4x4 matrix m.
	*/
	Matrix4x4 operator-(const Matrix4x4& m);

	/**@brief Subtracts the two given 4x4 matrices and returns a Matrix4x4 object with the result.
	*/
	Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2);

	/**@brief Multiplies the given 4x4 matrix with the given scalar and returns a Matrix4x4 object with the result.
	*/
	Matrix4x4 operator*(const Matrix4x4& m, const float& k);

	/**@brief Multiplies the the given scalar with the given 4x4 matrix and returns a Matrix4x4 object with the result.
	*/
	Matrix4x4 operator*(const float& k, const Matrix4x4& m);

	/**@brief Multiplies the two given 4x4 matrices and returns a Matrix4x4 object with the result.
	*/
	Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2);

	/**@brief Multiplies the given 4x4 matrix with the given 4D vector and returns a Vector4D object with the result.
	* The vector is a column vector.
	*/
	Vector4D operator*(const Matrix4x4& m, const Vector4D& v);

	/**@brief Multiplies the given 4D vector with the given 4x4 matrix and returns a Vector4D object with the result.
	* The vector is a row vector.
	*/
	Vector4D operator*(const Vector4D& a, const Matrix4x4& m);

	/**@brief Returns the tranpose of the given matrix m.
	*/
	Matrix4x4 transpose(const Matrix4x4& m);

	/**@brief Construct a 4x4 translation matrix with the given floats and post-multiply it by the given matrix.
	* cm = cm * translate
	*/
	Matrix4x4 translate(const Matrix4x4& cm, float x, float y, float z);

	/**@brief Construct a 4x4 scaling matrix with the given floats and post-multiply it by the given matrix.
	* cm = cm * scale
	*/
	Matrix4x4 scale(const Matrix4x4& cm, float x, float y, float z);

	/**@brief Construct a 4x4 rotation matrix with the given angle (in degrees) and axis (x, y, z) and post-multiply it by the given matrix.
	* cm = cm * rotate.\n
	*/
	Matrix4x4 rotate(const Matrix4x4& cm, float angle, float x, float y, float z);

	/**@brief Returns the determinant of the given matrix.
	*/
	double det(const Matrix4x4& m);

	/**@brief Returns the cofactor of the given row and col using the given matrix.
	*/
	double cofactor(const Matrix4x4& m, unsigned int row, unsigned int col);

	/**@brief Returns the adjoint of the given matrix.
	*/
	Matrix4x4 adjoint(const Matrix4x4& m);

	/**@brief Returns the inverse of the given matrix.
	* If the matrix is noninvertible/singular, the identity matrix is returned.
	*/
	Matrix4x4 inverse(const Matrix4x4& m);



#if defined(_DEBUG)
	void print(const Matrix4x4& m);
#endif















	/** @class Quaternion ""
	*
	*	The datatype for the components is float.\n
	*
	*/
	class Quaternion
	{
	public:

		/**@brief Default Constructor.
		*
		*	Constructs an identity quaternion.
		*/
		Quaternion();

		/**@brief Overloaded Constructor.
		*
		*	Constructs a quaternion with the given values.
		*/
		Quaternion(float scalar, float x, float y, float z);

		/**@brief Overloaded Constructor.
		*
		*	Constructs a quaternion with the given values.
		*/
		Quaternion(float scalar, const Vector3D& v);

		/**@brief Overloaded Constructor.
		*
		*	Constructs a quaternion with the given values in the 4D vector.\n
		*	The x value in the 4D vector should be the scalar.
		*	The y, z and w value in the 4D vector should be the axis.
		*/
		Quaternion(const Vector4D& v);

		/**@brief Returns a reference to the scalar component of the quaternion.
		*/
		float& scalar();

		/**@brief Returns a const reference to the scalar component of the quaternion.
		*/
		const float& scalar() const;

		/**@brief Returns a reference to the x value of the vector component in the quaternion.
		*/
		float& x();

		/**@brief Returns a const reference to the x value of the vector component in the quaternion.
		*/
		const float& x() const;

		/**@brief Returns a reference to the y value of the vector component in the quaternion.
		*/
		float& y();

		/**@brief Returns a const reference to the y value of the vector component in the quaternion.
		*/
		const float& y() const;

		/**@brief Returns a reference to the z value of the vector component in the quaternion.
		*/
		float& z();

		/**@brief Returns a const reference to the z value of the vector component in the quaternion.
		*/
		const float& z() const;

		/**@brief Returns the vector component of the quaternion.
		*/
		Vector3D vector();

		/**@brief Adds this quaternion to quaterion q and stores the result in this quaternion.
		*/
		Quaternion& operator+=(const Quaternion& q);

		/**@brief Subtracts this quaternion by quaterion q and stores the result in this quaternion.
		*/
		Quaternion& operator-=(const Quaternion& q);

		/**@brief Multiplies this quaternion by flaot k and stores the result in this quaternion.
		*/
		Quaternion& operator*=(float k);

		/**@brief Multiplies this quaternion by quaterion q and stores the result in this quaternion.
		*/
		Quaternion& operator*=(const Quaternion& q);


	private:

		float m_scalar;
		float m_x;
		float m_y;
		float m_z;
	};

	/**@brief Returns a quaternion that has the result of q1 + q2.
	*/
	Quaternion operator+(const Quaternion& q1, const Quaternion& q2);

	/**@brief Returns a quaternion that has the result of -q.
	*/
	Quaternion operator-(const Quaternion& q);

	/**@brief Returns a quaternion that has the result of q1 - q2.
	*/
	Quaternion operator-(const Quaternion& q1, const Quaternion& q2);

	/**@brief Returns a quaternion that has the result of k * q.
	*/
	Quaternion operator*(float k, const Quaternion& q);

	/**@brief Returns a quaternion that has the result of q * k.
	*/
	Quaternion operator*(const Quaternion& q, float k);

	/**@brief Returns a quaternion that has the result of q1 * q2.
	*/
	Quaternion operator*(const Quaternion& q1, const Quaternion& q2);


	/**@brief Returns true if quaternion q is a zero quaternion, false otherwise.
	*/
	bool isZeroQuaternion(const Quaternion& q);

	/**@brief Returns true if quaternion q is an identity quaternion, false otherwise.
	*/
	bool isIdentity(const Quaternion& q);

	/**@brief Returns the conjugate of quaternion q.
	*/
	Quaternion conjugate(const Quaternion& q);

	/**@brief Returns the length of quaternion q.
	*/
	float length(const Quaternion& q);

	/**@brief Normalizes quaternion q and returns the normalized quaternion.
	* If q is the zero quaternion then q is returned.
	*/
	Quaternion normalize(const Quaternion& q);

	/**@brief Returns the invese of quaternion q.
	*  If q is the zero quaternion then q is returned.
	*/
	Quaternion inverse(const Quaternion& q);

	/**@brief Returns a quaternion from the axis-angle rotation representation.
	*  The angle should be given in degrees.
	*/
	Quaternion rotationQuaternion(float angle, float x, float y, float z);

	/**@brief Returns a quaternion from the axis-angle rotation representation.
	* The angle should be given in degrees.
	*/
	Quaternion rotationQuaternion(float angle, const Vector3D& axis);

	/**@brief Returns a quaternion from the axis-angle rotation representation.
	* The x value in the 4D vector should be the angle(in degrees).\n
	* The y, z and w value in the 4D vector should be the axis.
	*/
	Quaternion rotationQuaternion(const Vector4D& angAxis);

	/**@brief Returns a matrix from the given quaterion for column vector-matrix multiplication.
	* Quaternion q should be a unit quaternion.
	*/
	Matrix4x4 quaternionRotationMatrixCol(const Quaternion& q);

	/**@brief Returns a matrix from the given quaterion for row vector-matrix multiplication.
	* Quaternion q should be a unit quaternion.
	*/
	Matrix4x4 quaternionRotationMatrixRow(const Quaternion& q);

#if defined(_DEBUG)
	void print(const Quaternion& q);
#endif

}
