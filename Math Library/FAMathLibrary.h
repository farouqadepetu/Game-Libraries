/** @file MathLibrary.h
*	@brief File that has namespace FAMath. Withn the namespace are the classes Vector2, Vector3, Vector4, Matrix4x4 and Quaternion.
*/

/**
*										FAMATH_H FILE
*										AUTHOR: FAROUQ ADEPETU
*/

#pragma once

/** @namespace FAMath
*	@brief Has Vector2, Vector3, Vector4, and Matrix4x4 classes
*/
namespace FAMath
{
	class Vector2;
	class Vector3;
	class Vector4;

	/** @class Vector2 ""
	*	@brief A vector class used for 2D vectors/points and their manipulations.
	*	
	*	The datatype for the components is double
	*/
	class Vector2
	{
	public:

		/** @name Constructors
		*	Constructors for class FAMath::Vector2
		*/
		///@{

		/**@brief Default Constructor.
		* 
		*	Creates a new 2D vector/point with the components initialized to 0.0.
		*	
		*	@see Vector2(double x, double y);
		*/
		Vector2();

		/**@brief Overloaded Constructor.
		* 
		*	Creates a new 2D vector/point with the components initialized to the arguments.
		*/
		Vector2(double x, double y);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 2D vector/point using v's x and y coordinates.
		*/
		Vector2(const Vector3& v);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 2D vector/point using v's x and y coordinates
		*
		*/
		Vector2(const Vector4& v);

		///@}

		/** @name Getters
		*	Getters for class FAMath::Vector2
		*/
		///@{

		/**@brief Returns the value of the x-coordinate.
		* 
		*   
		*	@see double y()
		*/
		double x() const;

		/**@brief Returns the value of the y-coordinate.
		*	
		*	@see double x()
		*/
		double y() const;

		///@}

		/** @name Setters
		*	Setters for class FAMath::Vector2.
		*/
		///@{

		/**@brief Sets the x-coordinate of the 2D vector/point.
		* 
		*	@see void setY(double y)
		*/
		void setX(double x);

		/**@brief Sets the y-coordinate of the 2D vector/point.
		* 
		*	@see void setX(double x)
		*/
		void setY(double y);

		///@}

		/**@brief Returns true if all the components of the 2D vector equals to zero, false otherwise.
		*/
		bool isZeroVector() const;

		/**@brief 2D vector addition through operator += overloading.
		*
		*	@return A reference to the current Vector2 object\n
		*	That has the result of the current Vector2 object + Vector2 object b.
		*/
		Vector2& operator+=(const Vector2& b);

		/**@brief 2D vector subtraction through operator -= overloading.
		*
		*	@return A reference to the current Vector2 object\n
		*	That has the result of the current Vector2 object - Vector2 object b.
		*/
		Vector2& operator-=(const Vector2& b);

		/**@brief 2D vector multiplication by a scalar through operator *= overloading.
		*
		*	@return A reference to the current vector object\n
		*	That has the result of the current Vector2 object * scalar.
		*/
		Vector2& operator*=(const double& scalar);

		/**@brief 2D vector division by a scalar through operator /= overloading.
		*	
		*	Throws a std::invalid_argument if scalar is zero
		* 
		*	@return A reference to the current vector object\n
		*	That has the result of the current Vector2 object / scalar.
		*/
		Vector2& operator/=(const double& scalar);

		/**@brief Assignment Operator Overloading
		*
		*	Stores the x and y values of Vector3s v into the x and y values of this Vector2.
		*/
		void operator=(const Vector3& v);

		/**@brief Assignment Operator Overloading
		*
		*	Stores the x and y values of Vector4s v into the x and y values of this Vector2.
		*/
		void operator=(const Vector4& v);

		/**@brief 2D vector addition through operator + overloading.
		*	
		*	@return A new Vector2 object that has the result of a + b.
		*/
		friend Vector2 operator+(const Vector2& a, const Vector2& b);

		/**@brief 2D vector subtraction through operator - overloading.
		*
		*	@return A new Vector2 object that has the result of a - b.
		*/
		friend Vector2 operator-(const Vector2& a, const Vector2& b);

		/**@brief 2D vector multiplication by a scalar through operator * overloading.
		*	
		*	Called when you do v * scalar, where v is a Vector2 object.
		* 
		*	@return A new Vector2 object that has the result of v * scalar.
		*/
		friend Vector2 operator*(const Vector2& v, const double& scalar);

		/**@brief 2D vector multiplication by a scalar through operator * overloading.
		*
		*	Called when you do scalar * v, where v is a Vector2 object
		*
		*	@return A new Vector2 object that has the result of scalar * v.
		*/
		friend Vector2 operator*(const double& scalar, const Vector2& v);

		/**@brief 2D vector divison by a scalar through operator / overloading.
		*
		*	@return A new Vector2 object that has the result of v/scalar.
		*/
		friend Vector2 operator/(const Vector2& v, const double& scalar);

		/**@brief Compares two 2D vectors through operator == overloading.
		*
		*	@return True if a equals to b.
		*	@return False otherwise.
		*/
		friend bool operator==(const Vector2& a, const Vector2& b);

		/**@brief Compares two 2D vectors through operator != overloading.
		*
		*	@return False if a equals to b.
		*	@return True otherwise.
		*/
		friend bool operator!=(const Vector2& a, const Vector2& b);
		
		/**@brief Magnitude/Length of a 2D vector.
		* 
		*	@return The length of Vector2 object v.
		*/
		friend double length(const Vector2& v);

		/**@brief Normalizes a 2D vector
		*
		*	If v is a zero vector then the zero vector is returned.
		* 
		*	@return A Vector2 object that has the result v / |v| which is a unit vector.
		*/
		friend Vector2 normalize(const Vector2& v);

		/**@brief Distane between two 2D points.
		* 
		*	@return The distance between two 2D points.
		*/
		friend double distance(const Vector2& a, const Vector2& b);

		/**@brief Dot Product.
		* 
		*	@return The value of a dot b.	
		*/
		friend double dotProduct(const Vector2& a, const Vector2& b);

		/**@brief Angle between two 2D vectors.
		*
		*	a and b should be unit vectors before using them as arguments.
		* 
		*	@return The angle between two 2D vectors.
		*/
		friend double angle(const Vector2& a, const Vector2& b);

		friend void print(Vector2 v);

	private:
		//components of a 2D Vector
		double m_x;
		double m_y;
	};










	/** @class Vector3 ""
	*	@brief A vector class used for 3D vectors/points and their manipulations.
	*
	*	The datatype for the components is double.
	*/
	class Vector3
	{
	public:

		/** @name Constructors
		*	Constructors for class FAMath::Vector3.
		*/
		///@{

		/**@brief Default Constructor.
		*
		*	Creates a new 3D vector/point with the components initialized to 0.0.
		* 
		*	@see Vector3(double x, double y);
		*/
		Vector3();

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 3D vector/point with the components initialized to the arguments.
		*/
		Vector3(double x, double y, double z);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 3D vector/point using v's x and y coordinates and sets z to 0.0.
		*/
		Vector3(const Vector2& v);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 3D vector/point using v's x and y coordinates and sets this Vector3 z component to the given z
		*/
		Vector3(const Vector2& v, const double& z);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 3D vector/point using v's x, y and z coordinates.
		*/
		Vector3(const Vector4& v);

		///@}

		/** @name Getters
		*	Getters for class FAMath::Vector3.
		*/
		///@{

		/**@brief Returns the value of the x-coordinate.
		*
		*	@see double y()
		*	@see double z()
		*/
		double x() const;

		/**@brief Returns the value of the y-coordinate.
		*
		*	@see double x()
		*	@see double z()
		*/
		double y() const;

		/**@brief Returns the value of the z-coordinate.
		*
		*	@see double x()
		*	@see double y()
		*/
		double z() const;

		///@}

		/** @name Setters
		*	Setters for class FAMath::Vector3.
		*/
		///@{

		/**@brief Sets the x-coordinate of the 3D vector/point.
		* 
		*	@see void setY(double y)
		*	@see void setZ(double z)
		*/
		void setX(double x);

		/**@brief Sets the y-coordinate of the 3D vector/point.
		* 
		*	@see void setX(double x)
		*	@see void setZ(double z)
		*/
		void setY(double y);

		/**@brief Sets the z-coordinate of the 3D vector/point.
		*
		*	@see void setX(double x)
		*	@see void setY(double y)
		*/
		void setZ(double z);

		///@}

		/**@brief Returns true if all the components of the 3D vector equals to zero, false otherwise.
		*/
		bool isZeroVector() const;

		/**@brief 3D vector addition through operator += overloading.
		*
		*	@return A reference to the current Vector3 object\n
		*	That has the result of the current Vector3 object + Vector3 object b.
		*/
		Vector3& operator+=(const Vector3& b);

		/**@brief 3D vector subtraction through operator -= overloading.
		*
		*	@return A reference to the current Vector3 object\n
		*	That has the result of the current Vector3 object - Vector3 object b.
		*/
		Vector3& operator-=(const Vector3& b);

		/**@brief 3D vector multiplication by a scalar through operator *= overloading.
		*
		*	@return A reference to the current vector object\n
		*	That has the result of the current Vector3 object * scalar.
		*/
		Vector3& operator*=(const double& scalar);

		/**@brief 3D vector division by a scalar through operator /= overloading.
		*
		*	Throws a std::invalid_argument if scalar is zero
		*
		*	@return A reference to the current vector object\n
		*	That has the result of the current Vector3 object / scalar.
		*/
		Vector3& operator/=(const double& scalar);

		/**@brief Assignment Operator Overloading
		*
		*	Stores the x and y values of Vector2s v into the x and y values of this Vector3 and sets z = 0.0.
		*/
		void operator=(const Vector2& v);

		/**@brief Assignment Operator Overloading
		*
		*	Stores the x, y and z values of Vector4s v into the x, y and z values of this Vector3.
		*/
		void operator=(const Vector4& v);

		/**@brief 3D vector addition through operator + overloading.
		*
		*	@return A new Vector3 object that has the result of a + b.
		*/
		friend Vector3 operator+(const Vector3& a, const Vector3& b);

		/**@brief 3D vector subtraction through operator - overloading.
		*
		*	@return A new Vector3 object that has the result of a - b.
		*/
		friend Vector3 operator-(const Vector3& a, const Vector3& b);

		/**@brief 3D vector multiplication by a scalar through operator * overloading.
		*
		*	Called when you do v * scalar, where v is a Vector3 object.
		*
		*	@return A new Vector3 object that has the result of v * scalar.
		*/
		friend Vector3 operator*(const Vector3& v, const double& scalar);

		/**@brief 3D vector multiplication by a scalar through operator * overloading.
		*
		*	Called when you do scalar * v, where v is a Vector3 object.
		*
		*	@return A new Vector3 object that has the result of scalar * v.
		*/
		friend Vector3 operator*(const double& scalar, const Vector3& v);

		/**@brief 3D vector divison by a scalar through operator / overloading.
		*
		*	@return A new Vector3 object that has the result of v/scalar.
		*/
		friend Vector3 operator/(const Vector3& v, const double& scalar);

		/**@brief Compares two 3D vectors through operator == overloading.
		*
		*	@return True if a equals to b.
		*	@return False otherwise.
		*/
		friend bool operator==(const Vector3& a, const Vector3& b);

		/**@brief Compares two 3D vectors through operator != overloading.
		*
		*	@return False if a equals to b.
		*	@return True otherwise.
		*/
		friend bool operator!=(const Vector3& a, const Vector3& b);

		/**@brief Magnitude/Length of a 3D vector.
		*
		*	@return The length of Vector3 object v.
		*/
		friend double length(const Vector3& v);

		/**@brief Normalizes a 3D vector
		*
		*	If v is a zero vector then the zero vector is returned.
		* 
		*	@return A Vector3 object that is a unit vector(has a length of 1).
		*/
		friend Vector3 normalize(const Vector3& v);

		/**@brief Distane between two 3D points.
		* 
		*	@return The distance between two 3D points.
		*/
		friend double distance(const Vector3& a, const Vector3& b);

		/**@brief Dot Product.
		*
		*	@return The value of a dot b.
		*/
		friend double dotProduct(const Vector3& a, const Vector3& b);

		/**@brief Angle between two 3D vectors.
		*
		*	a and b should be unit vectors before using them as arguments.
		*
		*	@return The angle between two 3D vectors.
		*/
		friend double angle(const Vector3& a, const Vector3& b);

		/**@brief Cross Product.
		*
		*	@return A Vector3 object that is perpendicular to a and b.
		*/
		friend Vector3 crossProduct(const Vector3& a, const Vector3& b);

		friend void print(Vector3 v);

	private:
		//components of a 3D Vector
		double m_x;
		double m_y;
		double m_z;
	};













	/** @class Vector4 ""
	*	@brief A vector class used for 4D vectors/points and their manipulations.
	*
	*	The datatype for the components is double
	*/
	class Vector4
	{
	public:

		/** @name Constructors
		*	Constructors for class FAMath::Vector4,
		*/
		///@{

		/**@brief Default Constructor.
		*
		*	Creates a new 4D vector/point with the components initialized to 0.0.
		*
		*	@see Vector4(double x, double y);
		*/
		Vector4();

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector/point with the components initialized to the arguments,
		*/
		Vector4(double x, double y, double z, double w);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector/point using v's x and y coordinates and sets z and w to 0.0.
		*/
		Vector4(const Vector2& v);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector/point using v's x and y coordinates and z and w values.
		*/
		Vector4(const Vector2& v, const double& z, const double& w);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector / point using v's x and y coordinates and the given z and sets w to 0.0.
		*/
		Vector4(const Vector2& v, const double& z);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector/point using v's x, y and z coordinates and sets w to 0.0.
		*/
		Vector4(const Vector3& v);

		/**@brief Overloaded Constructor.
		*
		*	Creates a new 4D vector/point using v's x, y and z coordinates and the given w value.
		*/
		Vector4(const Vector3& v, const double& w);


		///@}

		/** @name Getters
		*	Getters for class FAMath::Vector4.
		*/
		///@{

		/**@brief Returns the value of the x-coordinate.
		*	@see double y()
		*	@see double z()
		*	@see double w()
		*/
		double x() const;

		/**@brief Returns the value of the y-coordinate.
		*
		*	@see double x()
		*	@see double z()
		*	@see double w()
		*/
		double y() const;

		/**@brief Returns the value of the z-coordinate.
		*
		*	@see double x()
		*	@see double y()
		*	@see double w()
		*/
		double z() const;

		/**@brief Returns the value of the w-coordinate.
		*
		*	@see double x()
		*	@see double y()
		*	@see double z()
		*/
		double w() const;

		///@}

		/** @name Setters
		*	Setters for class FAMath::Vector4.
		*/
		///@{

		/**@brief Sets the x-coordinate of the 4D vector/point.
		* 
		*	@see void setY(double y)
		*	@see void setZ(double z)
		*/
		void setX(double x);

		/**@brief Sets the y-coordinate of the 4D vector/point.
		*
		*	@see void setX(double x)
		*	@see void setZ(double z)
		*/
		void setY(double y);

		/**@brief Sets the z-coordinate of the 4D vector/point.
		*
		*	@see void setX(double x)
		*	@see void setY(double y)
		*/
		void setZ(double z);

		/**@brief Sets the w-coordinate of the 4D vector/point.
		*
		*	@see void setX(double x)
		*	@see void setY(double y)
		*/
		void setW(double w);

		///@}

		/**@brief Returns true if all the components of the 4D vector equals to zero, false otherwise.
		*/
		bool isZeroVector() const;

		/**@brief 4D vector addition through operator += overloading.
		*
		*	@return A reference to the current Vector4 object\n
		*	That has the result of the current Vector4 object + Vector4 object b.
		*/
		Vector4& operator+=(const Vector4& b);

		/**@brief 4D vector subtraction through operator -= overloading.
		*
		*	@return A reference to the current Vector4 object\n
		*	That has the result of the current Vector4 object - Vector4 object b.
		*/
		Vector4& operator-=(const Vector4& b);

		/**@brief 4D vector multiplication by a scalar through operator *= overloading.
		*
		*	@return A reference to the current vector object\n
		*	That has the result of the current Vector4 object * scalar.
		*/
		Vector4& operator*=(const double& scalar);

		/**@brief 4D vector division by a scalar through operator /= overloading.
		*
		*	Throws a std::invalid_argument if scalar is zero
		*
		*	@return A reference to the current vector object\n
		*	That has the result of the current Vector4 object / scalar.
		*/
		Vector4& operator/=(const double& scalar);

		/**@brief Assignment Operator Overloading
		*
		*	Stores the x and y values of Vector2s v into the x and y values of this Vector4 and sets w and to 0.0.
		*/
		void operator=(const Vector2& v);

		/**@brief Assignment Operator Overloading
		*
		*	Stores the x, y and z values of Vector3s v into the x, y and z values of this Vector4 and sets w to 0.0.
		*/
		void operator=(const Vector3& v);

		/**@brief 4D vector addition through operator + overloading.
		*
		*	@return A new Vector4 object that has the result of a + b.
		*/
		friend Vector4 operator+(const Vector4& a, const Vector4& b);

		/**@brief 4D vector subtraction through operator - overloading.
		*
		*	@return A new Vector4 object that has the result of a - b.
		*/
		friend Vector4 operator-(const Vector4& a, const Vector4& b);

		/**@brief 4D vector multiplication by a scalar through operator * overloading.
		*
		*	Called when you do v * scalar, where v is a Vector4 object.
		*
		*	@return A new Vector4 object that has the result of v * scalar.
		*/
		friend Vector4 operator*(const Vector4& v, const double& scalar);

		/**@brief 4D vector multiplication by a scalar through operator * overloading.
		*
		*	Called when you do scalar * v, where v is a Vector4 object.
		*
		*	@return A new Vector4 object that has the result of scalar * v.
		*/
		friend Vector4 operator*(const double& scalar, const Vector4& v);

		/**@brief 4D vector divison by a scalar through operator / overloading.
		*
		*	@return A new Vector4 object that has the result of v/scalar.
		*/
		friend Vector4 operator/(const Vector4& v, const double& scalar);

		/**@brief Compares two 4D vectors through operator == overloading.
		*
		*	@return True if a equals to b.
		*	@return False otherwise.
		*/
		friend bool operator==(const Vector4& a, const Vector4& b);

		/**@brief Compares two 4D vectors through operator != overloading.
		*
		*	@return False if a equals to b.
		*	@return True otherwise.
		*/
		friend bool operator!=(const Vector4& a, const Vector4& b);

		/**@brief Magnitude/Length of a 4D vector.
		*
		*	@return The length of Vector4 object v.
		*/
		friend double length(const Vector4& v);

		/**@brief Normalizes a 4D vector
		*
		*	@return A Vector4 object that is a unit vector(has a length of 1).
		*/
		friend Vector4 normalize(const Vector4& v);

		/**@brief Distane between two 4D points.
		*
		*	If v is the zero vector then the zero vector is returned.
		*	@return The distance between two 4D points.
		*/
		friend double distance(const Vector4& a, const Vector4& b);

		/**@brief Dot Product.
		*
		*	@return The value of a dot b.
		*
		*/
		friend double dotProduct(const Vector4& a, const Vector4& b);

		/**@brief Angle between two 4D vectors.
		* 
		*	a and b should be unit vectors before using them as arguments.
		*
		*	@return The angle between two 4D vectors.
		*
		*/
		friend double angle(const Vector4& a, const Vector4& b);

		friend void print(Vector4 v);

	private:
		//components of a 4D Vector
		double m_x;
		double m_y;
		double m_z;
		double m_w;
	};















	/** @class Matrix4x4 ""
	*	@brief A matrix class used for 4x4 matrices and their manipulations.
	*
	*	The datatype for the matrix is double
	* 
	*	The 4x4 matrix is treated as a row-major matrix\n
	*	the constructors and other functions take in the data as row major format
	*	
	*	Internally the data is stored in column-major order
	*/
	class Matrix4x4
	{
	public:

		/** @name Constructors
		*	Constructors for class FAMath::Matrix4x4.
		*/
		///@{

		/**@brief Default Constructor.
		* 
		*	Constructs an identity matrix.
		*/
		Matrix4x4();

		/**@brief Overloaded Constructor.
		* 
		*	Constructors a 4x4 matrix from the given std::array.
		*/
		Matrix4x4(const double* values);


		/**@brief Overloaded Constructor.
		*
		*	Constructs a 4x4 matrix from the specified 16 elements.\n
		*	The elements are specified in row-major order.
		*/
		Matrix4x4(double m11, double m12, double m13, double m14,
			double m21, double m22, double m23, double m24,
			double m31, double m32, double m33, double m34,
			double m41, double m42, double m43, double m44);
		///@}

		/**@name Getter
		*	Getter for class FAMath::Matrix4x4.
		*/
		///@{

		/**@brief Returns the elements of column index as a Vector4.
		* 
		*	Throws an std::out_of_range if given index > 3.
		*/
		Vector4 column(const unsigned int& index) const;

		/**@brief Returns the elements of row index as a Vector4.
		* 
		*	Throws an std::out_of_range if given index > 3
		*/
		Vector4 row(const unsigned int& index) const;

		/**@brief Returns a pointer to the raw data of this matrix.
		*
		*	The raw data is stored in column-major format.
		*/
		double* data();

		/**@brief Returns a constant pointer to the raw data of this matrix.
		* 
		*	The raw data is stored in column-major format.
		*/
		const double* data() const;

		/**@brief Returns a constant pointer to the raw data of this matrix.
		*	
		*	The raw data is stored in column-major format.
		*/
		const double* constData() const;
		///@}

		/**@name Setter
		*	Setter for class FAMath::Matrix4x4.
		*/
		///@{
		
		/**@brief Sets the element at [row][col] to the specified value.
		*	
		*	Throws an std::out_of_range if row or col > 3.
		*/
		void set(const unsigned int& row, const unsigned int& col, const double& value);
		
		/**@brief Sets the matrix to the identity matrix.
		*/
		void setToIdentity();


		/**@brief Sets all of the elements of the matrix to value.
		*/
		void fill(const double& value);

		/**@brief Sets the elements of the column index to the components of the Vector4 object value.
		* 
		*	Throws a std::out_of_range exception if index > 3.
		*/
		void setColumn(const unsigned int& index, const Vector4& value);

		/**@brief Sets the elements of the row index to the components of the Vector4 object value.
		*
		*	Throws a std::out_of_range exception if index > 3.
		*/
		void setRow(const unsigned int& index, const Vector4& value);

		///@}

		/**@brief Returns true if the matrix is the identity matrix, false otherwise.
		*/
		bool isIdentity() const;

		/**@brief Returns this matrix, transposed about its diagonal.
		*/
		Matrix4x4 transposed() const;

		/**@name Operator Overloading Member Functions
		*	Operator Overloading Member Functions for class FAMath::Matrix4x4.
		*/
		///@{

		/**@brief Returns a constant reference to the element at position [row][col] in this matrix.
		*
		*	Throws an std::out_of_range if row or col > 3.
		*/
		const double& operator()(const unsigned int& row, const unsigned int& col) const;

		/**@brief Returns a reference to the element at position [row][col] in this matrix,
		*
		*	Throws an std::out_of_range if row or col > 3.
		*/
		double& operator()(const unsigned int& row, const unsigned int& col);

		/**@brief Add a 4x4 matrix with another 4x4 matrix through overloading operator +=.
		*
		*	@return a reference to the current Matrix4x4 object with the result of the current Matrix4x4 object + m.
		*/
		Matrix4x4& operator+=(const Matrix4x4& m);

		/**@brief Subrtact a 4x4 matrix with another 4x4 matrix through overloading operator -=.
		*
		*	@return a reference to the current Matrix4x4 object with the result of the current Matrix4x4 object - m.
		*/
		Matrix4x4& operator-=(const Matrix4x4& m);

		/**@brief Multiplying a 4x4 matrix with a scalar through overloading operator *=.
		*
		*	@return a reference to the current Matrix4x4 object with the result of the current Matrtix4x4 object * scalar.
		*/
		Matrix4x4& operator*=(const double& scalar);

		/**@brief Multiplies two 4x4 matrices through overloading operator *=.
		*
		*	@return a reference to the current Matrix4x4 object with the result of the current Matrtix4x4 object * scalar.
		*/
		Matrix4x4& operator*=(const Matrix4x4& m);

		/**@brief Divides a 4x4 matrix with a scalar through overloading operator *=.
		*
		*	@return a reference to the current Matrix4x4 object with the result of the current Matrtix4x4 object / scalar.
		* 
		*	Throws an invalid_argument exception if scalar is 0.0.
		*/
		Matrix4x4& operator/=(const double& scalar);

		/**@brief Multiplies this matrix by another that rotates angle degrees about vector v.
		*/
		void rotate(double angle, const Vector3& v);

		/**@brief Multiplies this matrix by another that rotates angle degrees about vector(x, y, z).
		*/
		void rotate(const double& angle, const double& x, const double& y, const double& z);

		/**@brief Multiplies this matrix by another that rotates the coordinates using the quaternion matrix.
		* 
		* (x, y, z) in v is the axis you want to rotate around normalized. The w value is the angle in degrees.
		*/
		void rotateUsingQuaternion(const Vector4& v);

		/**@brief Multiplies this matrix by another that rotates the coordinates using the quaternion matrix.
		*
		*	v is the axis to rotate around normalized.
		* 
		*	The angle should be given in degrees.
		*/
		void rotateUsingQuaternion(const double& angle, const Vector3& v);

		/**@brief Multiplies this matrix by another that rotates the coordinates using the quaternion matrix.
		*
		*	(x, y, z) is the axis to rotate around normalized.
		* 
		*	The angle should be given in degrees.
		*/
		void rotateUsingQuaternion(const double& angle, const double& x, const double& y, const double& z);

		/**@brief Multiplies this matrix by another that scales the coordinates by the components of vector v.
		*/
		void scale(const Vector3& v);

		/**@brief Multiplies this matrix by another that scales the coordi nates by the components x and y.
		*/
		void scale(const double& x, const double& y);

		/**@brief Multiplies this matrix by another that scales the coordinates by the components x, y and z.
		*
		*/
		void scale(const double& x, const double& y, const double& z);

		/**@brief Multiplies this matrix by another that scales the coordinates by the specified factor.
		*/
		void scale(const double& factor);

		/**@brief Multiplies this matrix by another that scales the coordinates by the specified factor along vector v.
		*/
		void scale(const Vector3& v, const double& factor);

		/**@brief Multiplies this matrix by another that translates coordinates by the components of v.
		*/
		void translate(const Vector3& v);

		/**@brief Multiplies this matrix by another that translates coordinates by the components of x and y.
		*/
		void translate(const double& x, const double& y);

		/**@brief Multiplies this matrix by another that translates coordinates by the components of x, y and z.
		*/
		void translate(const double& x, const double& y, const double& z);

		/**@brief Multiplies this matrix by another that applies an othrographic projection for a window with lower left corner
		* (left, bottom), upper right corner (right, top) and the specified near and far clipping planes.
		*/
		void ortho(const double& left, const double& right, const double& bottom, const double& top, const double& near, const double& far);

		/**@brief Multiplies this matrix by another that applies a persepective projection. 
		*
		*	The fov is the vertical angle in degrees.
		*	Aspect ratio is the aspect ratio of your window.
		*	Near and far are the distances from the viewer to the corresponding planes.
		*/
		void perspective(const double& fov, const double& aspectRatio, const double& near, const double& far);

		/**@brief Returns the determinant of this matrix.
		*/
		double determinant() const;

		///@}


		/**@brief Adds two 4x4 matrices overloading operator +.
		*
		*	@return a Matrix4x4 object with the result of m1 + m2.
		* 
		*/
		friend Matrix4x4 operator+(const Matrix4x4& m1, const Matrix4x4& m2);

		/**@brief Subtracts two 4x4 matrices overloading operator +.
		*
		*	@return a Matrix4x4 object with the result of m1 - m2.
		*
		*/
		friend Matrix4x4 operator-(const Matrix4x4& m1, const Matrix4x4& m2);

		/**@brief Negates the 4x4 matrix through overloading operator -.
		*
		*	@return a Matrix4x4 object with the result of -m.
		*
		*/
		friend Matrix4x4 operator-(Matrix4x4& m);
	
		/**@brief Multiplying a 4x4 matrix with a scalar through overloading operator *.
		*	
		*	@return a Matrix4x4 object with the result of m1 * scalar.
		*/
		friend Matrix4x4 operator*(const Matrix4x4& m1, const double& scalar);

		/**@brief Multiplying a 4x4 matrix with a scalar through overloading operator *.
		* 
		*	@return a Matrix4x4 object with the result of scalar * m1.
		*/
		friend Matrix4x4 operator*(const double& scalar, const Matrix4x4& m1);

		/**@brief Multiplies two 4x4 matrices through overloading operator *.
		*
		*	@return a Matrix4x4 object with the result of m1 * m2.
		*/
		friend Matrix4x4 operator*(const Matrix4x4& m1, const Matrix4x4& m2);

		/**@brief Multiplies a 4x4 matrix with a column vector(4x1) through overloading operator *.
		*
		*	@return a Vector4 object with the result of m * vec.
		*/
		friend Vector4 operator*(const Matrix4x4& m, const Vector4& vec);

		/**@brief Multiplies a 4x4 matrix with a row vector(1x4) through overloading operator *.
		*
		*	@return a Vector4 object with the result of vec * m.
		*/
		friend Vector4 operator*(const Vector4& vec, const Matrix4x4& m);


		/**@brief Divides a 4x4 matrix with a scalar through overloading operator.
		*	
		*	@return a Matrix4x4 object with the result of m1 / scalar.
		* 
		*	Throws an invalid_argument exception if scalar is 0.0.
		*/
		friend Matrix4x4 operator/(const Matrix4x4& m1, const double& scalar);

		/**@brief Returns true if m1 is identical to m2, false otherwise.
		*/
		friend bool operator==(const Matrix4x4& m1, const Matrix4x4& m2);

		/**@brief Returns false if m1 is identical to m2, true otherwise.
		*/
		friend bool operator!=(const Matrix4x4& m1, const Matrix4x4& m2);

		/**@brief Returns the inverse of the matrix m.
		* If the matrix can't be inverted then the identity matrix is returned.
		*/
		friend Matrix4x4 inverse(const Matrix4x4& m);

		friend void print(const Matrix4x4& m);

	private:
		//A static array of 16 doubles
		//1st row for m_matrix is at indices: 0 4 8 12
		//2nd row for m_matrix is at indices: 1 5 9 13
		//3rd row for m_matrix is at indices: 2 6 10 14
		//4th row for m_matrix is at indices: 3 7 11 15
		//1st column for m_matrix is at indices: 0 1 2 3
		//2nd column for m_matrix is at indices: 4 5 6 7
		//3rd column for m_matrix is at indices: 8 9 10 11
		//4th column for m_matrix is at indices: 12 13 14 15
		double m_matrix[16];
		
		/**@brief Returns the minor of a given row and column of this matrix
		*/
		double minor(const unsigned int& row, const unsigned int& col) const;

		/**@brief Returns the adjoint of this matrix
		*/
		Matrix4x4 adjoint() const;
	};





	/** @class Quaternion ""
	*	@brief A quaternion class used to represent rotations.
	*
	*	A quaternion consists of a scalar to represent rotation angle and a 3D vector to represent an axis.
	*	The classes uses a double to respresent the scalar and a 3D vector to represent the axis.
	*/
	class Quaternion
	{
	public:

		/** @name Constructors
		* Constructors for class FAMath::Quaternion
		*/
		///@{

		/**@brief Default Constructor.
		* 
		* Creates a new Quaternion with scalar value = 1 and the 3D vector = (0, 0, 0).
		*/
		Quaternion();

		/**@brief Overloaded Constructor.
		*
		* Creates a new Quaternion with scalar value equal to w and
		* the 3D vector equal to v.
		*/
		Quaternion(const double& w, const Vector3& v);

		/**@brief Overloaded Constructor.
		*
		* Creates a new Quaternion with scalar value equal to w and
		* the 3D vector equal to (x, y ,z).
		*/
		Quaternion(const double& w, const double& x, const double& y, const double& z);

		///@}

		/** @name Getters and Setters
		* Getters and Setters for class FAMath::Quaternion
		*/
		///@{

		/**@brief Returns the scalar component of the quaternion.
		*/
		double scalar() const;

		/**@brief Returns the 3D vector component of the quaternion.
		*/
		Vector3 vector() const;

		/**@brief Returns the x component of the quaternion's 3D vector.
		*/
		double x() const;

		/**@brief Returns the y component of the quaternion's 3D vector.
		*/
		double y() const;

		/**@brief Returns the z component of the quaternion's 3D vector.
		*/
		double z() const;

		/**@brief Sets the quaternion values.
		* The scalar value equals to w and the 3D vector equal to v.
		*/
		void setQuaternion(const double& w, const Vector3& v);

		/**@brief Sets the quaternion values.
		* The scalar value equals to w and the 3D vector equals to (x, y, z).
		*/
		void setQuaternion(const double& w, double& x, double& y, double& z);

		/**@brief Sets the scalar value in the quaternion to w.
		*/
		void setScalar(const double& w);

		/**@brief Sets the 3D vector in the quaternion to v.
		*/
		void setVector(const Vector3& v);

		/**@brief Sets the 3D vector in the quaternion to (x, y, z).
		*/
		void setVector(const double& x, const double& y, const double& z);

		/**@brief Sets the x component of the 3D vector in the quaternion to the given x.
		*/
		void setX(const double& x);

		/**@brief Sets the y component of the 3D vector in the quaternion to the given y.
		*/
		void setY(const double& y);

		/**@brief Sets the z component of the 3D vector in the quaternion to the given z.
		*/
		void setZ(const double& z);
		///@}

		/**@brief Returns true if the quaternion scalar value is equal to 0 and all the components of the 3D vector is equal to zero, false otherwise.
		*/
		bool isZeroQuaternion() const;

		/**@brief Multiplies this quaternion by a scalar and returns a reference to this quaternion.
		*/
		Quaternion& operator*=(const double& k);

		/**@brief Multiplies this quaternion by q and returns a reference to this quaternion.
		*/
		Quaternion& operator*=(const Quaternion& q);

		/**@brief Adds this quaternion to the given quaternion q and returns a reference to this quaternion.
		*/
		Quaternion& operator+=(const Quaternion& q);

		/**@brief Subtracts this quaternion from the given quaternion q and returns a reference to this quaternion.
		*/
		Quaternion& operator-=(const Quaternion& q);

		/**@brief Creates a rotation matrix from this quaternion.
		* Normalize the quaternion before using this function.
		*/
		Matrix4x4 toRotationMatrix();

		/**@brief Negates the scalar value and each componenet in the 3D vector of q.
		* 
		* @return a Quaternion object that has the result of -q.
		*/
		friend Quaternion operator-(const Quaternion& q);

		/**@brief Returns the magnitude of the quaternion.
		*/
		friend double length(const Quaternion& q);

		/**@brief Normalizes the quaternion q.
		*
		* If the quaternion is a zero quaternion then the zero quaternion is returned.
		* 
		* @return a Quaternion object that has the result q / |q|.
		*/
		friend Quaternion normalize(const Quaternion& q);

		/**@brief Returns the conjugate of quaternion q.
		*/
		friend Quaternion conjugate(const Quaternion& q);

		/**@brief Returns the inverse of quaternion q.
		* If the quaternion q is the zero quaternion, then the zero quaternion is returned.
		*/
		friend Quaternion inverse(const Quaternion& q);

		/**@brief Returns the product of q1 and q2 using quaternion multiplication.
		*/
		friend Quaternion operator*(const Quaternion& q1, const Quaternion& q2);

		/**@brief Returns a Quaternion object that has the result of q * k.
		*/
		friend Quaternion operator*(const Quaternion& q, const double& k);

		/**@brief Returns a Quaternion object that has the result of k * q.
		*/
		friend Quaternion operator*(const double& k, const Quaternion& q);

		/**@brief Rotates 3D vector v by quaternion q to produce a new 3D vector.
		*/
		friend Vector3 operator*(const Quaternion& q, const Vector3& v);

		/**@brief Returns a Quaternion object that is the sum of q1 and q2.
		*/
		friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2);

		/**@brief Returns a Quaternion object that has the result of q1 - q2;
		*/
		friend Quaternion operator-(const Quaternion& q1, const Quaternion& q2);

		/**@brief Return true if q1 and q2 are equal, false otherwise.
		*/
		friend bool operator==(const Quaternion& q1, const Quaternion& q2);

		/**@brief Return true if q1 and q2 aren't equal, false otherwise.
		*/
		friend bool operator!=(const Quaternion& q1, const Quaternion& q2);

		/**@brief Returns the dot product between q1 and q2.
		*/
		friend double dotProduct(const Quaternion& q1, const Quaternion& q2);

		/**@brief Spherically Interpolates between rotations q1 and q2.
		* 
		* t should be between 0 and 1. If t < 0 q1 will be returned. If t > 1 q2 will be returned.
		*/
		friend Quaternion slerp(const Quaternion& q1, const Quaternion& q2, const double& t);


		friend void print(const Quaternion& q);


	private:
		//Scalar value of the quaternion
		double m_w;
		
		//3D vector of the quaternion
		Vector3 m_v;

	};
}

