#include <iostream>
#include <assert.h>
#include <limits>
#include <iomanip>
#include "MathLibraryTest.h"
#include "FAUtil.h"


#define MIN -10000.0
#define MAX 10000.0
#define LOOPNUM 1000000
#define EPSILON 1e-7
#define PI 3.14159265

//VECTOR2 TESTING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void Vector2DefaultConstructorTest()
{
	std::cout << "Vector2 Default Constructor Test Start\n" << std::endl;

	FAMath::Vector2* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector2();
		assert(FAUtil::compareDoubles(a->x(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), 0.0, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Default Constructor Test Finished\n" << std::endl;
}

void Vector2OverloadedConstructor1Test()
{
	std::cout << "Vector2 Overloaded Constructor 1 Test Start\n" << std::endl;

	FAMath::Vector2* a{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Vector2(x, y);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));

		delete a;
	}

	std::cout << "Vector2 Overloaded Constructor 1 Test Finished\n" << std::endl;
}

void Vector2OverloadedConstructor2Test()
{
	std::cout << "Vector2 Overloaded Constructor 2 Test Start\n" << std::endl;

	FAMath::Vector2* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector3(x, y, z);
		a = new FAMath::Vector2(*b);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector2 Overloaded Constructor 2 Test Finished\n" << std::endl;
}

void Vector2OverloadedConstructor3Test()
{
	std::cout << "Vector2 Overloaded Constructor 3 Test Start\n" << std::endl;

	FAMath::Vector2* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector4(x, y, z, w);
		a = new FAMath::Vector2(*b);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector2 Overloaded Constructor 3 Test Finished\n" << std::endl;
}

void Vector2AssignmentOperator1Test()
{
	std::cout << "Vector2 Assignment Operator 1 Test Start\n" << std::endl;

	FAMath::Vector2* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector3(x, y, z);
		a = new FAMath::Vector2();
		*a = *b;
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector2 Assignment Operator 1 Test Finished\n" << std::endl;
}

void Vector2AssignmentOperator2Test()
{
	std::cout << "Vector2 Assignment Operator 2 Test Start\n" << std::endl;

	FAMath::Vector2* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector4(x, y, z, w);
		a = new FAMath::Vector2();
		*a = *b;
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector2 Assignment Operator 2 Test Finished\n" << std::endl;
}

void Vector2SettersTest()
{
	std::cout << "Vector2 Setters Test Start\n" << std::endl;

	double x = 0.0;
	double y = 0.0;
	FAMath::Vector2* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2();
		a->setX(x);
		a->setY(y);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Setters Test Finished\n" << std::endl;
}

void Vector2AddtionTest()
{
	std::cout << "Vector2 Addition Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a = nullptr;
	FAMath::Vector2* b = nullptr;
	FAMath::Vector2 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);
		c = *a + *b;
		assert(FAUtil::compareDoubles(c.x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 + y2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Addition Test Finished\n" << std::endl;
}

void Vector2Addtion2Test()
{
	std::cout << "Vector2 Addition 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a = nullptr;
	FAMath::Vector2* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);
		*a += *b;
		assert(FAUtil::compareDoubles(a->x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 + y2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Addition 2 Test Finished\n" << std::endl;
}

void Vector2SubtractionTest()
{
	std::cout << "Vector2 Subtraction Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a = nullptr;
	FAMath::Vector2* b = nullptr;
	FAMath::Vector2 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);
		c = *a - *b;
		assert(FAUtil::compareDoubles(c.x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 - y2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Subtraction Test Finished\n" << std::endl;
}

void Vector2Subtraction2Test()
{
	std::cout << "Vector2 Subtraction 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a = nullptr;
	FAMath::Vector2* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);
		*a -= *b;
		assert(FAUtil::compareDoubles(a->x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 - y2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Subtraction 2 Test Finished\n" << std::endl;
}

void Vector2MultiplicationByAScalarTest()
{
	std::cout << "Vector2 Multiplication By A Scalar Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = *a * randomF;
		assert(FAUtil::compareDoubles(b.x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Multiplication By A Scalar Test Finished\n" << std::endl;
}

void Vector2MultiplicationByAScalar2Test()
{
	std::cout << "Vector2 Multiplication By A Scalar 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = randomF * *a;
		assert(FAUtil::compareDoubles(b.x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Multiplication By A Scalar 2 Test Finished\n" << std::endl;
}

void Vector2MultiplicationByAScalar3Test()
{
	std::cout << "Vector2 Multiplication By A Scalar 3 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector2* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		*a *= randomF;
		assert(FAUtil::compareDoubles(a->x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Multiplication By A Scalar 3 Test Finished\n" << std::endl;
}

void Vector2DivisionByAScalarTest()
{
	std::cout << "Vector2 Division By A Scalar Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		while (FAUtil::compareDoubles(randomF, 0.0, EPSILON))
		{
			randomF = FAUtil::randomRangeD(MIN, MAX);
		}
		a = new FAMath::Vector2(x1, y1);
		b = *a / randomF;
		assert(FAUtil::compareDoubles(b.x(), x1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 / randomF, EPSILON));
	}

	std::cout << "Vector2 Division By A Scalar Test Finished\n" << std::endl;
}

void Vector2DivisionByAScalar2Test()
{
	std::cout << "Vector2 Division By A Scalar 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector2* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		while (FAUtil::compareDoubles(randomF, 0.0, EPSILON))
		{
			randomF = FAUtil::randomRangeD(MIN, MAX);
		}
		a = new FAMath::Vector2(x1, y1);
		*a /= randomF;
		assert(FAUtil::compareDoubles(a->x(), x1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 / randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Division By A Scalar 2 Test Finished\n" << std::endl;
}

void Vector2EqualTest()
{
	std::cout << "Vector2 Equal Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a = nullptr;
	FAMath::Vector2* b = nullptr;
	FAMath::Vector2* c = nullptr;
	FAMath::Vector2* d = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x1, y1);

		c = new FAMath::Vector2(x1, y1);
		d = new FAMath::Vector2(x2, y2);

		assert(*a == *b);
		assert(!(*c == *d));

		delete a;
		delete b;
		delete c;
		delete d;
	}

	std::cout << "Vector2 Equal Test Finished\n" << std::endl;
}

void Vector2ZeroVectorTest()
{
	std::cout << "Vector2 Zero Vector Test Start\n" << std::endl;

	FAMath::Vector2* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector2();
		assert(a->isZeroVector());
		a->setX(FAUtil::randomRangeD(MIN, MAX) + 1.0);
		assert(!a->isZeroVector());
		delete a;
	}

	std::cout << "Vector2 Zero Vector Test Finished\n" << std::endl;
}

void Vector2LengthTest()
{
	std::cout << "Vector2 Length Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	FAMath::Vector2* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		assert(FAUtil::compareDoubles(length(*a), sqrt(x1 * x1 + y1 * y1), EPSILON));
		delete a;
	}

	std::cout << "Vector2 Length Test Finished\n" << std::endl;
}

void Vector2NormalizeTest()
{
	std::cout << "Vector2 Normalize Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector2();
		b = normalize(*a);
		assert(b.isZeroVector());
		delete a;

		x1 = FAUtil::randomRangeD(MIN, MAX) + 1.0;
		y1 = FAUtil::randomRangeD(MIN, MAX) + 1.0;
		a = new FAMath::Vector2(x1, y1);
		b = normalize(*a);
		assert(FAUtil::compareDoubles(length(b), 1.0, EPSILON));
		delete a;
	}

	std::cout << "Vector2 Normalize Test Finished\n" << std::endl;
}

void Vector2DistanceTest()
{
	std::cout << "Vector2 Distance Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);
		assert(FAUtil::compareDoubles(distance(*a, *b), sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)), EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Distance Test Finished\n" << std::endl;
}

void Vector2DotProductTest()
{
	std::cout << "Vector2 Dot Product Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2= FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);
		assert(FAUtil::compareDoubles(dotProduct(*a, *b), a->x() * b->x() + a->y() * b->y(), EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Dot Product Test Finished\n" << std::endl;
}

void Vector2AngleTest()
{
	std::cout << "Vector2 Angle Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double lenA{ 0.0 };
	double lenB{ 0.0 };
	double dP{ 0.0 };
	double ang{ 0.0 };
	FAMath::Vector2* a{ nullptr };
	FAMath::Vector2* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector2(x1, y1);
		b = new FAMath::Vector2(x2, y2);

		lenA = sqrt(x1 * x1 + y1 * y1);
		lenB = sqrt(x2 * x2 + y2 * y2);
		x1 = x1 / lenA;
		y1 = y1 / lenA;
		x2 = x2 / lenB;
		y2 = y2 / lenB;
		dP = x1 * x2 + y1 * y2;
		if (dP < -1.0)
		{
			dP = -1.0;
		}
		else if (dP > 1.0)
		{
			dP = 1.0;
		}
		ang = acos(dP) * 180 / PI;
		assert(FAUtil::compareDoubles(angle(normalize(*a), normalize(*b)), ang, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector2 Angle Test Finished\n" << std::endl;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------









//VECTOR3 TESTING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void Vector3DefaultConstructorTest()
{
	std::cout << "Vector3 Default Constructor Test Start\n" << std::endl;

	FAMath::Vector3* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector3();
		assert(FAUtil::compareDoubles(a->x(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Default Constructor Test Finished\n" << std::endl;
}

void Vector3OverloadedConstructor1Test()
{
	std::cout << "Vector3 Overloaded Constructor 1 Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Vector3(x, y, z);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));

		delete a;

	}

	std::cout << "Vector3 Overloaded Constructor 1 Test Finished\n" << std::endl;
}

void Vector3OverloadedConstructor2Test()
{
	std::cout << "Vector3 Overloaded Constructor 2 Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector3(*b);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector3 Overloaded Constructor 2 Test Finished\n" << std::endl;
}

void Vector3OverloadedConstructor3Test()
{
	std::cout << "Vector3 Overloaded Constructor 3 Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector3(*b, z);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector3 Overloaded Constructor 3 Test Finished\n" << std::endl;
}

void Vector3OverloadedConstructor4Test()
{
	std::cout << "Vector3 Overloaded Constructor 4 Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector4(x, y, z, w);
		a = new FAMath::Vector3(*b);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector3 Overloaded Constructor 4 Test Finished\n" << std::endl;
}

void Vector3AssignmentOperator1Test()
{
	std::cout << "Vector3 Assignment Operator 1 Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector3();
		*a = *b;
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector3 Assignment Operator 1 Test Finished\n" << std::endl;
}

void Vector3AssignmentOperator2Test()
{
	std::cout << "Vector3 Assignment Operator 2 Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector4(x, y, z, w);
		a = new FAMath::Vector3();
		*a = *b;
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector3 Assignment Operator 2 Test Finished\n" << std::endl;
}

void Vector3SettersTest()
{
	std::cout << "Vector3 Setters Test Start\n" << std::endl;

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	FAMath::Vector3* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3();
		a->setX(x);
		a->setY(y);
		a->setZ(z);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Setters Test Finished\n" << std::endl;
}

void Vector3AddtionTest()
{
	std::cout << "Vector3 Addition Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a = nullptr;
	FAMath::Vector3* b = nullptr;
	FAMath::Vector3 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		c = *a + *b;
		assert(FAUtil::compareDoubles(c.x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 + y2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z1 + z2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Addition Test Finished\n" << std::endl;
}

void Vector3Addtion2Test()
{
	std::cout << "Vector3 Addition 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a = nullptr;
	FAMath::Vector3* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		*a += *b;
		assert(FAUtil::compareDoubles(a->x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 + y2, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 + z2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Addition 2 Test Finished\n" << std::endl;
}

void Vector3SubtractionTest()
{
	std::cout << "Vector3 Subtraction Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a = nullptr;
	FAMath::Vector3* b = nullptr;
	FAMath::Vector3 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		c = *a - *b;
		assert(FAUtil::compareDoubles(c.x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 - y2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z1 - z2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Subtraction Test Finished\n" << std::endl;
}

void Vector3Subtraction2Test()
{
	std::cout << "Vector3 Subtraction 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a = nullptr;
	FAMath::Vector3* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		*a -= *b;
		assert(FAUtil::compareDoubles(a->x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 - y2, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 - z2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Subtraction 2 Test Finished\n" << std::endl;
}

void Vector3MultiplicationByAScalarTest()
{
	std::cout << "Vector3 Multiplication By A Scalar Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = *a * randomF;
		assert(FAUtil::compareDoubles(b.x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), z1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Multiplication By A Scalar Test Finished\n" << std::endl;
}

void Vector3MultiplicationByAScalar2Test()
{
	std::cout << "Vector3 Multiplication By A Scalar 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = randomF * *a;
		assert(FAUtil::compareDoubles(b.x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), z1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Multiplication By A Scalar 2 Test Finished\n" << std::endl;
}

void Vector3MultiplicationByAScalar3Test()
{
	std::cout << "Vector3 Multiplication By A Scalar 3 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector3* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		*a *= randomF;
		assert(FAUtil::compareDoubles(a->x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Multiplication By A Scalar 3 Test Finished\n" << std::endl;
}

void Vector3DivisionByAScalarTest()
{
	std::cout << "Vector3 Division By A Scalar Test Start\n" << std::endl;

	std::cout << "Vector3 Multiplication By A Scalar Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = *a / randomF;
		assert(FAUtil::compareDoubles(b.x(), x1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), z1 / randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Division By A Scalar Test Finished\n" << std::endl;
}

void Vector3DivisionByAScalar2Test()
{
	std::cout << "Vector3 Division By A Scalar 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector3* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		while (FAUtil::compareDoubles(randomF, 0.0, EPSILON))
		{
			randomF = FAUtil::randomRangeD(MIN, MAX);
		}
		a = new FAMath::Vector3(x1, y1, z1);
		*a /= randomF;
		assert(FAUtil::compareDoubles(a->x(), x1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 / randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Division By A Scalar 2 Test Finished\n" << std::endl;
}

void Vector3EqualTest()
{
	std::cout << "Vector3 Equal Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a = nullptr;
	FAMath::Vector3* b = nullptr;
	FAMath::Vector3* c = nullptr;
	FAMath::Vector3* d = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x1, y1, z1);

		c = new FAMath::Vector3(x1, y1, z1);
		d = new FAMath::Vector3(x2, y2, z2);

		assert(*a == *b);
		assert(!(*c == *d));

		delete a;
		delete b;
		delete c;
		delete d;
	}

	std::cout << "Vector3 Equal Test Finished\n" << std::endl;
}

void Vector3ZeroVectorTest()
{
	std::cout << "Vector3 Zero Vector Test Start\n" << std::endl;

	FAMath::Vector3* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector3();
		assert(a->isZeroVector());
		a->setX(FAUtil::randomRangeD(MIN, MAX) + 1.0);
		assert(!a->isZeroVector());
		delete a;
	}

	std::cout << "Vector3 Zero Vector Test Finished\n" << std::endl;
}

void Vector3LengthTest()
{
	std::cout << "Vector3 Length Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	FAMath::Vector3* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		assert(FAUtil::compareDoubles(length(*a), sqrt(x1 * x1 + y1 * y1 + z1 * z1), EPSILON));
		delete a;
	}

	std::cout << "Vector3 Length Test Finished\n" << std::endl;
}

void Vector3NormalizeTest()
{
	std::cout << "Vector3 Normalize Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector3();
		b = normalize(*a);
		assert(b.isZeroVector());
		delete a;

		x1 = FAUtil::randomRangeD(MIN, MAX) + 1.0;
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = normalize(*a);
		assert(FAUtil::compareDoubles(length(b), 1.0, EPSILON));
		delete a;
	}

	std::cout << "Vector3 Normalize Test Finished\n" << std::endl;
}

void Vector3DistanceTest()
{
	std::cout << "Vector3 Distance Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		assert(FAUtil::compareDoubles(distance(*a, *b), sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) *(z1 - z2)), EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Distance Test Finished\n" << std::endl;
}

void Vector3DotProductTest()
{
	std::cout << "Vector3 Dot Product Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		assert(FAUtil::compareDoubles(dotProduct(*a, *b), a->x() * b->x() + a->y() * b->y() + a->z() * b->z(), EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Dot Product Test Finished\n" << std::endl;
}

void Vector3AngleTest()
{
	std::cout << "Vector3 Angle Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double lenA{ 0.0 };
	double lenB{ 0.0 };
	double dP{ 0.0 };
	double ang{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);

		lenA = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
		lenB = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
		x1 = x1 / lenA;
		y1 = y1 / lenA;
		z1 = z1 / lenA;
		x2 = x2 / lenB;
		y2 = y2 / lenB;
		z2 = z2 / lenB;
		dP = x1 * x2 + y1 * y2 + z1 * z2;
		if (dP < -1.0)
		{
			dP = -1.0;
		}
		else if (dP > 1.0)
		{
			dP = 1.0;
		}
		ang = acos(dP) * 180 / PI;
		assert(FAUtil::compareDoubles(angle(normalize(*a), normalize(*b)), ang, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Angle Test Finished\n" << std::endl;
}

void Vector3CrossProductTest()
{
	std::cout << "Vector3 Cross Product Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	FAMath::Vector3* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	FAMath::Vector3 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector3(x1, y1, z1);
		b = new FAMath::Vector3(x2, y2, z2);
		c = crossProduct(*a, *b);
		assert(FAUtil::compareDoubles(c.x(), y1 * z2 - z1 * y2 , EPSILON));
		assert(FAUtil::compareDoubles(c.y(), z1 * x2 - x1 * z2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), x1 * y2 - y1 * x2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector3 Cross Product Test Finished\n" << std::endl;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------







//Vector4 TESTING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void Vector4DefaultConstructorTest()
{
	std::cout << "Vector4 Default Constructor Test Start\n" << std::endl;

	FAMath::Vector4* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector4();
		assert(FAUtil::compareDoubles(a->x(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), 0.0, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Default Constructor Test Finished\n" << std::endl;
}

void Vector4OverloadedConstructor1Test()
{
	std::cout << "Vector4 Overloaded Constructor 1 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w= FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Vector4(x, y, z, w);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w, EPSILON));

		delete a;
	}

	std::cout << "Vector4 Overloaded Constructor 1 Test Finished\n" << std::endl;
}

void Vector4OverloadedConstructor2Test()
{
	std::cout << "Vector4 Overloaded Constructor 2 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector4(*b);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), 0.0, EPSILON));

		delete a;
	}

	std::cout << "Vector4 Overloaded Constructor 2 Test Finished\n" << std::endl;
}

void Vector4OverloadedConstructor3Test()
{
	std::cout << "Vector4 Overloaded Constructor 3 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector4(*b, z);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), 0.0, EPSILON));

		delete a;
	}

	std::cout << "Vector4 Overloaded Constructor 3 Test Finished\n" << std::endl;
}

void Vector4OverloadedConstructor4Test()
{
	std::cout << "Vector4 Overloaded Constructor 4 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector4(*b, z, w);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w, EPSILON));

		delete a;
	}

	std::cout << "Vector4 Overloaded Constructor 4 Test Finished\n" << std::endl;
}

void Vector4OverloadedConstructor5Test()
{
	std::cout << "Vector4 Overloaded Constructor 5 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector3(x, y, z);
		a = new FAMath::Vector4(*b);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), 0.0, EPSILON));

		delete a;
	}

	std::cout << "Vector4 Overloaded Constructor 5 Test Finished\n" << std::endl;
}

void Vector4OverloadedConstructor6Test()
{
	std::cout << "Vector4 Overloaded Constructor 6 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector3(x, y, z);
		a = new FAMath::Vector4(*b, w);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w, EPSILON));

		delete a;
	}

	std::cout << "Vector4 Overloaded Constructor 6 Test Finished\n" << std::endl;
}

void Vector4AssignmentOperator1Test()
{
	std::cout << "Vector4 Assignment Operator 1 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector2* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector2(x, y);
		a = new FAMath::Vector4();
		*a = *b;
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), 0.0, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector4 Assignment Operator 1 Test Finished\n" << std::endl;
}

void Vector4AssignmentOperator2Test()
{
	std::cout << "Vector4 Assignment Operator 2 Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new FAMath::Vector3(x, y, z);
		a = new FAMath::Vector4();
		*a = *b;
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), 0.0, EPSILON));

		delete a;
		delete b;
	}

	std::cout << "Vector4 Assignment Operator 2 Test Finished\n" << std::endl;
}

void Vector4SettersTest()
{
	std::cout << "Vector4 Setters Test Start\n" << std::endl;

	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double w = 0.0;
	FAMath::Vector4* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		w = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4();
		a->setX(x);
		a->setY(y);
		a->setZ(z);
		a->setW(w);
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Setters Test Finished\n" << std::endl;
}

void Vector4AddtionTest()
{
	std::cout << "Vector4 Addition Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a = nullptr;
	FAMath::Vector4* b = nullptr;
	FAMath::Vector4 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);
		c = *a + *b;
		assert(FAUtil::compareDoubles(c.x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 + y2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z1 + z2, EPSILON));
		assert(FAUtil::compareDoubles(c.w(), w1 + w2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Addition Test Finished\n" << std::endl;
}

void Vector4Addtion2Test()
{
	std::cout << "Vector4 Addition 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a = nullptr;
	FAMath::Vector4* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);
		*a += *b;
		assert(FAUtil::compareDoubles(a->x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 + y2, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 + z2, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w1 + w2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Addition 2 Test Finished\n" << std::endl;
}

void Vector4SubtractionTest()
{
	std::cout << "Vector4 Subtraction Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a = nullptr;
	FAMath::Vector4* b = nullptr;
	FAMath::Vector4 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);
		c = *a - *b;
		assert(FAUtil::compareDoubles(c.x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 - y2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z1 - z2, EPSILON));
		assert(FAUtil::compareDoubles(c.w(), w1 - w2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Subtraction Test Finished\n" << std::endl;
}

void Vector4Subtraction2Test()
{
	std::cout << "Vector4 Subtraction 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a = nullptr;
	FAMath::Vector4* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);
		*a -= *b;
		assert(FAUtil::compareDoubles(a->x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 - y2, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 - z2, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w1 - w2, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Subtraction 2 Test Finished\n" << std::endl;
}

void Vector4MultiplicationByAScalarTest()
{
	std::cout << "Vector4 Multiplication By A Scalar Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = *a * randomF;
		assert(FAUtil::compareDoubles(b.x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), z1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.w(), w1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Multiplication By A Scalar Test Finished\n" << std::endl;
}

void Vector4MultiplicationByAScalar2Test()
{
	std::cout << "Vector4 Multiplication By A Scalar 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = randomF * *a;
		assert(FAUtil::compareDoubles(b.x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), z1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.w(), w1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Multiplication By A Scalar 2 Test Finished\n" << std::endl;
}

void Vector4MultiplicationByAScalar3Test()
{
	std::cout << "Vector4 Multiplication By A Scalar 3 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector4* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		*a *= randomF;
		assert(FAUtil::compareDoubles(a->x(), x1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 * randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w1 * randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Multiplication By A Scalar 3 Test Finished\n" << std::endl;
}

void Vector4DivisionByAScalarTest()
{
	std::cout << "Vector4 Division By A Scalar Test Start\n" << std::endl;

	std::cout << "Vector4 Multiplication By A Scalar Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = *a / randomF;
		assert(FAUtil::compareDoubles(b.x(), x1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), y1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), z1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(b.w(), w1 / randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Division By A Scalar Test Finished\n" << std::endl;
}

void Vector4DivisionByAScalar2Test()
{
	std::cout << "Vector4 Division By A Scalar 2 Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double randomF{ 0.0 };
	FAMath::Vector4* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		randomF = FAUtil::randomRangeD(MIN, MAX);
		while (FAUtil::compareDoubles(randomF, 0.0, EPSILON))
		{
			randomF = FAUtil::randomRangeD(MIN, MAX);
		}
		a = new FAMath::Vector4(x1, y1, z1, w1);
		*a /= randomF;
		assert(FAUtil::compareDoubles(a->x(), x1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z1 / randomF, EPSILON));
		assert(FAUtil::compareDoubles(a->w(), w1 / randomF, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Division By A Scalar 2 Test Finished\n" << std::endl;
}

void Vector4EqualTest()
{
	std::cout << "Vector4 Equal Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a = nullptr;
	FAMath::Vector4* b = nullptr;
	FAMath::Vector4* c = nullptr;
	FAMath::Vector4* d = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX); 
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x1, y1, z1, w1);

		c = new FAMath::Vector4(x1, y1, z1, w1);
		d = new FAMath::Vector4(x2, y2, z2, w2);

		assert(*a == *b);
		assert(!(*c == *d));

		delete a;
		delete b;
		delete c;
		delete d;
	}

	std::cout << "Vector4 Equal Test Finished\n" << std::endl;
}

void Vector4ZeroVectorTest()
{
	std::cout << "Vector4 Zero Vector Test Start\n" << std::endl;

	FAMath::Vector4* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector4();
		assert(a->isZeroVector());
		a->setX(FAUtil::randomRangeD(MIN, MAX) + 1.0);
		assert(!a->isZeroVector());
		delete a;
	}

	std::cout << "Vector4 Zero Vector Test Finished\n" << std::endl;
}

void Vector4LengthTest()
{
	std::cout << "Vector4 Length Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	FAMath::Vector4* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		assert(FAUtil::compareDoubles(length(*a), sqrt(x1 * x1 + y1 * y1 + z1 * z1 + w1 * w1), EPSILON));
		delete a;
	}

	std::cout << "Vector4 Length Test Finished\n" << std::endl;
}

void Vector4NormalizeTest()
{
	std::cout << "Vector4 Normalize Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4 b;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Vector4();
		b = normalize(*a);
		assert(b.isZeroVector());
		delete a;

		x1 = FAUtil::randomRangeD(MIN, MAX) + 1.0;
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = normalize(*a);
		assert(FAUtil::compareDoubles(length(b), 1.0, EPSILON));
		delete a;
	}

	std::cout << "Vector4 Normalize Test Finished\n" << std::endl;
}

void Vector4DistanceTest()
{
	std::cout << "Vector4 Distance Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);
		assert(FAUtil::compareDoubles(distance(*a, *b), sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2) + (w1 - w2) * (w1 - w2)), EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Distance Test Finished\n" << std::endl;
}

void Vector4DotProductTest()
{
	std::cout << "Vector4 Dot Product Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);
		assert(FAUtil::compareDoubles(dotProduct(*a, *b), a->x() * b->x() + a->y() * b->y() + a->z() * b->z() + a->w() * b->w(), EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Dot Product Test Finished\n" << std::endl;
}

void Vector4AngleTest()
{
	std::cout << "Vector4 Angle Test Start\n" << std::endl;

	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w1{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };
	double w2{ 0.0 };
	double lenA{ 0.0 };
	double lenB{ 0.0 };
	double dP{ 0.0 };
	double ang{ 0.0 };
	FAMath::Vector4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Vector4(x1, y1, z1, w1);
		b = new FAMath::Vector4(x2, y2, z2, w2);

		lenA = sqrt(x1 * x1 + y1 * y1 + z1 * z1 + w1 * w1);
		lenB = sqrt(x2 * x2 + y2 * y2 + z2 * z2 + w2 * w2);
		x1 = x1 / lenA;
		y1 = y1 / lenA;
		z1 = z1 / lenA;
		w1 = w1 / lenA;
		x2 = x2 / lenB;
		y2 = y2 / lenB;
		z2 = z2 / lenB;
		w2 = w2 / lenB;
		dP = x1 * x2 + y1 * y2 + z1 * z2 + w1 * w2;
		if (dP < -1.0)
		{
			dP = -1.0;
		}
		else if (dP > 1.0)
		{
			dP = 1.0;
		}
		ang = acos(dP) * 180 / PI;
		assert(FAUtil::compareDoubles(angle(normalize(*a), normalize(*b)), ang, EPSILON));
		delete a;
		delete b;
	}

	std::cout << "Vector4 Angle Test Finished\n" << std::endl;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------






//MATRIX4x4 TESTING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------

void Matrix4x4DefaultConstructorTest()
{
	std::cout << "Matrix4x4 Default Constructor Test Start\n" << std::endl;

	FAMath::Matrix4x4* a;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();

		assert(a->isIdentity());

		delete a;
	}
	std::cout << "Matrix4x4 Default Constructor Test Finished\n" << std::endl;
}

void Matrix4x4OverloadedConstructor1Test()
{
	std::cout << "Matrix4x4 Overloaded Constructor 1 Test Start\n" << std::endl;

	FAMath::Matrix4x4* b{ nullptr };
	const double* c{ nullptr };
	double* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			a[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		b = new FAMath::Matrix4x4(a);
		c = b->constData();

		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(a[i * 4], c[i], EPSILON));
			assert(FAUtil::compareDoubles(a[i * 4 + 1], c[i + 4], EPSILON));
			assert(FAUtil::compareDoubles(a[i * 4 + 2], c[i + 8], EPSILON));
			assert(FAUtil::compareDoubles(a[i * 4 + 3], c[i + 12], EPSILON));
		}

		delete[] a;
		delete b;
	}

	std::cout << "Matrix4x4 Overloaded Constructor 1 Test Finished\n" << std::endl;
}

void Matrix4x4OverloadedConstructor2Test()
{
	std::cout << "Matrix4x4 Overloaded Constructor 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* b{ nullptr };
	const double* c{ nullptr };;
	double* a{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			a[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		b = new FAMath::Matrix4x4(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]);
		c = b->constData();

		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(a[i * 4], c[i], EPSILON));
			assert(FAUtil::compareDoubles(a[i * 4 + 1], c[i + 4], EPSILON));
			assert(FAUtil::compareDoubles(a[i * 4 + 2], c[i + 8], EPSILON));
			assert(FAUtil::compareDoubles(a[i * 4 + 3], c[i + 12], EPSILON));
		}

		delete[] a;
		delete b;
	}

	std::cout << "Matrix4x4 Overloaded Constructor 2 Test Finished\n" << std::endl;
}

void Matrix4x4ColumnTest()
{
	std::cout << "Matrix4x4 Column Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	FAMath::Vector4 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		for (unsigned int j = 0; j < 4; ++j)
		{
			c = a->column(j);
			assert(FAUtil::compareDoubles(c.x(), b[j], EPSILON));
			assert(FAUtil::compareDoubles(c.y(), b[j + 4], EPSILON));
			assert(FAUtil::compareDoubles(c.z(), b[j + 8], EPSILON));
			assert(FAUtil::compareDoubles(c.w(), b[j + 12], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Column Test Finished\n" << std::endl;
}

void Matrix4x4RowTest()
{
	std::cout << "Matrix4x4 Row Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	FAMath::Vector4 c;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		for (unsigned int j = 0; j < 4; ++j)
		{
			c = a->row(j);
			assert(FAUtil::compareDoubles(c.x(), b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c.y(), b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c.z(), b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c.w(), b[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Row Test Finished\n" << std::endl;
}

void Matrix4x4SetTest()
{
	std::cout << "Matrix4x4 Set Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4();

		for (unsigned int j = 0; j < 4; ++j)
		{	
			for (unsigned int k = 0; k < 4; ++k)
			{
				a->set(j, k, b[j * 4 + k]);
			}
		}

		c = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j], b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 4], b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 8], b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 12], b[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Set Test Finished \n" << std::endl;
}

void Matrix4x4SetToIdentityTest()
{
	std::cout << "Matrix4x4 Set To Identity Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	double c[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		a->setToIdentity();

		assert(a->isIdentity());

		d = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], c[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], c[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], c[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], c[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Set To Identity Test Finished\n" << std::endl;
}

void Matrix4x4SetColumnTest()
{
	std::cout << "Matrix4x4 Set Column Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	FAMath::Vector4 c;
	int rI{ 0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		b = new FAMath::Vector4(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		rI = FAUtil::randomRange(0, 3);

		a->setColumn(rI, *b);

		assert(a->column(rI) == *b);

		delete a;
		delete b;
	}

	std::cout << "Matrix4x4 Set Column Test Finished\n" << std::endl;
}

void Matrix4x4SetRowTest()
{
	std::cout << "Matrix4x4 Set Row Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	FAMath::Vector4 c;
	int rI{ 0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		b = new FAMath::Vector4(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		rI = FAUtil::randomRange(0, 3);

		a->setRow(rI, *b);

		c = a->row(rI);

		assert(a->row(rI) == *b);

		delete a;
		delete b;
	}

	std::cout << "Matrix4x4 Set Row Test Finished\n" << std::endl;
}

void Matrix4x4Fill()
{
	std::cout << "Matrix4x4 Fill Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	const double* b{ nullptr };
	double rF = 0;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		rF = FAUtil::randomRangeD(MIN, MAX);

		a->fill(rF);
		b = a->constData();

		for (unsigned int i = 0; i < 16; ++i)
		{
			assert(FAUtil::compareDoubles(b[i], rF, EPSILON));
		}
		

		delete a;
	}

	std::cout << "Matrix4x4 Fill Test Finished\n" << std::endl;
}

void Matrix4x4TransposedTest()
{
	std::cout << "Matrix4x4 Transposed Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		c = a->transposed().constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j * 4], b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 1], b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 2], b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 3], b[j * 4 + 3], EPSILON));
		}

		c = a->transposed().transposed().constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j], b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 4], b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 8], b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 12], b[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Transposed Test Finished\n" << std::endl;

}

void Matrix4x4IsIdentityTest()
{
	std::cout << "Matrix4x4 Is Identity Test Start\n" << std::endl;

	FAMath::Matrix4x4* a = nullptr;
	double* b = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}


		a = new FAMath::Matrix4x4(b);
		assert(!a->isIdentity());

		a->setToIdentity();
		assert(a->isIdentity());

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Is Identity Test Finished\n" << std::endl;
}

void Matrix4x4AdditionTest()
{
	std::cout << "Matrix4x4 Addition Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4* b{ nullptr };
	double* c{ nullptr };
	double* d{ nullptr };
	const double* e{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(c);
		b = new FAMath::Matrix4x4(d);

		*a += *b;

		e = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(e[j], c[j * 4] + d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 4], c[j * 4 + 1] + d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 8], c[j * 4 + 2] + d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 12], c[j * 4 + 3] + d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] c;
		delete[] d;
	}
	std::cout << "Matrix4x4 Addition Test Finished\n" << std::endl;
}

void Matrix4x4Addition2Test()
{
	std::cout << "Matrix4x4 Addition 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4* b{ nullptr };
	double* c{ nullptr };
	double* d{ nullptr };
	const double* e{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(c);
		b = new FAMath::Matrix4x4(d);

		e = (*a + *b).constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(e[j], c[j * 4] + d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 4], c[j * 4 + 1] + d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 8], c[j * 4 + 2] + d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 12], c[j * 4 + 3] + d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] c;
		delete[] d;
	}

	std::cout << "Matrix4x4 Addition 2 Test Finished\n" << std::endl;
}

void Matrix4x4SubtractionTest()
{
	std::cout << "Matrix4x4 Subtraction Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4* b{ nullptr };
	double* c{ nullptr };
	double* d{ nullptr };
	const double* e{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(c);
		b = new FAMath::Matrix4x4(d);

		*a -= *b;

		e = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(e[j], c[j * 4] - d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 4], c[j * 4 + 1] - d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 8], c[j * 4 + 2] - d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 12], c[j * 4 + 3] - d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] c;
		delete[] d;
	}
	std::cout << "Matrix4x4 Subtraction Test Finished\n" << std::endl;
}

void Matrix4x4Subtraction2Test()
{
	std::cout << "Matrix4x4 Subtraction 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4* b{ nullptr };
	double* c{ nullptr };
	double* d{ nullptr };
	const double* e{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(c);
		b = new FAMath::Matrix4x4(d);

		e = (*a - *b).constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(e[j], c[j * 4] - d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 4], c[j * 4 + 1] - d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 8], c[j * 4 + 2] - d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(e[j + 12], c[j * 4 + 3] - d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] c;
		delete[] d;
	}
	std::cout << "Matrix4x4 Subtraction 2 Test Finished\n" << std::endl;
}

void Matrix4x4NegationTest()
{
	std::cout << "Matrix4x4 Negation Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4 b;
	double* c{ nullptr };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(c);

		b = -(*a);

		d = b.constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], -c[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], -c[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], -c[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], -c[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] c;
	}

	std::cout << "Matrix4x4 Negation Test Finished\n" << std::endl;
}

void Matrix4x4MultiplicationByAScalarTest()
{
	std::cout << "Matrix4x4 Multiplication by A Scalar Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	double c{ 0.0 };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		c = FAUtil::randomRangeD(MIN, MAX);

		*a *= c;

		d = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], b[j * 4] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], b[j * 4 + 1] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], b[j * 4 + 2] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], b[j * 4 + 3]* c, EPSILON));
		}

		delete a;
		delete[] b;
	}
	std::cout << "Matrix4x4 Multiplication by A Scalar Test Finished\n" << std::endl;
}

void Matrix4x4MultiplicationByAScalar2Test()
{
	std::cout << "Matrix4x4 Multiplication by A Scalar 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	double c{ 0.0 };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		c = FAUtil::randomRangeD(MIN, MAX);

		d = (*a * c).constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], b[j * 4] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], b[j * 4 + 1] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], b[j * 4 + 2] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], b[j * 4 + 3] * c, EPSILON));
		}

		delete a;
		delete[] b;
	}
	std::cout << "Matrix4x4 Multiplication by A Scalar 2 Test Finished\n" << std::endl;
}

void Matrix4x4MultiplicationByAScalar3Test()
{
	std::cout << "Matrix4x4 Multiplication by A Scalar 3 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	double c{ 0.0 };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		c = FAUtil::randomRangeD(MIN, MAX);

		d = (c * *a).constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], b[j * 4] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], b[j * 4 + 1] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], b[j * 4 + 2] * c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], b[j * 4 + 3] * c, EPSILON));
		}

		delete a;
		delete[] b;
	}
	std::cout << "Matrix4x4 Multiplication by A Scalar 3 Test Finished\n" << std::endl;
}

void Matrix4x4DivisionByAScalarTest()
{
	std::cout << "Matrix4x4 Division by A Scalar Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	double c{ 0.0 };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		c = FAUtil::randomRangeD(MIN, MAX);

		*a /= c;

		d = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], b[j * 4] / c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], b[j * 4 + 1] / c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], b[j * 4 + 2] / c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], b[j * 4 + 3] / c, EPSILON));
		}

		delete a;
		delete[] b;
	}
	std::cout << "Matrix4x4 Division by A Scalar Test Finished\n" << std::endl;
}

void Matrix4x4DivisionByAScalar2Test()
{
	std::cout << "Matrix4x4 Division by A Scalar 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	double c{ 0.0 };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);

		c = FAUtil::randomRangeD(MIN, MAX);

		d = (*a / c).constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], b[j * 4] / c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], b[j * 4 + 1] / c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], b[j * 4 + 2] / c, EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], b[j * 4 + 3] / c, EPSILON));
		}

		delete a;
		delete[] b;
	}
	std::cout << "Matrix4x4 Division by A Scalar 2 Test Finished\n" << std::endl;
}

void Matrix4x4MatrixMultiplicationTest()
{
	std::cout << "Matrix4x4 Matrix Multiplication Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4* b{ nullptr };
	double* c{ nullptr };
	double* d{ nullptr };
	double e[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	const double* f{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		for (unsigned int j = 0; j < 4; ++j)
		{
			unsigned int first = j * 4;
			unsigned int second = j * 4 + 1;
			unsigned int third = j * 4 + 2;
			unsigned int fourth = j * 4 + 3;

			e[first] = (c[first] * d[0]) + (c[second] * d[4]) + (c[third] * d[8]) + (c[fourth] * d[12]);
			e[second] = (c[first] * d[1]) + (c[second] * d[5]) + (c[third] * d[9]) + (c[fourth] * d[13]);
			e[third] = (c[first] * d[2]) + (c[second] * d[6]) + (c[third] * d[10]) + (c[fourth] * d[14]);
			e[fourth] = (c[first] * d[3]) + (c[second] * d[7]) + (c[third] * d[11]) + (c[fourth] * d[15]);
		}


		a = new FAMath::Matrix4x4(c);
		b = new FAMath::Matrix4x4(d);

		*a *= *b;

		f = a->constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(f[j], e[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 4], e[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 8], e[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 12], e[j * 4 + 3], EPSILON));
		}

		for (unsigned int j = 0; j < 16; ++j)
		{
			e[j] = 0.0;
		}
			

		delete a;
		delete b;
		delete[] c;
		delete[] d;
	}

	std::cout << "Matrix4x4 Matrix Multiplication Test Finished\n" << std::endl;
}

void Matrix4x4MatrixMultiplication2Test()
{
	std::cout << "Matrix4x4 Matrix Multiplication 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Matrix4x4* b{ nullptr };
	double* c{ nullptr };
	double* d{ nullptr };
	double e[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	const double* f{ nullptr };
	const double* g{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		c = new double[16];
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			c[i] = FAUtil::randomRangeD(MIN, MAX);
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		for (unsigned int j = 0; j < 4; ++j)
		{
			unsigned int first = j * 4;
			unsigned int second = j * 4 + 1;
			unsigned int third = j * 4 + 2;
			unsigned int fourth = j * 4 + 3;

			e[first] = (c[first] * d[0]) + (c[second] * d[4]) + (c[third] * d[8]) + (c[fourth] * d[12]);
			e[second] = (c[first] * d[1]) + (c[second] * d[5]) + (c[third] * d[9]) + (c[fourth] * d[13]);
			e[third] = (c[first] * d[2]) + (c[second] * d[6]) + (c[third] * d[10]) + (c[fourth] * d[14]);
			e[fourth] = (c[first] * d[3]) + (c[second] * d[7]) + (c[third] * d[11]) + (c[fourth] * d[15]);
		}

		a = new FAMath::Matrix4x4(c);
		b = new FAMath::Matrix4x4(d);

		f = (*a * *b).constData();
		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(f[j], e[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 4], e[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 8], e[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 12], e[j * 4 + 3], EPSILON));
		}

		f = (*a * *b).transposed().constData();
		*a = a->transposed();
		*b = b->transposed();
		g = (*b * *a).constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(f[j], g[j], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 4], g[j + 4], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 8], g[j + 8], EPSILON));
			assert(FAUtil::compareDoubles(f[j + 12], g[j + 12], EPSILON));
		}

		for (unsigned int j = 0; j < 16; ++j)
		{
			e[j] = 0.0;
		}

		delete a;
		delete b;
		delete[] c;
		delete[] d;
	}

	std::cout << "Matrix4x4 Matrix Multiplication 2 Test Finished\n" << std::endl;
}

void Matrix4x4MatrixMultiplicationWithColumnVector()
{
	std::cout << "Matrix4x4 Multiplication with Column Vector Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	FAMath::Vector4 c;
	double* d{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(d);
		b = new FAMath::Vector4(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));

		x = d[0] * b->x() + d[1] * b->y() + d[2] * b->z() + d[3] * b->w();
		y = d[4] * b->x() + d[5] * b->y() + d[6] * b->z() + d[7] * b->w();
		z = d[8] * b->x() + d[9] * b->y() + d[10] * b->z() + d[11] * b->w();
		w = d[12] * b->x() + d[13] * b->y() + d[14] * b->z() + d[15] * b->w();

		c = *a * *b;

		assert(FAUtil::compareDoubles(c.x(), x, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z, EPSILON));
		assert(FAUtil::compareDoubles(c.w(), w, EPSILON));

		delete a;
		delete b;
		delete[] d;
	}

	std::cout << "Matrix4x4 Multiplication with Column Vector Test Finished\n" << std::endl;
}

void Matrix4x4MatrixMultiplicationWithRowVector()
{
	std::cout << "Matrix4x4 Multiplication with Row Vector Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	FAMath::Vector4 c;
	double* d{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		d = new double[16];

		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(d);
		b = new FAMath::Vector4(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));

		x = d[0] * b->x() + d[4] * b->y() + d[8] * b->z() + d[12] * b->w();
		y = d[1] * b->x() + d[5] * b->y() + d[9] * b->z() + d[13] * b->w();
		z = d[2] * b->x() + d[6] * b->y() + d[10] * b->z() + d[14] * b->w();
		w = d[3] * b->x() + d[7] * b->y() + d[11] * b->z() + d[15] * b->w();

		c = *b * *a;

		assert(FAUtil::compareDoubles(c.x(), x, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z, EPSILON));
		assert(FAUtil::compareDoubles(c.w(), w, EPSILON));

		delete a;
		delete b;
		delete[] d;
	}

	std::cout << "Matrix4x4 Multiplication with Row Vector Test Finished\n" << std::endl;
}


void Matrix4x4RotationTest()
{
	std::cout << "Matrix4x4 Rotation Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double* d{ nullptr };
	const double* e{ nullptr };
	double angle{ 0.0 };
	double rad{ 0.0 };
	double c{ 0.0 };
	double s{ 0.0 };
	double omc{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		b = new FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		d = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = 0.0;
		}

		angle = FAUtil::randomRangeD(0, 360);
		rad = angle * PI / 180.0;
		c = cos(rad);
		s = sin(rad);
		omc = 1 - c;

		d[0] = b->x() * b->x() * omc + c;
		d[1] = b->x() * b->y() * omc + b->z() * s;
		d[2] = b->x() * b->z() * omc - b->y() * s;
		d[4] = b->x() * b->y() * omc - b->z() * s;
		d[5] = b->y() * b->y() * omc + c;
		d[6] = b->y() * b->z() * omc + b->x() * s;
		d[8] = b->x() * b->z() * omc + b->y() * s;
		d[9] = b->y() * b->z() * omc - b->x() * s;
		d[10] = b->z() * b->z() * omc + c;
		d[15] = 1.0;

		a->rotate(angle, *b);
		e = a->constData();


		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(e[j * 4], d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(e[j * 4 + 1], d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(e[j * 4 + 2], d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(e[j * 4 + 3], d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] d;
	}

	std::cout << "Matrix4x4 Rotation Test Finished\n" << std::endl;
}

void Matrix4x4Rotation2Test()
{
	std::cout << "Matrix4x4 Rotation 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* d{ nullptr };
	const double* e{ nullptr };
	double angle{ 0.0 };
	double rad{ 0.0 };
	double c{ 0.0 };
	double s{ 0.0 };
	double omc{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		d = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = 0.0;
		}

		angle = FAUtil::randomRangeD(0, 360);
		rad = angle * PI / 180.0;
		c = cos(rad);
		s = sin(rad);
		omc = 1 - c;

		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		d[0] = x * x * omc + c;
		d[1] = x * y * omc + z * s;
		d[2] = x * z * omc - y * s;
		d[4] = x * y * omc - z * s;
		d[5] = y * y * omc + c;
		d[6] = y * z * omc + x * s;
		d[8] = x * z * omc + y * s;
		d[9] = y * z * omc - x * s;
		d[10] = z * z * omc + c;
		d[15] = 1.0;

		a->rotate(angle, x, y, z);
		e = a->constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(e[j * 4], d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(e[j * 4 + 1], d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(e[j * 4 + 2], d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(e[j * 4 + 3], d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] d;
	}

	std::cout << "Matrix4x4 Rotation 2 Test Finished\n" << std::endl;
}

void Matrix4x4RotationUsingQuaternionTest()
{
	std::cout << "Matrix4x4 Rotation Uisng Quaternion Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector4* b{ nullptr };
	const double* c{ nullptr };
	double* d = { nullptr };
	double angle{ 0.0 };
	double s{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		d = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = 0.0;
		}

		angle = FAUtil::randomRangeD(0, 360);
	
		b = new FAMath::Vector4(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), angle);

		a = new FAMath::Matrix4x4();

		a->rotateUsingQuaternion(*b);
		c = a->constData();

		b->setX(b->x() * sin((angle / 2) * PI / 180.0));
		b->setY(b->y() * sin((angle / 2) * PI / 180.0));
		b->setZ(b->z() * sin((angle / 2) * PI / 180.0));
		b->setW(cos((angle / 2) * PI / 180.0));

		d[0] = 1 - 2 * b->y() * b->y() - 2 * b->z() * b->z();
		d[1] = 2 * b->x() * b->y() + 2 * b->w() * b->z();
		d[2] = 2 * b->x() * b->z() - 2 * b->w() * b->y();

		d[4] = 2 * b->x() * b->y() - 2 * b->w() * b->z();
		d[5] = 1 - 2 * b->x() * b->x() - 2 * b->z() * b->z();
		d[6] = 2 * b->y() * b->z() + 2 * b->w() * b->x();

		d[8] = 2 * b->x() * b->z() + 2 * b->w() * b->y();
		d[9] = 2 * b->y() * b->z() - 2 * b->w() * b->x();
		d[10] = 1 - 2 * b->x() * b->x() - 2 * b->y() * b->y();
		d[15] = 1.0;

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j * 4], d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 1], d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 2], d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 3], d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] d;
	}

	std::cout << "Matrix4x4 Rotation Uisng Quaternion Test Finished\n" << std::endl;
}

void Matrix4x4RotationUsingQuaternion2Test()
{
	std::cout << "Matrix4x4 Rotation Uisng Quaternion 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	const double* c{ nullptr };
	double* d = { nullptr };
	double angle{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		d = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = 0.0;
		}

		angle = FAUtil::randomRangeD(0, 360);

		b = new FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));

		a = new FAMath::Matrix4x4();

		a->rotateUsingQuaternion(angle, *b);
		c = a->constData();

		b->setX(b->x() * sin((angle / 2) * PI / 180.0));
		b->setY(b->y() * sin((angle / 2) * PI / 180.0));
		b->setZ(b->z() * sin((angle / 2) * PI / 180.0));
		w = cos((angle / 2) * PI / 180.0);

		d[0] = 1 - 2 * b->y() * b->y() - 2 * b->z() * b->z();
		d[1] = 2 * b->x() * b->y() + 2 * w * b->z();
		d[2] = 2 * b->x() * b->z() - 2 * w * b->y();

		d[4] = 2 * b->x() * b->y() - 2 * w * b->z();
		d[5] = 1 - 2 * b->x() * b->x() - 2 * b->z() * b->z();
		d[6] = 2 * b->y() * b->z() + 2 * w * b->x();

		d[8] = 2 * b->x() * b->z() + 2 * w * b->y();
		d[9] = 2 * b->y() * b->z() - 2 * w * b->x();
		d[10] = 1 - 2 * b->x() * b->x() - 2 * b->y() * b->y();
		d[15] = 1.0;

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j * 4], d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 1], d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 2], d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 3], d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] d;
	}

	std::cout << "Matrix4x4 Rotation Uisng Quaternion 2 Test Finished\n" << std::endl;
}

void Matrix4x4RotationUsingQuaternion3Test()
{
	std::cout << "Matrix4x4 Rotation Uisng Quaternion 3 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	const double* c{ nullptr };
	double* d = { nullptr };
	double angle{ 0.0 };
	double w{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		d = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			d[i] = 0.0;
		}

		angle = FAUtil::randomRangeD(0, 360);

		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		a = new FAMath::Matrix4x4();

		a->rotateUsingQuaternion(angle, x, y, z);
		c = a->constData();

		x = x * sin((angle / 2) * PI / 180.0);
		y = y * sin((angle / 2) * PI / 180.0);
		z = z * sin((angle / 2) * PI / 180.0);
		w = cos((angle / 2) * PI / 180.0);

		d[0] = 1 - 2 * y * y - 2 * z * z;
		d[1] = 2 * x * y + 2 * w * z;
		d[2] = 2 * x * z - 2 * w * y;

		d[4] = 2 * x * y - 2 * w * z;
		d[5] = 1 - 2 * x * x - 2 * z * z;
		d[6] = 2 * y * z + 2 * w * x;

		d[8] = 2 * x * z + 2 * w * y;
		d[9] = 2 * y * z - 2 * w * x;
		d[10] = 1 - 2 * x * x - 2 * y * y;
		d[15] = 1.0;

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j * 4], d[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 1], d[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 2], d[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j * 4 + 3], d[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] d;
	}

	std::cout << "Matrix4x4 Rotation Uisng Quaternion 3 Test Finished\n" << std::endl;
}

void Matrix4x4ScaleTest()
{
	std::cout << "Matrix4x4 Scale Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double* c{ nullptr };
	const double* d{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		b = new FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));

		c = new double[16]{ b->x(), 0, 0, 0, 0, b->y(), 0, 0, 0, 0, b->z(), 0, 0, 0, 0, 1 };

		a->scale(*b);
		d = a->constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j], c[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 4], c[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 8], c[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(d[j + 12], c[j * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
		delete[] c;
	}

	std::cout << "Matrix4x4 Scale Test Finished\n" << std::endl;
}

void Matrix4x4Scale2Test()
{
	std::cout << "Matrix4x4 Scale 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();

		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		b = new double[16]{ x, 0, 0, 0, 0, y, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

		a->scale(x, y);
		c = a->constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j], b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 4], b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 8], b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 12], b[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Scale 2 Test Finished\n" << std::endl;
}

void Matrix4x4Scale3Test()
{
	std::cout << "Matrix4x4 Scale 3 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();

		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b = new double[16]{ x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1 };

		a->scale(x, y, z);
		c = a->constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j], b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 4], b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 8], b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 12], b[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Scale 3 Test Finished\n" << std::endl;
}

void Matrix4x4Scale4Test()
{
	std::cout << "Matrix4x4 Scale 4 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };
	double factor{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();

		factor = FAUtil::randomRangeD(MIN, MAX);

		b = new double[16]{ factor, 0, 0, 0, 0, factor, 0, 0, 0, 0, factor, 0, 0, 0, 0, 1 };

		a->scale(factor);
		c = a->constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(c[j], b[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 4], b[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 8], b[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[j + 12], b[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Scale 4 Test Finished\n" << std::endl;
}

void Matrix4x4Scale5Test()
{
	std::cout << "Matrix4x4 Scale 5 Test Start\n" << std::endl;


	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double* c{ nullptr };
	const double* d{ nullptr };
	double factor{ 0.0 };
	double fmone{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		b = new FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));

		factor = FAUtil::randomRangeD(MIN, MAX);
		fmone = factor - 1;

		c = new double[16];
		for (unsigned int j = 0; j < 16; ++j)
		{
			c[j] = 0.0;
		}

		c[0] = 1 + fmone * b->x() * b->x();
		c[1] = fmone * b->x() * b->y();
		c[2] = fmone * b->x() * b->z();
		c[4] = fmone * b->x() * b->y();
		c[5] = 1 + fmone * b->y() * b->y();
		c[6] = fmone * b->y() * b->z();
		c[8] = fmone * b->x() * b->z();
		c[9] = fmone * b->y() * b->z();
		c[10] = 1 + fmone * b->z() * b->z();
		c[15] = 1.0;

		a->scale(*b, factor);
		d = a->constData();

		for (unsigned int j = 0; j < 4; ++j)
		{
			assert(FAUtil::compareDoubles(d[j * 4], c[j * 4], EPSILON));
			assert(FAUtil::compareDoubles(d[j * 4 + 1], c[j * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(d[j * 4 + 2], c[j * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(d[j * 4 + 3], c[j * 4 + 3], EPSILON));
		}

		delete a;
		delete[] c;
	}

	std::cout << "Matrix4x4 Scale 5 Test Finished\n" << std::endl;
}

void Matrix4x4OthroTest()
{
	std::cout << "Matrix4x4 Ortho Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };
	double left{ 0.0 };
	double right{ 0.0 };
	double top{ 0.0 };
	double bottom{ 0.0 };
	double near{ 0.0 };
	double far{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		left = FAUtil::randomRangeD(MIN, MAX);
		right = FAUtil::randomRangeD(MIN, MAX);
		top = FAUtil::randomRangeD(MIN, MAX);
		bottom = FAUtil::randomRangeD(MIN, MAX);
		near = FAUtil::randomRangeD(MIN, MAX);
		far = FAUtil::randomRangeD(MIN, MAX);

		while (FAUtil::compareDoubles(left, right, EPSILON))
		{
			left = FAUtil::randomRangeD(MIN, MAX);
		}
		while (FAUtil::compareDoubles(bottom, top, EPSILON))
		{
			top = FAUtil::randomRangeD(MIN, MAX);
		}
		while (FAUtil::compareDoubles(near, far, EPSILON))
		{
			near = FAUtil::randomRangeD(MIN, MAX);
		}

		b = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = 0.0;
		}

		b[0] = 2 / (right - left);
		b[5] = 2 / (top - bottom);
		b[10] = -2 / (far - near);
		b[3] = -(right + left) / (right - left);
		b[7] = -(top + bottom) / (top - bottom);
		b[11] = -(far + near) / (far - near);
		b[15] = 1.0;

		a->ortho(left, right, bottom, top, near, far);
		c = a->constData();
		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(c[i], b[i * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[i + 4], b[i * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[i + 8], b[i * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[i + 12], b[i * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;
	}


	std::cout << "Matrix4x4 Ortho Test Finished\n" << std::endl;
}

void Matrix4x4PerspectiveTest()
{
	std::cout << "Matrix4x4 Perspective Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	const double* c{ nullptr };
	double  fov{ 0.0 };
	double ar{ 0.0 };
	double near{ 0.0 };
	double far{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();

		b = new double[16]{ 0.0 };
		
		fov = FAUtil::randomRangeD(30, 60);
		ar = FAUtil::randomRangeD(0.1, 1.8);
		near = FAUtil::randomRangeD(5, 90);
		far = FAUtil::randomRangeD(100, 1000);

		b[0] = (1 / ar) * (1 / tan((fov / 2) * PI / 180.0));
		b[5] = 1 / tan((fov / 2) * PI / 180.0);
		b[10] = far / (far - near);
		b[11] = 1.0;
		b[14] = -(far * near) / (far - near);

		a->perspective(fov, ar, near, far);
		c = a->constData();

		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(c[i * 4], b[i * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[i * 4 + 1], b[i * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[i * 4 + 2], b[i * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[i * 4 + 3], b[i * 4 + 3], EPSILON));
		}

		delete a;
		delete[] b;

	}

	std::cout << "Matrix4x4 Perspective Test Finished\n" << std::endl;
}

void MatrixDeterminantTest()
{
	std::cout << "Matrix4x4 Determinant Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		a = new FAMath::Matrix4x4(b);
		FAMath::Matrix4x4 c;
		c = a->transposed();
		double t{ c.determinant() };
		double r{ a->determinant() };

		if (!FAUtil::compareDoubles(t, r, EPSILON))
		{
			std::cout << std::setprecision(20) << fabs(t - r) << std::endl;
			std::cout << std::setprecision(20) << "Determinant = " << r << std::endl;
			std::cout << std::setprecision(20) << "Transposed Determinant = " << t << std::endl;
		}

		assert(FAUtil::compareDoubles(t, r, EPSILON));

		a->setToIdentity();
		double id = a->determinant();
		if (!FAUtil::compareDoubles(id, 1.0, EPSILON))
		{
			print(*a);
			std::cout << id << std::endl;
		}
		assert(FAUtil::compareDoubles(id, 1.0, EPSILON));

		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Determinant Test Finished\n" << std::endl;
}

void MatrixInverseTest()
{
	std::cout << "Matrix4x4 Inverse Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double* b{ nullptr };
	FAMath::Matrix4x4 inv;
	FAMath::Matrix4x4 t;
	FAMath::Matrix4x4 c;
	FAMath::Matrix4x4 d;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = new double[16];
		for (unsigned int i = 0; i < 16; ++i)
		{
			b[i] = FAUtil::randomRangeD(MIN, MAX);
		}

		d = inverse(d);
		assert(d.isIdentity());

		a = new FAMath::Matrix4x4(b);

		inv = inverse(*a);
		c = inv * *a;
		if (!c.isIdentity())
		{
			std::cout << "C = " << std::endl; print(c); std::cout << std::endl;
		}
		assert(c.isIdentity());
		
		c = *a * inv;
		if (!c.isIdentity())
		{
			std::cout << "C = " << std::endl; print(c); std::cout << std::endl;
		}
		assert(c.isIdentity());

		double invD = inv.determinant();
		double aD = a->determinant();
		aD = 1.0 / aD;
		if (!FAUtil::compareDoubles(aD, invD, EPSILON))
		{
			std::cout << std::setprecision(20) << std::endl;
			std::cout << "aD = " << aD << std::endl;
			std::cout << "invD = " << invD << std::endl;
		}
		assert(FAUtil::compareDoubles(aD, invD, EPSILON));


		inv = inv.transposed();
		t = a->transposed();
		c = inverse(t);
		if (c != inv)
		{
			std::cout << "C = " << std::endl; print(c); std::cout << std::endl;
			std::cout << "inv = " << std::endl; print(inv); std::cout << std::endl;
		}
		assert(c == inv);


		delete a;
		delete[] b;
	}

	std::cout << "Matrix4x4 Inverse Test Finished\n" << std::endl;
}

void Matrix4x4TranslationTest()
{
	std::cout << "Matrix4x4 Translation Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	FAMath::Vector3* b{ nullptr };
	double c[16];
	const double* d{ nullptr };

	for (unsigned int i = 0; i < 16; ++i)
	{
		c[i] = 0.0;
	}

	c[0] = 1.0;
	c[5] = 1.0;
	c[10] = 1.0;
	c[15] = 1.0;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		b = new FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));

		c[12] = b->x();
		c[13] = b->y();
		c[14] = b->z();

		a->translate(*b);
		d = a->constData();

		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(c[i * 4], d[i * 4], EPSILON));
			assert(FAUtil::compareDoubles(c[i * 4 + 1], d[i * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(c[i * 4 + 2], d[i * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(c[i * 4 + 3], d[i * 4 + 3], EPSILON));
		}

		delete a;
		delete b;
	}

	std::cout << "Matrix4x4 Translation Test Finished\n" << std::endl;
}
void Matrix4x4Translation2Test()
{
	std::cout << "Matrix4x4 Translation 2 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double b[16];
	const double* c{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };

	for (unsigned int i = 0; i < 16; ++i)
	{
		b[i] = 0.0;
	}

	b[0] = 1.0;
	b[5] = 1.0;
	b[10] = 1.0;
	b[15] = 1.0;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);

		b[12] = x;
		b[13] = y;

		a->translate(x, y);
		c = a->constData();

		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(b[i * 4], c[i * 4], EPSILON));
			assert(FAUtil::compareDoubles(b[i * 4 + 1], c[i * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(b[i * 4 + 2], c[i * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(b[i * 4 + 3], c[i * 4 + 3], EPSILON));
		}

		delete a;
	}

	std::cout << "Matrix4x4 Translation 2 Test Finished\n" << std::endl;
}
void Matrix4x4Translation3Test()
{
	std::cout << "Matrix4x4 Translation 3 Test Start\n" << std::endl;

	FAMath::Matrix4x4* a{ nullptr };
	double b[16];
	const double* c{ nullptr };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < 16; ++i)
	{
		b[i] = 0.0;
	}

	b[0] = 1.0;
	b[5] = 1.0;
	b[10] = 1.0;
	b[15] = 1.0;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Matrix4x4();
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);

		b[12] = x;
		b[13] = y;
		b[14] = z;

		a->translate(x, y, z);
		c = a->constData();

		for (unsigned int i = 0; i < 4; ++i)
		{
			assert(FAUtil::compareDoubles(b[i * 4], c[i * 4], EPSILON));
			assert(FAUtil::compareDoubles(b[i * 4 + 1], c[i * 4 + 1], EPSILON));
			assert(FAUtil::compareDoubles(b[i * 4 + 2], c[i * 4 + 2], EPSILON));
			assert(FAUtil::compareDoubles(b[i * 4 + 3], c[i * 4 + 3], EPSILON));
		}

		delete a;
	}

	std::cout << "Matrix4x4 Translation 3 Test Finished\n" << std::endl;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------------------


//QUATERNION TESTING
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
void QuaternionDefaultConstructor()
{
	std::cout << "Quaternion Default Constructor Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Quaternion();
		assert(FAUtil::compareDoubles(a->scalar(), 1.0, EPSILON));
		assert(FAUtil::compareDoubles(a->vector().x(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->vector().y(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->vector().z(), 0.0, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Default Constructor Test Finished\n" << std::endl;
}

void QuaternionOverloadedConstructor1Test()
{
	std::cout << "Quaternion Overloaded Constructor 1 Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double w{ 0.0 };
	FAMath::Vector3 b;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		b = FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		a = new FAMath::Quaternion(w, b);
		assert(FAUtil::compareDoubles(a->scalar(), w, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), b.x(), EPSILON));
		assert(FAUtil::compareDoubles(a->y(), b.y(), EPSILON));
		assert(FAUtil::compareDoubles(a->z(), b.z(), EPSILON));
		delete a;
	}

	std::cout << "Quaternion Overloaded Constructor 1 Test Finished\n" << std::endl;
}

void QuaternionOverloadedConstructor2Test()
{
	std::cout << "Quaternion Overloaded Constructor 2 Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion(w, x, y, z);
		assert(FAUtil::compareDoubles(a->scalar(), w, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Overloaded Constructor 2 Test Finished\n" << std::endl;
}

void QuaternionSetQuaternion1Test()
{
	std::cout << "Quaternion Set Quaternion 1 Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double w{ 0.0 };
	FAMath::Vector3 b;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		b = FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		a = new FAMath::Quaternion();
		a->setQuaternion(w, b);
		assert(FAUtil::compareDoubles(a->scalar(), w, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), b.x(), EPSILON));
		assert(FAUtil::compareDoubles(a->y(), b.y(), EPSILON));
		assert(FAUtil::compareDoubles(a->z(), b.z(), EPSILON));
		delete a;
	}

	std::cout << "Quaternion Set Quaternion 1 Test Finished\n" << std::endl;
}

void QuaternionSetQuaternion2Test()
{
	std::cout << "Quaternion Set Quaternion 2 Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion();
		a->setQuaternion(w, x, y, z);
		assert(FAUtil::compareDoubles(a->scalar(), w, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Set Quaternion 2 Test Finished\n" << std::endl;
}

void QuaternionSetScalarTest()
{
	std::cout << "Quaternion Set Scalar Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double w{ 0.0 };
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion();
		a->setScalar(w);
		assert(FAUtil::compareDoubles(a->scalar(), w, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));

		delete a;
	}

	std::cout << "Quaternion Set Scalar Test Finished\n" << std::endl;
}

void QuaternionSetVector1Test()
{
	std::cout << "Quaternion Set Vector 1 Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	FAMath::Vector3 b;
	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		b = FAMath::Vector3(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		a = new FAMath::Quaternion();
		a->setVector(b);
		assert(FAUtil::compareDoubles(a->scalar(), 1.0, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), b.x(), EPSILON));
		assert(FAUtil::compareDoubles(a->y(), b.y(), EPSILON));
		assert(FAUtil::compareDoubles(a->z(), b.z(), EPSILON));

		delete a;
	}

	std::cout << "Quaternion Set Vector 1 Test Finished\n" << std::endl;
}

void QuaternionSetVector2Test()
{
	std::cout << "Quaternion Set Vector 2 Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion();
		a->setVector(x, y, z);
		assert(FAUtil::compareDoubles(a->scalar(), 1.0, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Set Vector 2 Test Finished\n" << std::endl;
}

void QuaternionSetXTest()
{
	std::cout << "Quaternion Set X Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double x{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		x = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion();
		a->setX(x);
		assert(FAUtil::compareDoubles(a->scalar(), 1.0, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), x, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Set X Test Finished\n" << std::endl;
}

void QuaternionSetYTest()
{
	std::cout << "Quaternion Set Y Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double y{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		y = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion();
		a->setY(y);
		assert(FAUtil::compareDoubles(a->scalar(), 1.0, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->y(), y, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), 0.0, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Set Y Test Finished\n" << std::endl;
}

void QuaternionSetZTest()
{
	std::cout << "Quaternion Set Z Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion();
		a->setZ(z);
		assert(FAUtil::compareDoubles(a->scalar(), 1.0, EPSILON));
		assert(FAUtil::compareDoubles(a->x(), 0.0 , EPSILON));
		assert(FAUtil::compareDoubles(a->y(), 0.0, EPSILON));
		assert(FAUtil::compareDoubles(a->z(), z, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Set Z Test Finished\n" << std::endl;
}

void QuaternionNegationTest()
{
	std::cout << "Quaternion Negation Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	FAMath::Quaternion b;
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion(w, x, y, z);
		b = -(*a);
		assert(FAUtil::compareDoubles(b.scalar(), -a->scalar(), EPSILON));
		assert(FAUtil::compareDoubles(b.x(), -a->x(), EPSILON));
		assert(FAUtil::compareDoubles(b.y(), -a->y(), EPSILON));
		assert(FAUtil::compareDoubles(b.z(), -a->z(), EPSILON));
		delete a;
	}

	std::cout << "Quaternion Negation Test Finished\n" << std::endl;
}

void QuaternionZeroQuaterionTest()
{
	std::cout << "Quaternion Zero Quaternion Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Quaternion();
		assert(!a->isZeroQuaternion());
		a->setScalar(0.0);
		assert(a->isZeroQuaternion());
		delete a;
	}

	std::cout << "Quaternion Zero Quaternion Test Finished\n" << std::endl;
}

void QuaternionLengthTest()
{
	std::cout << "Quaternion Length Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double mag{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion(w, x, y, z);
		mag = length(*a);
		assert(FAUtil::compareDoubles(mag, sqrt(a->scalar() * a->scalar() + a->x() * a->x() + a->y() * a->y() + a->z() * a->z()), EPSILON));
		delete a;
	}

	std::cout << "Quaternion Length Test Finished\n" << std::endl;
}

void QuaternionNormalizeTest()
{
	std::cout << "Quaternion Normalize Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	FAMath::Quaternion b;
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double mag{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Quaternion(0.0, 0.0, 0.0, 0.0);
		b = normalize(*a);
		assert(b.isZeroQuaternion());
		delete a;

		w = FAUtil::randomRangeD(MIN, MAX) + 1.0;
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion(w, x, y, z);
		*a = normalize(*a);
		mag = length(*a);
		assert(FAUtil::compareDoubles(mag, 1.0, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Normalize Test Finished\n" << std::endl;
}

void QuaternionConjugateTest()
{
	std::cout << "Quaternion Conjugate Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	FAMath::Quaternion b;
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion(w, x, y, z);
		b = conjugate(*a);
		assert(FAUtil::compareDoubles(b.scalar(), a->scalar(), EPSILON));
		assert(FAUtil::compareDoubles(b.x(), -a->x(), EPSILON));
		assert(FAUtil::compareDoubles(b.y(), -a->y(), EPSILON));
		assert(FAUtil::compareDoubles(b.z(), -a->z(), EPSILON));
		delete a;
	}

	std::cout << "Quaternion Conjugate Test Finished\n" << std::endl;
}

void QuaternionInverseTest()
{
	std::cout << "Quaternion Inverse Test Start\n" << std::endl;

	FAMath::Quaternion* a = nullptr;
	FAMath::Quaternion b;
	FAMath::Quaternion c;
	double mag{ 0.0 };
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		a = new FAMath::Quaternion(0.0, 0.0, 0.0, 0.0);
		b = conjugate(*a);
		assert(b.isZeroQuaternion());
		delete a;

		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		a = new FAMath::Quaternion(w, x, y, z);
		b = conjugate(*a);
		mag = length(*a);
		c = inverse(*a);
		assert(FAUtil::compareDoubles(c.scalar(), b.scalar() / mag, EPSILON));
		assert(FAUtil::compareDoubles(c.x(), b.x() / mag, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), b.y() / mag, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), b.z() / mag, EPSILON));
		delete a;
	}

	std::cout << "Quaternion Inverse Test Finished\n" << std::endl;
}

void QuaternionMultiplication1Test()
{
	std::cout << "Quaternion Multiplication 1 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;
	FAMath::Quaternion c;
	FAMath::Quaternion d;

	double mag1{ 0.0 };
	double mag2{ 0.0 };
	double w1{ 0.0 };
	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w2{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w1, x1, y1, z1);
		b.setQuaternion(w2, x2, y2, z2);

		d = b * a;
		c = inverse(b) * inverse(a);
		mag1 = length(a) * length(b);
		a *= b;
		mag2 = length(a);
		a = inverse(a);

		assert(a != d);
		assert(c == a);
		assert(FAUtil::compareDoubles(mag1, mag2, EPSILON));
	}


	std::cout << "Quaternion Multiplication 1 Test Finished\n" << std::endl;
}

void QuaternionMultiplication2Test()
{
	std::cout << "Quaternion Multiplication 2 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;
	FAMath::Quaternion c;
	FAMath::Quaternion d;

	double mag1{ 0.0 };
	double mag2{ 0.0 };
	double w1{ 0.0 };
	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w2{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w1, x1, y1, z1);
		b.setQuaternion(w2, x2, y2, z2);

		c = a * b;
		d = b * a;

		mag1 = length(a) * length(b);
		mag2 = length(c);

		assert(c != d);
		assert(FAUtil::compareDoubles(mag1, mag2, EPSILON));

		c = inverse(c);
		d = inverse(b) * inverse(a);

		assert(c == d);
	}


	std::cout << "Quaternion Multiplication 2 Test Finished\n" << std::endl;
}

void QuaternionAddition1Test()
{
	std::cout << "Quaternion Addition 1 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;

	double w1{ 0.0 };
	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w2{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w1, x1, y1, z1);
		b.setQuaternion(w2, x2, y2, z2);

		a += b;

		assert(FAUtil::compareDoubles(a.scalar(), w1 + w2, EPSILON));
		assert(FAUtil::compareDoubles(a.x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(a.y(), y1 + y2, EPSILON));
		assert(FAUtil::compareDoubles(a.z(), z1 + z2, EPSILON));
	}

	std::cout << "Quaternion Addition 1 Test Finished\n" << std::endl;
}

void QuaternionAddition2Test()
{
	std::cout << "Quaternion Addition 2 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;
	FAMath::Quaternion c;

	double w1{ 0.0 };
	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w2{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w1, x1, y1, z1);
		b.setQuaternion(w2, x2, y2, z2);

		c = a + b;

		assert(FAUtil::compareDoubles(c.scalar(), w1 + w2, EPSILON));
		assert(FAUtil::compareDoubles(c.x(), x1 + x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 + y2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z1 + z2, EPSILON));
	}

	std::cout << "Quaternion Addition 2 Test Finished\n" << std::endl;
}

void QuaternionSubtraction1Test()
{
	std::cout << "Quaternion Subtraction 1 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;

	double w1{ 0.0 };
	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w2{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w1, x1, y1, z1);
		b.setQuaternion(w2, x2, y2, z2);

		a -= b;

		assert(FAUtil::compareDoubles(a.scalar(), w1 - w2, EPSILON));
		assert(FAUtil::compareDoubles(a.x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(a.y(), y1 - y2, EPSILON));
		assert(FAUtil::compareDoubles(a.z(), z1 - z2, EPSILON));
	}

	std::cout << "Quaternion Subraction 1 Test Finished\n" << std::endl;
}

void QuaternionSubtraction2Test()
{
	std::cout << "Quaternion Subtraction 2 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;
	FAMath::Quaternion c;

	double w1{ 0.0 };
	double x1{ 0.0 };
	double y1{ 0.0 };
	double z1{ 0.0 };
	double w2{ 0.0 };
	double x2{ 0.0 };
	double y2{ 0.0 };
	double z2{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w1 = FAUtil::randomRangeD(MIN, MAX);
		x1 = FAUtil::randomRangeD(MIN, MAX);
		y1 = FAUtil::randomRangeD(MIN, MAX);
		z1 = FAUtil::randomRangeD(MIN, MAX);
		w2 = FAUtil::randomRangeD(MIN, MAX);
		x2 = FAUtil::randomRangeD(MIN, MAX);
		y2 = FAUtil::randomRangeD(MIN, MAX);
		z2 = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w1, x1, y1, z1);
		b.setQuaternion(w2, x2, y2, z2);

		c = a - b;

		assert(FAUtil::compareDoubles(c.scalar(), w1 - w2, EPSILON));
		assert(FAUtil::compareDoubles(c.x(), x1 - x2, EPSILON));
		assert(FAUtil::compareDoubles(c.y(), y1 - y2, EPSILON));
		assert(FAUtil::compareDoubles(c.z(), z1 - z2, EPSILON));
	}

	std::cout << "Quaternion Subtraction 2 Test Finished\n" << std::endl;
}


void QuaternionMultiplicationByAScalar1Test()
{
	std::cout << "Quaternion Multiplication By A Scalar 1 Test Start\n" << std::endl;

	FAMath::Quaternion a;

	double k{ 0.0 };
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		k = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w, x, y, z);

		a *= k;

		assert(FAUtil::compareDoubles(a.scalar(), k * w, EPSILON));
		assert(FAUtil::compareDoubles(a.x(), k * x, EPSILON));
		assert(FAUtil::compareDoubles(a.y(), k * y, EPSILON));
		assert(FAUtil::compareDoubles(a.z(), k * z, EPSILON));
	
	}

	std::cout << "Quaternion Multiplication By A Scalar 1 Test Finished\n" << std::endl;
}

void QuaternionMultiplicationByAScalar2Test()
{
	std::cout << "Quaternion Multiplication By A Scalar 2 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;

	double k{ 0.0 };
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		k = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w, x, y, z);

		b = a * k;

		assert(FAUtil::compareDoubles(b.scalar(), k * w, EPSILON));
		assert(FAUtil::compareDoubles(b.x(), k * x, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), k * y, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), k * z, EPSILON));

	}

	std::cout << "Quaternion Multiplication By A Scalar 2 Test Finished\n" << std::endl;
}

void QuaternionMultiplicationByAScalar3Test()
{
	std::cout << "Quaternion Multiplication By A Scalar 3 Test Start\n" << std::endl;

	FAMath::Quaternion a;
	FAMath::Quaternion b;

	double k{ 0.0 };
	double w{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		w = FAUtil::randomRangeD(MIN, MAX);
		x = FAUtil::randomRangeD(MIN, MAX);
		y = FAUtil::randomRangeD(MIN, MAX);
		z = FAUtil::randomRangeD(MIN, MAX);
		k = FAUtil::randomRangeD(MIN, MAX);

		a.setQuaternion(w, x, y, z);

		b = k * a;

		assert(FAUtil::compareDoubles(b.scalar(), k * w, EPSILON));
		assert(FAUtil::compareDoubles(b.x(), k * x, EPSILON));
		assert(FAUtil::compareDoubles(b.y(), k * y, EPSILON));
		assert(FAUtil::compareDoubles(b.z(), k * z, EPSILON));

	}

	std::cout << "Quaternion Multiplication By A Scalar 3 Test Finished\n" << std::endl;
}

void QuaternionToRotationMatrixTest()
{
	std::cout << "Quaternion To Rotation Matrix Test Start\n" << std::endl;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		FAMath::Vector4 point(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		double ang = FAUtil::randomRangeD(0.0, 360.0);
		FAMath::Vector3 axis(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		axis = normalize(axis);

		FAMath::Matrix4x4 rot1;
		rot1.rotate(ang, axis);

		ang = (ang / 2.0f) * PI / 180.0;
		axis = axis * sin(ang);
		double w = cos(ang);
		FAMath::Quaternion q(w, axis);
		normalize(q);
		FAMath::Matrix4x4 rot2 = q.toRotationMatrix();

		FAMath::Vector4 r1 = rot1 * point;
		FAMath::Vector4 r2 = rot2 * point;

		if (r1 != r2)
		{
			std::cout << i << std::endl;
			std::cout << "r1 = "; print(r1);
			std::cout << "r2 = "; print(r2);
		}
		assert(r1 == r2);
	}

	std::cout << "Quaternion To Rotation Matrix Test Finished\n" << std::endl;

}

void QuaternionDotProductTest()
{
	std::cout << "Quaternion Dot Product Test Start\n" << std::endl;

	for (unsigned int i = 0; i < LOOPNUM; ++i)
	{
		FAMath::Quaternion a(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));
		FAMath::Quaternion b(FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX), FAUtil::randomRangeD(MIN, MAX));;
		FAMath::Quaternion c;

		a = normalize(a);
		b = normalize(b);
		c = b * inverse(a);

		assert(FAUtil::compareDoubles(c.scalar(), dotProduct(a, b), EPSILON));
	}

	std::cout << "Quaternion Dot Product Test Finished\n" << std::endl;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------