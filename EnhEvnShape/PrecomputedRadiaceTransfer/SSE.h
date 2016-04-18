#pragma once

#include <smmintrin.h>


#if defined(WIN32)

#include <intrin.h>

typedef unsigned __int64 uint64;
typedef __int64 int64;

__forceinline uint64 __rdpmc(int i) {
	return __readpmc(i);
}

__forceinline int __bsf(int v) {
	unsigned long r = 0; _BitScanForward(&r,v); return r;
}

__forceinline int __bsr(int v) {
	unsigned long r = 0; _BitScanReverse(&r,v); return r;
}

__forceinline int __btc(int v, int i) {
	long r = v; _bittestandcomplement(&r,i); return r;
}

__forceinline int __bts(int v, int i) {
	long r = v; _bittestandset(&r,i); return r;
}

__forceinline int __btr(int v, int i) {
	long r = v; _bittestandreset(&r,i); return r;
}

#if defined(_WIN64)

__forceinline size_t __bsf(size_t v) {
	size_t r = 0; _BitScanForward64((unsigned long*)&r,v); return r;
}

__forceinline size_t __bsr(size_t v) {
	size_t r = 0; _BitScanReverse64((unsigned long*)&r,v); return r;
}

__forceinline size_t __btc(size_t v, size_t i) {
	//size_t r = v; _bittestandcomplement64((__int64*)&r,i); return r;
	return v ^ (size_t(1) << i); // faster than using intrinsics, as intrinsic goes through memory
}

__forceinline size_t __bts(size_t v, size_t i) {
	__int64 r = v; _bittestandset64(&r,i); return r;
}

__forceinline size_t __btr(size_t v, size_t i) {
	__int64 r = v; _bittestandreset64(&r,i); return r;
}

#endif

#if defined(_WIN64)

typedef __int64 atomic_t;

__forceinline int64 atomic_add(volatile int64* m, const int64 v) {
	return _InterlockedExchangeAdd64(m,v);
}

__forceinline int64 atomic_xchg(volatile int64 *p, int64 v) {
	return _InterlockedExchange64((volatile long long *)p, v);
}

__forceinline int64 atomic_cmpxchg(volatile int64* m, const int64 v, const int64 c) {
	return _InterlockedCompareExchange64(m,v,c);
}

#else

typedef __int32 int32;

typedef int32 atomic_t;

__forceinline int32 atomic_add(volatile int32* p, const int32 v) {
	return _InterlockedExchangeAdd((volatile long*)p,v);
}

__forceinline int32 atomic_xchg(volatile int32 *p, int32 v) {
	return _InterlockedExchange((volatile long*)p, v);
}

__forceinline int32 atomic_cmpxchg(volatile int32* p, const int32 v, const int32 c) {
	return _InterlockedCompareExchange((volatile long*)p,v,c);
}

#endif

#endif


#include <cmath>
#include <utility>
#include <limits>



#pragma warning(disable: 4244)
#pragma warning(disable: 4018)
#pragma warning(disable: 4267)
#pragma warning(disable: 4996)
#pragma warning(disable: 4800)

// C++ support
#include <cstdio>
#include <malloc.h>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <fstream>

// STL include files
#include <cstring>
#include <vector>
using std::vector;
#include <algorithm>
using std::min;
using std::max;
using std::swap;
using std::sort;


namespace EDX
{


	namespace Constants
	{
		static struct Null {
		} EDX_NULL;

		static struct True {
			__forceinline operator bool() const { return true; }
		} EDX_TRUE;

		static struct False {
			__forceinline operator bool() const { return false; }
		} EDX_FALSE;

		static struct Emp {
		} EDX_EMP;

		static struct Full {
		} EDX_FULL;
	}

	namespace Math
	{
		static struct Zero
		{
			__forceinline operator          double    () const { return 0.0; }
			__forceinline operator          float     () const { return 0.0f; }
			__forceinline operator          long long () const { return 0; }
			__forceinline operator unsigned long long () const { return 0; }
			__forceinline operator          long      () const { return 0; }
			__forceinline operator unsigned long      () const { return 0; }
			__forceinline operator          int       () const { return 0; }
			__forceinline operator unsigned int       () const { return 0; }
			__forceinline operator          short     () const { return 0; }
			__forceinline operator unsigned short     () const { return 0; }
			__forceinline operator          char      () const { return 0; }
			__forceinline operator unsigned char      () const { return 0; }
		} EDX_ZERO;

		static struct One
		{
			__forceinline operator          double    () const { return 1.0; }
			__forceinline operator          float     () const { return 1.0f; }
			__forceinline operator          long long () const { return 1; }
			__forceinline operator unsigned long long () const { return 1; }
			__forceinline operator          long      () const { return 1; }
			__forceinline operator unsigned long      () const { return 1; }
			__forceinline operator          int       () const { return 1; }
			__forceinline operator unsigned int       () const { return 1; }
			__forceinline operator          short     () const { return 1; }
			__forceinline operator unsigned short     () const { return 1; }
			__forceinline operator          char      () const { return 1; }
			__forceinline operator unsigned char      () const { return 1; }
		} EDX_ONE;

		static struct NegInf
		{
			__forceinline operator          double    () const { return -std::numeric_limits<double>::infinity(); }
			__forceinline operator          float     () const { return -std::numeric_limits<float>::infinity(); }
			__forceinline operator          long long () const { return std::numeric_limits<long long>::min(); }
			__forceinline operator unsigned long long () const { return std::numeric_limits<unsigned long long>::min(); }
			__forceinline operator          long      () const { return std::numeric_limits<long>::min(); }
			__forceinline operator unsigned long      () const { return std::numeric_limits<unsigned long>::min(); }
			__forceinline operator          int       () const { return std::numeric_limits<int>::min(); }
			__forceinline operator unsigned int       () const { return std::numeric_limits<unsigned int>::min(); }
			__forceinline operator          short     () const { return std::numeric_limits<short>::min(); }
			__forceinline operator unsigned short     () const { return std::numeric_limits<unsigned short>::min(); }
			__forceinline operator          char      () const { return std::numeric_limits<char>::min(); }
			__forceinline operator unsigned char      () const { return std::numeric_limits<unsigned char>::min(); }

		} EDX_NEG_INFINITY;

		static struct PosInf
		{
			__forceinline operator          double    () const { return std::numeric_limits<double>::infinity(); }
			__forceinline operator          float     () const { return std::numeric_limits<float>::infinity(); }
			__forceinline operator          long long () const { return std::numeric_limits<long long>::max(); }
			__forceinline operator unsigned long long () const { return std::numeric_limits<unsigned long long>::max(); }
			__forceinline operator          long      () const { return std::numeric_limits<long>::max(); }
			__forceinline operator unsigned long      () const { return std::numeric_limits<unsigned long>::max(); }
			__forceinline operator          int       () const { return std::numeric_limits<int>::max(); }
			__forceinline operator unsigned int       () const { return std::numeric_limits<unsigned int>::max(); }
			__forceinline operator          short     () const { return std::numeric_limits<short>::max(); }
			__forceinline operator unsigned short     () const { return std::numeric_limits<unsigned short>::max(); }
			__forceinline operator          char      () const { return std::numeric_limits<char>::max(); }
			__forceinline operator unsigned char      () const { return std::numeric_limits<unsigned char>::max(); }
		} EDX_INFINITY, EDX_POS_INFINITY;

		static struct NaN
		{
			__forceinline operator double() const { return std::numeric_limits<double>::quiet_NaN(); }
			__forceinline operator float () const { return std::numeric_limits<float>::quiet_NaN(); }
		} EDX_NAN;

		static struct Epsilon
		{
			__forceinline operator double() const { return std::numeric_limits<double>::epsilon(); }
			__forceinline operator float () const { return std::numeric_limits<float>::epsilon(); }
		} EDX_EPSILON;

		static struct Pi
		{
			__forceinline operator double() const { return 3.14159265358979323846; }
			__forceinline operator float () const { return 3.14159265358979323846f; }
		} EDX_PI;

		static struct InvPi
		{
			__forceinline operator double() const { return 0.31830988618379069122; }
			__forceinline operator float () const { return 0.31830988618379069122f; }
		} EDX_INV_PI;

		static struct TwoPi
		{
			__forceinline operator double() const { return 6.283185307179586232; }
			__forceinline operator float () const { return 6.283185307179586232f; }
		} EDX_TWO_PI;

		static struct PiOverTwo
		{
			__forceinline operator double() const { return 1.57079632679489661923; }
			__forceinline operator float () const { return 1.57079632679489661923f; }
		} EDX_PI_2;

		static struct InvTwoPi
		{
			__forceinline operator double() const { return 0.15915494309189534561; }
			__forceinline operator float () const { return 0.15915494309189534561f; }
		} EDX_INV_2PI;

		static struct FourPi
		{
			__forceinline operator double() const { return 12.566370614359172464; }
			__forceinline operator float () const { return 12.566370614359172464f; }
		} EDX_FOUR_PI;

		static struct PiOverFour
		{
			__forceinline operator double() const { return 0.785398163397448309616; }
			__forceinline operator float () const { return 0.785398163397448309616f; }
		} EDX_PI_4;

		static struct InvFourPi
		{
			__forceinline operator double() const { return 0.079577471545947672804; }
			__forceinline operator float () const { return 0.079577471545947672804f; }
		} EDX_INV_4PI;

		static struct Step {
		} EDX_STEP;
	}

	template<unsigned int N, class T>
	class Vec
	{
	public:
		T v[N];

		inline const T& operator [] (const size_t idx) const { assert(idx < 1); return v[idx]; }
		inline		 T& operator [] (const size_t idx)		 { assert(idx < 1); return v[idx]; }
	};



	template<class T>
	class Vec<2, T>
	{
	public:
		union
		{
			struct { T x, y; };
			struct { T u, v; };
		};

	public:
		//----------------------------------------------------------------------------------------------
		// Constructors, Assignment & Cast Operators
		//----------------------------------------------------------------------------------------------
		inline Vec()
			: x(Math::EDX_ZERO), y(Math::EDX_ZERO) {}

		inline Vec(const Vec& vCopyFrom)
			: x(vCopyFrom.x), y(vCopyFrom.y) {}

		inline Vec(const T& tVal)
			: x(tVal), y(tVal) {}

		inline Vec(const T& tx, const T& ty)
			: x(tx), y(ty) {}

		template<class T1>
		inline Vec(const Vec<2, T1>& vConvertFrom)
			: x(T(vConvertFrom.x)), y(T(vConvertFrom.y)) {}

		template<class T1>
		inline Vec& operator = (const Vec<2, T1>& vOther)
		{
			x = T(vOther.x); y = T(vOther.y);
			return *this;
		}

		~Vec() {}

		inline T Sum() const { return x + y; }
		inline T Product() const { return x * y; }

		inline const T& operator [] (const size_t idx) const { assert(idx < 2); return (&x)[idx]; }
		inline		 T& operator [] (const size_t idx)		 { assert(idx < 2); return (&x)[idx]; }

		//----------------------------------------------------------------------------------------------
		// Unary Operators
		//----------------------------------------------------------------------------------------------
		inline Vec operator + () const { return Vec(+x, +y); }
		inline Vec operator - () const { return Vec(-x, -y); }

		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		inline Vec operator + (const Vec& rhs) const { return Vec(x + rhs.x, y + rhs.y); }
		inline Vec operator - (const Vec& rhs) const { return Vec(x - rhs.x, y - rhs.y); }
		inline Vec operator * (const Vec& rhs) const { return Vec(x * rhs.x, y * rhs.y); }
		inline Vec operator * (const T& rhs) const { return Vec(x * rhs, y * rhs); }
		inline Vec operator / (const Vec& rhs) const { return Vec(x / rhs.x, y / rhs.y); }
		inline Vec operator / (const T& rhs) const { return Vec(x / rhs, y / rhs); }

		//----------------------------------------------------------------------------------------------
		// Assignment Operators
		//----------------------------------------------------------------------------------------------
		inline const Vec& operator += (const Vec& rhs) { x += rhs.x; y += rhs.y; return *this; }
		inline const Vec& operator -= (const Vec& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
		inline const Vec& operator *= (const T& rhs) { x *= rhs; y *= rhs; return *this; }
		inline const Vec& operator /= (const T& rhs) { x /= rhs; y /= rhs; return *this; }
		inline const Vec& operator *= (const Vec& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
		inline const Vec& operator /= (const Vec& rhs) { x /= rhs.x; y /= rhs.y; return *this; }

		//----------------------------------------------------------------------------------------------
		// Comparison Operators
		//----------------------------------------------------------------------------------------------
		inline bool operator == (const Vec& rhs) const { return x == rhs.x && y == rhs.y; }
		inline bool operator != (const Vec& rhs) const { return x != rhs.x || y != rhs.y; }

		static const Vec ZERO;
		static const Vec UNIT_SCALE;
		static const Vec UNIT_X;
		static const Vec UNIT_Y;
		static const Vec UNIT[2];
	};

	template<class T> const Vec<2, T> Vec<2, T>::ZERO(Math::EDX_ZERO);
	template<class T> const Vec<2, T> Vec<2, T>::UNIT_SCALE(Math::EDX_ONE);
	template<class T> const Vec<2, T> Vec<2, T>::UNIT_X(Math::EDX_ONE, Math::EDX_ZERO);
	template<class T> const Vec<2, T> Vec<2, T>::UNIT_Y(Math::EDX_ZERO, Math::EDX_ONE);
	template<class T> const Vec<2, T> Vec<2, T>::UNIT[2] =
	{
		Vec<2, T>(Math::EDX_ONE, Math::EDX_ZERO),
		Vec<2, T>(Math::EDX_ZERO, Math::EDX_ONE)
	};

	//----------------------------------------------------------------------------------------------
	// Binary Operators
	//----------------------------------------------------------------------------------------------
	template<class T> inline Vec<2, T> operator * (const T& lhs, const Vec<2, T>& rhs) { return rhs * lhs; }

	//----------------------------------------------------------------------------------------------
	// Output Operators
	//----------------------------------------------------------------------------------------------
	template<typename T> inline std::ostream& operator << (std::ostream& out, const Vec<2, T>& rhs)
	{
		return out << "(" << rhs.x << ", " << rhs.y << ")";
	}

	typedef Vec<2, float>	Vector2;
	typedef Vec<2, int>		Vector2i;
	typedef Vec<2, bool>	Vector2b;

	template<class T>
	class Vec<3, T>
	{
	public:
		T x, y, z;

	public:
		//----------------------------------------------------------------------------------------------
		// Constructors, Assignment & Cast Operators
		//----------------------------------------------------------------------------------------------
		inline Vec()
			: x(Math::EDX_ZERO), y(Math::EDX_ZERO), z(Math::EDX_ZERO) {}

		inline Vec(const Vec& vCopyFrom)
			: x(vCopyFrom.x), y(vCopyFrom.y), z(vCopyFrom.z) {}



		inline Vec(const T& tVal)
			: x(tVal), y(tVal), z(tVal) {}

		inline Vec(const T& tx, const T& ty, const T& tz)
			: x(tx), y(ty), z(tz) {}

		//		inline Vec(const Vec<2, T>& vCopyFrom, T val = T(Math::EDX_ZERO))
		//			: x(vCopyFrom.x), y(vCopyFrom.y), z(val) {}

		template<class T1>
		inline Vec(const Vec<3, T1>& vConvertFrom)
			: x(T(vConvertFrom.x)), y(T(vConvertFrom.y)), z(T(vConvertFrom.z)) {}

		template<class T1>
		inline Vec(const Vec<2, T1>& vConvertFrom, T1 val)
			: x(T(vConvertFrom.x)), y(T(vConvertFrom.y)), z(T(val)) {}

		template<class T1>
		inline Vec& operator = (const Vec<3, T1>& vOther)
		{
			x = T(vOther.x); y = T(vOther.y); z = T(vOther.z);
			return *this;
		}

		~Vec() {}

		inline T Sum() const { return x + y + z; }
		inline T Product() const { return x * y * z; }

		inline const T& operator [] (const size_t idx) const { assert(idx < 3); return (&x)[idx]; }
		inline		 T& operator [] (const size_t idx)		 { assert(idx < 3); return (&x)[idx]; }

		//----------------------------------------------------------------------------------------------
		// Unary Operators
		//----------------------------------------------------------------------------------------------
		inline Vec operator + () const { return Vec(+x, +y, +z); }
		inline Vec operator - () const { return Vec(-x, -y, -z); }

		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		inline Vec operator + (const Vec& rhs) const { return Vec(x + rhs.x, y + rhs.y, z + rhs.z); }
		inline Vec operator - (const Vec& rhs) const { return Vec(x - rhs.x, y - rhs.y, z - rhs.z); }
		inline Vec operator * (const Vec& rhs) const { return Vec(x * rhs.x, y * rhs.y, z * rhs.z); }
		inline Vec operator * (const T& rhs) const { return Vec(x * rhs, y * rhs, z * rhs); }
		inline Vec operator / (const Vec& rhs) const { return Vec(x / rhs.x, y / rhs.y, z / rhs.z); }
		inline Vec operator / (const T& rhs) const { return Vec(x / rhs, y / rhs, z / rhs); }

		//----------------------------------------------------------------------------------------------
		// Assignment Operators
		//----------------------------------------------------------------------------------------------
		inline const Vec& operator += (const Vec& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
		inline const Vec& operator -= (const Vec& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
		inline const Vec& operator *= (const T& rhs) { x *= rhs; y *= rhs; z *= rhs; return *this; }
		inline const Vec& operator /= (const T& rhs) { x /= rhs; y /= rhs; z /= rhs; return *this; }
		inline const Vec& operator *= (const Vec& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }
		inline const Vec& operator /= (const Vec& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; return *this; }

		//----------------------------------------------------------------------------------------------
		// Comparison Operators
		//----------------------------------------------------------------------------------------------
		inline bool operator == (const Vec& rhs) const { return x == rhs.x && y == rhs.y && z == rhs.z; }
		inline bool operator != (const Vec& rhs) const { return x != rhs.x || y != rhs.y || z != rhs.z; }

		static const Vec ZERO;
		static const Vec UNIT_SCALE;
		static const Vec UNIT_X;
		static const Vec UNIT_Y;
		static const Vec UNIT_Z;
		static const Vec UNIT[3];
	};

	template<class T> const Vec<3, T> Vec<3, T>::ZERO(Math::EDX_ZERO);
	template<class T> const Vec<3, T> Vec<3, T>::UNIT_SCALE(Math::EDX_ONE);
	template<class T> const Vec<3, T> Vec<3, T>::UNIT_X(Math::EDX_ONE, Math::EDX_ZERO, Math::EDX_ZERO);
	template<class T> const Vec<3, T> Vec<3, T>::UNIT_Y(Math::EDX_ZERO, Math::EDX_ONE, Math::EDX_ZERO);
	template<class T> const Vec<3, T> Vec<3, T>::UNIT_Z(Math::EDX_ZERO, Math::EDX_ZERO, Math::EDX_ONE);
	template<class T> const Vec<3, T> Vec<3, T>::UNIT[3] =
	{
		Vec<3, T>(Math::EDX_ONE, Math::EDX_ZERO, Math::EDX_ZERO),
		Vec<3, T>(Math::EDX_ZERO, Math::EDX_ONE, Math::EDX_ZERO),
		Vec<3, T>(Math::EDX_ZERO, Math::EDX_ZERO, Math::EDX_ONE)
	};

	//----------------------------------------------------------------------------------------------
	// Binary Operators
	//----------------------------------------------------------------------------------------------
	template<class T> inline Vec<3, T> operator * (const T& lhs, const Vec<3, T>& rhs) { return rhs * lhs; }

	//----------------------------------------------------------------------------------------------
	// Output Operators
	//----------------------------------------------------------------------------------------------
	template<typename T> inline std::ostream& operator << (std::ostream& out, const Vec<3, T>& rhs)
	{
		return out << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
	}

	template<typename T> inline ::std::istream &operator >> (::std::istream &is, Vec<3, T> &v)
	{
		using namespace ::std;
		char c1 = 0, c2 = 0;

		is >> c1;
		if (c1 == '(') {
			is >> v.x >> ws >> c2;
			if (c2 == ',')
				is >> v.y >> ws >> c2;
			else
				is.setstate(ios::failbit);
			if (c2 == ',')
				is >> v.z >> ws >> c2;
			else
				is.setstate(ios::failbit);
		}

		if (c1 == '(' && c2 != ')')
			is.setstate(ios::failbit);

		return is;
	}








	typedef Vec<3, float>	Vector3;
	typedef Vec<3, int>		Vector3i;
	typedef Vec<3, bool>	Vector3b;

	template<class T>
	class Vec<4, T>
	{
	public:
		T x, y, z, w;

	public:
		//----------------------------------------------------------------------------------------------
		// Constructors, Assignment & Cast Operators
		//----------------------------------------------------------------------------------------------
		inline Vec()
			: x(Math::EDX_ZERO), y(Math::EDX_ZERO), z(Math::EDX_ZERO), w(Math::EDX_ZERO) {}

		inline Vec(const Vec& vCopyFrom)
			: x(vCopyFrom.x), y(vCopyFrom.y), z(vCopyFrom.z), w(vCopyFrom.w) {}

		inline Vec(const T& tVal)
			: x(tVal), y(tVal), z(tVal), w(tVal) {}

		inline Vec(const T& tx, const T& ty, const T& tz, const T& tw)
			: x(tx), y(ty), z(tz), w(tw) {}

		inline Vec(const Vec<2, T>& vCopyFrom, T valZ = T(Math::EDX_ZERO), T valW = T(Math::EDX_ZERO))
			: x(vCopyFrom.x), y(vCopyFrom.y), z(valZ), w(valW) {}

		inline Vec(const Vec<3, T>& vCopyFrom, T val = T(Math::EDX_ZERO))
			: x(vCopyFrom.x), y(vCopyFrom.y), z(vCopyFrom.z), w(val) {}

		template<class T1>
		inline Vec(const Vec<4, T1>& vConvertFrom)
			: x(T(vConvertFrom.x)), y(T(vConvertFrom.y)), z(T(vConvertFrom.z)), w(T(vConvertFrom.w)) {}

		template<class T1>
		inline Vec(const Vec<2, T1>& vConvertFrom, T1 valZ, T1 valW)
			: x(T(vConvertFrom.x)), y(T(vConvertFrom.y)), z(T(valZ)), w(T1(valW)) {}

		template<class T1>
		inline Vec(const Vec<3, T1>& vConvertFrom, T1 valW)
			: x(T(vConvertFrom.x)), y(T(vConvertFrom.y)), z(T(vConvertFrom.z)), w(T1(valW)) {}

		template<class T1>
		inline Vec& operator = (const Vec<4, T1>& vOther)
		{
			x = T(vOther.x); y = T(vOther.y); z = T(vOther.z); w = T(vOther.w);
			return *this;
		}

		~Vec() {}

		inline T Sum() const { return x + y + z + w; }
		inline T Product() const { return x * y * z * w; }
		inline Vec<3, T> xyz() const { return Vec<3, T>(x, y, z); }

		inline const T& operator [] (const size_t idx) const { assert(idx < 4); return (&x)[idx]; }
		inline		 T& operator [] (const size_t idx)		 { assert(idx < 4); return (&x)[idx]; }

		//----------------------------------------------------------------------------------------------
		// Unary Operators
		//----------------------------------------------------------------------------------------------
		inline Vec operator + () const { return Vec(+x, +y, +z, +w); }
		inline Vec operator - () const { return Vec(-x, -y, -z, -w); }

		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		inline Vec operator + (const Vec& rhs) const { return Vec(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w); }
		inline Vec operator - (const Vec& rhs) const { return Vec(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w); }
		inline Vec operator * (const Vec& rhs) const { return Vec(x * rhs.x, y * rhs.y, z * rhs.z, w * rhs.w); }
		inline Vec operator * (const T& rhs) const { return Vec(x * rhs, y * rhs, z * rhs, w * rhs); }
		inline Vec operator / (const Vec& rhs) const { return Vec(x / rhs.x, y / rhs.y, z / rhs.z, w / rhs.w); }
		inline Vec operator / (const T& rhs) const { return Vec(x / rhs, y / rhs, z / rhs, w / rhs); }

		//----------------------------------------------------------------------------------------------
		// Assignment Operators
		//----------------------------------------------------------------------------------------------
		inline const Vec& operator += (const Vec& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this; }
		inline const Vec& operator -= (const Vec& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this; }
		inline const Vec& operator *= (const T& rhs) { x *= rhs; y *= rhs; z *= rhs; w *= rhs; return *this; }
		inline const Vec& operator /= (const T& rhs) { x /= rhs; y /= rhs; z /= rhs; w /= rhs; return *this; }
		inline const Vec& operator *= (const Vec& rhs) { x *= rhs.x; y *= rhs.y; z *= rhs.z; w *= rhs.w; return *this; }
		inline const Vec& operator /= (const Vec& rhs) { x /= rhs.x; y /= rhs.y; z /= rhs.z; w /= rhs.w; return *this; }

		//----------------------------------------------------------------------------------------------
		// Comparison Operators
		//----------------------------------------------------------------------------------------------
		inline bool operator == (const Vec& rhs) const { return x == rhs.x && y == rhs.y && z == rhs.z && w == rhs.w; }
		inline bool operator != (const Vec& rhs) const { return x != rhs.x || y != rhs.y || z != rhs.z || w != rhs.w; }

		static const Vec ZERO;
		static const Vec UNIT_SCALE;
		static const Vec UNIT_X;
		static const Vec UNIT_Y;
		static const Vec UNIT_Z;
		static const Vec UNIT_W;
		static const Vec UNIT[4];
	};

	template<class T> const Vec<4, T> Vec<4, T>::ZERO(Math::EDX_ZERO);
	template<class T> const Vec<4, T> Vec<4, T>::UNIT_SCALE(Math::EDX_ONE);
	template<class T> const Vec<4, T> Vec<4, T>::UNIT_X(Math::EDX_ONE, Math::EDX_ZERO, Math::EDX_ZERO, EDX_ZERO);
	template<class T> const Vec<4, T> Vec<4, T>::UNIT_Y(Math::EDX_ZERO, Math::EDX_ONE, Math::EDX_ZERO, EDX_ZERO);
	template<class T> const Vec<4, T> Vec<4, T>::UNIT_Z(Math::EDX_ZERO, Math::EDX_ZERO, Math::EDX_ONE, EDX_ZERO);
	template<class T> const Vec<4, T> Vec<4, T>::UNIT_W(Math::EDX_ZERO, Math::EDX_ZERO, Math::EDX_ZERO, EDX_ONE);
	template<class T> const Vec<4, T> Vec<4, T>::UNIT[4] =
	{
		Vec<4, T>(Math::EDX_ONE, Math::EDX_ZERO, Math::EDX_ZERO, EDX_ZERO),
		Vec<4, T>(Math::EDX_ZERO, Math::EDX_ONE, Math::EDX_ZERO, EDX_ZERO),
		Vec<4, T>(Math::EDX_ZERO, Math::EDX_ZERO, Math::EDX_ONE, EDX_ZERO),
		Vec<4, T>(Math::EDX_ZERO, Math::EDX_ZERO, Math::EDX_ZERO, EDX_ONE)
	};

	//----------------------------------------------------------------------------------------------
	// Binary Operators
	//----------------------------------------------------------------------------------------------
	template<class T> inline Vec<4, T> operator * (const T& lhs, const Vec<4, T>& rhs) { return rhs * lhs; }

	//----------------------------------------------------------------------------------------------
	// Output Operators
	//----------------------------------------------------------------------------------------------
	template<typename T> inline std::ostream& operator << (std::ostream& out, const Vec<4, T>& rhs)
	{
		return out << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
	}

	typedef Vec<4, float>	Vector4;
	typedef Vec<4, int>		Vector4i;
	typedef Vec<4, bool>	Vector4b;







	typedef unsigned int	 uint;
	typedef unsigned char	 byte;
	typedef unsigned long	 dword;
	typedef			 __int8	 int8;
	typedef unsigned __int8	 uint8;
	typedef			 __int16 int16;
	typedef unsigned __int16 uint16;
	typedef			 __int32 int32;
	typedef unsigned __int32 uint32;
	typedef			 __int64 int64;
	typedef unsigned __int64 uint64;

	class Object
	{
	public:
		virtual ~Object()
		{
		}
	};

	namespace Math
	{
		template<class T>
		class DotProductKernel;
	}



	namespace Math
	{
		template <typename T>
		inline bool NumericValid(T num) { return !_isnan(num) && _finite(num); }

		template <class T1, class T2, class T3>
		inline T1 Clamp(T1 tVal, T2 tMin, T3 tMax)
		{
			if(tVal < tMin) return tMin;
			if(tVal > tMax) return tMax;
			return tVal;
		}
		inline float Saturate(float fVal) { return Clamp(fVal, 0.0f, 1.0f); }
		inline float Pow(float fVal, float fPow) { return powf(fVal, fPow); }
		//inline float Abs(float fVal) { return fabsf(fVal); }
		template <class T> inline T Abs(T tVal) { return tVal >= 0 ? tVal : -tVal; }
		inline float Sqrt(float fVal) { return sqrtf(fVal); }
		inline float ToRadians(float fDeg) { return (fDeg / 180.0f) * float(Math::EDX_PI); }
		inline float ToDegrees(float fRad) { return (fRad / float(Math::EDX_PI)) * 180.0f; }
		inline float Sin(float fVal) { return sinf(fVal); }
		inline float Cos(float fVal) { return cosf(fVal); }
		inline void SinCos(float fVal, float& fSin, float& fCos) { fSin = sinf(fVal); fCos = cosf(fVal); }
		inline float Tan(float fVal) { return tanf(fVal); }
		inline float Exp(float fVal) { return expf(fVal); }
		inline float Atan2(float fVal1, float fVal2) { return atan2f(fVal1, fVal2); }

		inline int FloorToInt(float fVal) { return (int)fVal; }
		inline int CeilToInt(float fVal)
		{
			int iVal = (int)fVal;
			if(fVal - iVal > 1e-5f)
			{
				return iVal + 1;
			}
			else
			{
				return iVal;
			}
		}
		inline int RoundToInt(float fVal) { return FloorToInt(fVal + 0.5f); }

		template <class T1, class T2>
		inline float LinStep(T1 tVal, T2 tMin, T2 tMax)
		{
			if(tMin == tMax)
			{
				return 0.0f;
			}
			return Saturate((tVal - tMin) / (tMax - tMin));
		}
		template<class T1, class T2> inline auto Lerp(T1 eMin, T2 eMax, float fVal) -> decltype(eMin + eMax) { return eMin * (1 - fVal) + eMax * fVal; }
		template<class T1, class T2> inline auto Max(T1 e1, T2 e2) -> decltype(e1 + e2)  { return e1 > e2 ? e1 : e2; }
		template<class T1, class T2> inline auto AbsMax(T1 e1, T2 e2) -> decltype(e1 + e2)  { return Abs(e1) > Abs(e2) ? e1 : e2; }
		template<class T1, class T2> inline auto Min(T1 e1, T2 e2) -> decltype(e1 + e2)  { return e1 < e2 ? e1 : e2; }
		template<class T1, class T2> inline auto AbsMin(T1 e1, T2 e2) -> decltype(e1 + e2)  { return Abs(e1) < Abs(e2) ? e1 : e2; }

		template<class T1, class T2, class T3, class T4>
		inline auto BiLerp(T1 e00, T2 e01, T3 e10, T4 e11, float fVal1, float fVal2) -> decltype(e00 + e01 + e10 + e11)
		{
			return Lerp(Lerp(e00, e01, fVal1), Lerp(e10, e11, fVal1), fVal2);
		}
		template<class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8>
		inline auto TriLerp(T1 e000, T2 e001, T3 e010, T4 e011, T5 e100, T6 e101, T7 e110, T8 e111, float fVal1, float fVal2, float fVal3) -> decltype(e000 + e001 + e010 + e011 + e100 + e101 + e110 + e111)
		{
			return Lerp(BiLerp(e000, e001, e010, e011, fVal1, fVal2), BiLerp(e100, e101, e110, e111, fVal1, fVal2), fVal3);
		}

		inline bool IsPowOfTwo(int iVal)
		{
			return (iVal & (iVal - 1)) == 0;
		}
		inline int RoundUpPowTwo(uint iVal)
		{
			iVal--;
			iVal |= iVal >> 1;
			iVal |= iVal >> 2;
			iVal |= iVal >> 4;
			iVal |= iVal >> 8;
			iVal |= iVal >> 16;
			return iVal + 1;
		}
		inline int RoundUpTo(unsigned int iVal, unsigned int iRound)
		{
			iVal += iRound - (iVal % iRound);
			return iVal;
		}

		template<uint Dim> class Pow2
		{
		public:
			static const int Value = 1 << Dim;
		};
		template<uint Dim> class Pow4
		{
		public:
			static const int Value = 1 << 2 * Dim;
		};

		inline float Xor(const float f1, const float f2)
		{ 
			return _mm_cvtss_f32(_mm_xor_ps(_mm_set_ss(f1), _mm_set_ss(f2)));
		}

		inline float SignMask(const float f)
		{
			return _mm_cvtss_f32(_mm_and_ps(_mm_set_ss(f),_mm_castsi128_ps(_mm_set1_epi32(0x80000000))));
		}

		template<class T>
		inline T MonoCubicLerp(T e0, T e1, T e2, T e3, float fLerp)
		{
			assert(abs(e2) < 1e8f);
			T eVal = Math::Lerp(e1, e2, fLerp);

			T eMinusCur = fLerp;
			T eMinusCurSqr = eMinusCur * eMinusCur;
			T eMinusCurCub = eMinusCurSqr * eMinusCur;

			T eK = e2 - e1;
			T eD1 = (e2 - e0) * 0.5f;
			T eD2 = (e3 - e1) * 0.5f;
			if(eK * eD1 < 0)
			{
				eD1 = 0;
			}
			if(eK * eD2 < 0)
			{
				eD2 = 0;
			}
			if(eK == 0)
			{
				eD1 = eD2 = 0;
			}

			T eA0 = e1;
			T eA1 = eD1;
			T eA2 = 3 * eK - 2 * eD1 - eD2;
			T eA3 = eD1 + eD2 - 2 * eK;
			T eRet = eA3 * eMinusCurCub + eA2 * eMinusCurSqr + eA1 * eMinusCur + eA0;
			//eRet = eA0 * eMinusCurCub + eA1 * eMinusCurSqr + eA2 * eMinusCur + eA3;
			assert(NumericValid(eRet));

			return eRet;
		}
	}

	const __m128 _mm_lookupmask_ps[16] = {
		_mm_castsi128_ps(_mm_set_epi32(0, 0, 0, 0)),
		_mm_castsi128_ps(_mm_set_epi32(0, 0, 0,-1)),
		_mm_castsi128_ps(_mm_set_epi32(0, 0,-1, 0)),
		_mm_castsi128_ps(_mm_set_epi32(0, 0,-1,-1)),
		_mm_castsi128_ps(_mm_set_epi32(0,-1, 0, 0)),
		_mm_castsi128_ps(_mm_set_epi32(0,-1, 0,-1)),
		_mm_castsi128_ps(_mm_set_epi32(0,-1,-1, 0)),
		_mm_castsi128_ps(_mm_set_epi32(0,-1,-1,-1)),
		_mm_castsi128_ps(_mm_set_epi32(-1, 0, 0, 0)),
		_mm_castsi128_ps(_mm_set_epi32(-1, 0, 0,-1)),
		_mm_castsi128_ps(_mm_set_epi32(-1, 0,-1, 0)),
		_mm_castsi128_ps(_mm_set_epi32(-1, 0,-1,-1)),
		_mm_castsi128_ps(_mm_set_epi32(-1,-1, 0, 0)),
		_mm_castsi128_ps(_mm_set_epi32(-1,-1, 0,-1)),
		_mm_castsi128_ps(_mm_set_epi32(-1,-1,-1, 0)),
		_mm_castsi128_ps(_mm_set_epi32(-1,-1,-1,-1))
	};

	// 4-wide SSE bool type.
	class BoolSSE
	{
	public:
		typedef BoolSSE Mask;                    // mask type for us
		enum   { size = 4 };                  // number of SIMD elements
		union  { __m128 m128; int32 v[4]; };  // data

	public:
		//----------------------------------------------------------------------------------------------
		// Constructors, Assignment & Cast Operators
		//----------------------------------------------------------------------------------------------
		__forceinline BoolSSE() {}

		__forceinline BoolSSE(const BoolSSE& copyFrom)
			: m128(copyFrom.m128) {}

		__forceinline BoolSSE& operator = (const BoolSSE& copyFrom)
		{
			m128 = copyFrom.m128;
			return *this;
		}

		__forceinline BoolSSE(const __m128& val)
			: m128(val) {}

		__forceinline operator const __m128&(void) const { return m128; }
		__forceinline operator const __m128i(void) const { return _mm_castps_si128(m128); }
		__forceinline operator const __m128d(void) const { return _mm_castps_pd(m128); }

		__forceinline BoolSSE     (bool  a)
			: m128(_mm_lookupmask_ps[(size_t(a) << 3) | (size_t(a) << 2) | (size_t(a) << 1) | size_t(a)]) {}
		__forceinline BoolSSE     (bool  a, bool  b) 
			: m128(_mm_lookupmask_ps[(size_t(b) << 3) | (size_t(a) << 2) | (size_t(b) << 1) | size_t(a)]) {}
		__forceinline BoolSSE     (bool  a, bool  b, bool  c, bool  d)
			: m128(_mm_lookupmask_ps[(size_t(d) << 3) | (size_t(c) << 2) | (size_t(b) << 1) | size_t(a)]) {}

		//----------------------------------------------------------------------------------------------
		// Constants
		//----------------------------------------------------------------------------------------------
		__forceinline BoolSSE(Constants::False) : m128(_mm_setzero_ps()) {}
		__forceinline BoolSSE(Constants::True) : m128(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128()))) {}

		//----------------------------------------------------------------------------------------------
		// Array Access
		//----------------------------------------------------------------------------------------------
		__forceinline bool operator [] (const size_t i) const { assert(i < 4); return (_mm_movemask_ps(m128) >> i) & 1; }
		__forceinline int32& operator [] (const size_t i) { assert(i < 4); return v[i]; }

		//----------------------------------------------------------------------------------------------
		// Unary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE operator ! () const { return _mm_xor_ps(*this, BoolSSE(Constants::EDX_TRUE)); }

		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE operator & (const BoolSSE& rhs) const { return _mm_and_ps(*this, rhs); }
		__forceinline const BoolSSE operator | (const BoolSSE& rhs) const { return _mm_or_ps (*this, rhs); }
		__forceinline const BoolSSE operator ^ (const BoolSSE& rhs) const { return _mm_xor_ps(*this, rhs); }

		//----------------------------------------------------------------------------------------------
		// Assignment Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE operator &= (const BoolSSE& rhs) { return *this = *this & rhs; }
		__forceinline const BoolSSE operator |= (const BoolSSE& rhs) { return *this = *this | rhs; }
		__forceinline const BoolSSE operator ^= (const BoolSSE& rhs) { return *this = *this ^ rhs; }

		//----------------------------------------------------------------------------------------------
		// Comparison Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE operator != (const BoolSSE& rhs) const { return _mm_xor_ps(*this, rhs); }
		__forceinline const BoolSSE operator == (const BoolSSE& rhs) const { return _mm_castsi128_ps(_mm_cmpeq_epi32(*this, rhs)); }

	};

	namespace SSE
	{
		//----------------------------------------------------------------------------------------------
		// Select
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE Select(const BoolSSE& m, const BoolSSE& t, const BoolSSE& f)
		{
			return _mm_blendv_ps(f, t, m); 
		}

		//----------------------------------------------------------------------------------------------
		// Movement/Shifting/Shuffling Functions
		//----------------------------------------------------------------------------------------------

		__forceinline const BoolSSE UnpackLow(const BoolSSE& lhs, const BoolSSE& rhs) { return _mm_unpacklo_ps(lhs, rhs); }
		__forceinline const BoolSSE UnpackHigh(const BoolSSE& lhs, const BoolSSE& rhs) { return _mm_unpackhi_ps(lhs, rhs); }

		template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const BoolSSE Shuffle(const BoolSSE& lhs)
		{
			return _mm_shuffle_epi32(lhs, _MM_SHUFFLE(i3, i2, i1, i0));
		}

		template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const BoolSSE Shuffle(const BoolSSE& lhs, const BoolSSE& rhs)
		{
			return _mm_Shuffle_ps(lhs, rhs, _MM_SHUFFLE(i3, i2, i1, i0));
		}

		template<> __forceinline const BoolSSE Shuffle<0, 0, 2, 2>(const BoolSSE& lhs) { return _mm_moveldup_ps(lhs); }
		template<> __forceinline const BoolSSE Shuffle<1, 1, 3, 3>(const BoolSSE& lhs) { return _mm_movehdup_ps(lhs); }
		template<> __forceinline const BoolSSE Shuffle<0, 1, 0, 1>(const BoolSSE& lhs) { return _mm_castpd_ps(_mm_movedup_pd (lhs)); }

		template<size_t dst, size_t src, size_t clr> __forceinline const BoolSSE Insert(const BoolSSE& lhs, const BoolSSE& rhs) { return _mm_insert_ps(lhs, rhs, (dst << 4) | (src << 6) | clr); }
		template<size_t dst, size_t src> __forceinline const BoolSSE Insert(const BoolSSE& lhs, const BoolSSE& rhs) { return insert<dst, src, 0>(lhs, rhs); }
		template<size_t dst>             __forceinline const BoolSSE Insert(const BoolSSE& lhs, const bool rhs) { return insert<dst,0>(lhs, BoolSSE(rhs)); }

		//----------------------------------------------------------------------------------------------
		// Reduction Operations
		//----------------------------------------------------------------------------------------------
		__forceinline size_t popcnt(const BoolSSE& lhs) { return bool(lhs[0])+bool(lhs[1])+bool(lhs[2])+bool(lhs[3]); }

		__forceinline bool ReduceAnd(const BoolSSE& lhs) { return _mm_movemask_ps(lhs) == 0xf; }
		__forceinline bool ReduceOr(const BoolSSE& lhs) { return _mm_movemask_ps(lhs) != 0x0; }
		__forceinline bool All(const BoolSSE& rhs) { return _mm_movemask_ps(rhs) == 0xf; }
		__forceinline bool Any(const BoolSSE& rhs) { return _mm_movemask_ps(rhs) != 0x0; }
		__forceinline bool None(const BoolSSE& rhs) { return _mm_movemask_ps(rhs) == 0x0; }

	}
	//----------------------------------------------------------------------------------------------
	// Output Operators
	//----------------------------------------------------------------------------------------------
	inline std::ostream& operator << (std::ostream& out, const BoolSSE& rhs)
	{
		return out << "<" << rhs[0] << ", " << rhs[1] << ", " << rhs[2] << ", " << rhs[3] << ">";
	}

	// 4-wide SSE integer type.
	class IntSSE
	{
	public:
		typedef BoolSSE Mask;          // mask type for us
		enum  { size = 4 };         // number of SIMD elements
		union { __m128i m128; int32 v[4]; }; // data

	public:
		//----------------------------------------------------------------------------------------------
		// Constructors, Assignment & Cast Operators
		//----------------------------------------------------------------------------------------------
		__forceinline IntSSE()
			: m128(_mm_setzero_si128()) {}

		__forceinline IntSSE(const IntSSE& copyFrom)
			: m128(copyFrom.m128) {}

		__forceinline IntSSE& operator = (const IntSSE& copyFrom)
		{
			m128 = copyFrom.m128;
			return *this;
		}

		__forceinline IntSSE(const __m128i& val)
			: m128(val) {}

		__forceinline operator const __m128i&(void) const { return m128; }
		__forceinline operator		 __m128i&(void)		  { return m128; }

		__forceinline explicit IntSSE (const int32* const piVal)
			: m128(_mm_loadu_si128((__m128i*)piVal)) {}

		//__forceinline ssei(int32 a) : m128(_mm_set1_epi32(a)) {}

		__forceinline IntSSE(const int32& a)
			: m128(_mm_shuffle_epi32(_mm_castps_si128(_mm_load_ss((float*)&a)), _MM_SHUFFLE(0, 0, 0, 0))) {}

		__forceinline IntSSE(int32 a, int32 b)
			: m128(_mm_set_epi32(b, a, b, a)) {}

		__forceinline IntSSE(int32 a, int32 b, int32 c, int32 d)
			: m128(_mm_set_epi32(d, c, b, a)) {}

		__forceinline explicit IntSSE(const __m128& val)
			: m128(_mm_cvtps_epi32(val)) {}

		//----------------------------------------------------------------------------------------------
		// Constants
		//----------------------------------------------------------------------------------------------
		__forceinline IntSSE(Math::Zero)		: m128(_mm_setzero_si128()) {}
		__forceinline IntSSE(Math::One)		: m128(_mm_set_epi32(1, 1, 1, 1)) {}
		__forceinline IntSSE(Math::PosInf)		: m128(_mm_set_epi32(Math::EDX_INFINITY,Math::EDX_INFINITY,Math::EDX_INFINITY,Math::EDX_INFINITY)) {}
		__forceinline IntSSE(Math::NegInf)		: m128(_mm_set_epi32(Math::EDX_NEG_INFINITY,Math::EDX_NEG_INFINITY,Math::EDX_NEG_INFINITY,Math::EDX_NEG_INFINITY)) {}
		__forceinline IntSSE(Math::Step)		: m128(_mm_set_epi32(3, 2, 1, 0)) {}

		//----------------------------------------------------------------------------------------------
		// Array Access
		//----------------------------------------------------------------------------------------------
		__forceinline const int32& operator [] (const size_t i) const { assert(i < 4); return v[i]; }
		__forceinline		int32& operator [] (const size_t i)		  { assert(i < 4); return v[i]; }

		//----------------------------------------------------------------------------------------------
		// Unary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const IntSSE operator + () const { return *this; }
		__forceinline const IntSSE operator - () const { return _mm_sub_epi32(_mm_setzero_si128(), m128); }

		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const IntSSE operator + (const IntSSE& rhs) const { return _mm_add_epi32(m128, rhs.m128); }
		__forceinline const IntSSE operator + (const int32& rhs) const { return *this + IntSSE(rhs); }
		__forceinline const IntSSE operator - (const IntSSE& rhs) const { return _mm_sub_epi32(m128, rhs.m128); }
		__forceinline const IntSSE operator - (const int32& rhs) const { return *this - IntSSE(rhs); }
		__forceinline const IntSSE operator * (const IntSSE& rhs) const { return _mm_mullo_epi32(m128, rhs.m128); }
		__forceinline const IntSSE operator * (const int32& rhs) const { return *this * IntSSE(rhs); }
		__forceinline const IntSSE operator & (const IntSSE& rhs) const { return _mm_and_si128(m128, rhs.m128); }
		__forceinline const IntSSE operator & (const int32& rhs) const { return *this & IntSSE(rhs); }
		__forceinline const IntSSE operator | (const IntSSE& rhs) const { return _mm_or_si128(m128, rhs.m128); }
		__forceinline const IntSSE operator | (const int32& rhs) const { return *this | IntSSE(rhs); }
		__forceinline const IntSSE operator ^ (const IntSSE& rhs) const { return _mm_xor_si128(m128, rhs.m128); }
		__forceinline const IntSSE operator ^ (const int32& rhs) const { return *this ^ IntSSE(rhs); }
		__forceinline const IntSSE operator << (const int32& n) const { return _mm_slli_epi32(m128, n); }
		__forceinline const IntSSE operator >> (const int32& n) const { return _mm_srai_epi32(m128, n); }

		//----------------------------------------------------------------------------------------------
		// Assignment Operators
		//----------------------------------------------------------------------------------------------
		__forceinline IntSSE& operator += (const IntSSE& rhs) { return *this = *this + rhs; }
		__forceinline IntSSE& operator += (const int32& rhs) { return *this = *this + rhs; }

		__forceinline IntSSE& operator -= (const IntSSE& rhs) { return *this = *this - rhs; }
		__forceinline IntSSE& operator -= (const int32& rhs) { return *this = *this - rhs; }

		__forceinline IntSSE& operator *= (const IntSSE& rhs) { return *this = *this * rhs; }
		__forceinline IntSSE& operator *= (const int32& rhs) { return *this = *this * rhs; }

		__forceinline IntSSE& operator &= (const IntSSE& rhs) { return *this = *this & rhs; }
		__forceinline IntSSE& operator &= (const int32& rhs) { return *this = *this & rhs; }

		__forceinline IntSSE& operator |= (const IntSSE& rhs) { return *this = *this | rhs; }
		__forceinline IntSSE& operator |= (const int32& rhs) { return *this = *this | rhs; }

		__forceinline IntSSE& operator <<= (const int32& rhs) { return *this = *this << rhs; }
		__forceinline IntSSE& operator >>= (const int32& rhs) { return *this = *this >> rhs; }

		//----------------------------------------------------------------------------------------------
		// Comparison Operators + Select
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE operator == (const IntSSE& rhs) const { return _mm_castsi128_ps(_mm_cmpeq_epi32 (m128, rhs.m128)); }
		__forceinline const BoolSSE operator == (const int32& rhs) const { return *this == IntSSE(rhs); }

		__forceinline const BoolSSE operator != (const IntSSE& rhs) const { return !(*this == rhs); }
		__forceinline const BoolSSE operator != (const int32& rhs) const { return *this != IntSSE(rhs); }

		__forceinline const BoolSSE operator < (const IntSSE& rhs) const { return _mm_castsi128_ps(_mm_cmplt_epi32 (m128, rhs.m128)); }
		__forceinline const BoolSSE operator < (const int32& rhs) const { return *this < IntSSE(rhs); }

		__forceinline const BoolSSE operator >= (const IntSSE& rhs) const { return !(*this < rhs); }
		__forceinline const BoolSSE operator >= (const int32& rhs) const { return *this >= IntSSE(rhs); }

		__forceinline const BoolSSE operator > (const IntSSE& rhs) const { return _mm_castsi128_ps(_mm_cmpgt_epi32 (m128, rhs.m128)); }
		__forceinline const BoolSSE operator > (const int32& rhs) const { return *this > IntSSE(rhs); }

		__forceinline const BoolSSE operator <= (const IntSSE& rhs) const { return !(*this > rhs); }
		__forceinline const BoolSSE operator <= (const int32& rhs) const { return *this <= IntSSE(rhs); }
	};

	//----------------------------------------------------------------------------------------------
	// Binary Operators
	//----------------------------------------------------------------------------------------------
	__forceinline const IntSSE operator + (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) + rhs; }
	__forceinline const IntSSE operator - (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) - rhs; }
	__forceinline const IntSSE operator * (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) * rhs; }
	__forceinline const IntSSE operator & (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) & rhs; }
	__forceinline const IntSSE operator | (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) | rhs; }
	__forceinline const IntSSE operator ^ (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) ^ rhs; }

	__forceinline const BoolSSE operator == (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) == rhs; }
	__forceinline const BoolSSE operator != (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) != rhs; }
	__forceinline const BoolSSE operator < (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) < rhs; }
	__forceinline const BoolSSE operator >= (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) >= rhs; }
	__forceinline const BoolSSE operator > (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) > rhs; }
	__forceinline const BoolSSE operator <= (const int32& lhs, const IntSSE& rhs) { return IntSSE(lhs) <= rhs; }

	namespace SSE
	{
		//----------------------------------------------------------------------------------------------
		// Comparison Operators + Select
		//----------------------------------------------------------------------------------------------
		__forceinline const IntSSE Abs(const IntSSE& lhs) { return _mm_abs_epi32(lhs.m128); }

		__forceinline const IntSSE srlhs (const IntSSE& lhs, const int32& rhs) { return _mm_srai_epi32(lhs.m128, rhs); }
		__forceinline const IntSSE srl (const IntSSE& lhs, const int32& rhs) { return _mm_srli_epi32(lhs.m128, rhs); }

		__forceinline const IntSSE Min(const IntSSE& lhs, const IntSSE& rhs) { return _mm_min_epi32(lhs.m128, rhs.m128); }
		__forceinline const IntSSE Min(const IntSSE& lhs, const int32& rhs) { return Min(lhs, IntSSE(rhs)); }
		__forceinline const IntSSE Min(const int32& lhs, const IntSSE& rhs) { return Min(IntSSE(lhs), rhs); }

		__forceinline const IntSSE Max(const IntSSE& lhs, const IntSSE& rhs) { return _mm_max_epi32(lhs.m128, rhs.m128); }
		__forceinline const IntSSE Max(const IntSSE& lhs, const int32& rhs) { return Max(lhs, IntSSE(rhs)); }
		__forceinline const IntSSE Max(const int32& lhs, const IntSSE& rhs) { return Max(IntSSE(lhs), rhs); }


		__forceinline const IntSSE Select(const BoolSSE& m, const IntSSE& lhs, const IntSSE& rhs)
		{
			return _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(rhs), _mm_castsi128_ps(lhs), m)); 
		}

		//----------------------------------------------------------------------------------------------
		// Movement/Shifting/Shuffling Functions
		//----------------------------------------------------------------------------------------------
		__forceinline IntSSE UnpackLow(const IntSSE& lhs, const IntSSE& rhs) { return _mm_castps_si128(_mm_unpacklo_ps(_mm_castsi128_ps(lhs.m128), _mm_castsi128_ps(rhs.m128))); }
		__forceinline IntSSE UnpackHigh(const IntSSE& lhs, const IntSSE& rhs) { return _mm_castps_si128(_mm_unpackhi_ps(_mm_castsi128_ps(lhs.m128), _mm_castsi128_ps(rhs.m128))); }

		template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const IntSSE Shuffle(const IntSSE& lhs)
		{
			return _mm_shuffle_epi32(lhs, _MM_SHUFFLE(i3, i2, i1, i0));
		}

		template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const IntSSE Shuffle(const IntSSE& lhs, const IntSSE& rhs)
		{
			return _mm_castps_si128(_mm_Shuffle_ps(_mm_castsi128_ps(lhs), _mm_castsi128_ps(rhs), _MM_SHUFFLE(i3, i2, i1, i0)));
		}

		template<> __forceinline const IntSSE Shuffle<0, 0, 2, 2>(const IntSSE& lhs) { return _mm_castps_si128(_mm_moveldup_ps(_mm_castsi128_ps(lhs))); }
		template<> __forceinline const IntSSE Shuffle<1, 1, 3, 3>(const IntSSE& lhs) { return _mm_castps_si128(_mm_movehdup_ps(_mm_castsi128_ps(lhs))); }
		template<> __forceinline const IntSSE Shuffle<0, 1, 0, 1>(const IntSSE& lhs) { return _mm_castpd_si128(_mm_movedup_pd (_mm_castsi128_pd(lhs))); }

		template<size_t src> __forceinline int Extract(const IntSSE& rhs) { return _mm_extract_epi32(rhs, src); }
		template<size_t dst> __forceinline const IntSSE Insert(const IntSSE& lhs, const int32 rhs) { return _mm_insert_epi32(lhs, rhs, dst); }

		//----------------------------------------------------------------------------------------------
		// Reductions
		//----------------------------------------------------------------------------------------------
		__forceinline const IntSSE VReduceMin(const IntSSE& v) { IntSSE h = Min(Shuffle<1,0,3,2>(v),v); return Min(Shuffle<2,3,0,1>(h),h); }
		__forceinline const IntSSE VReduceMax(const IntSSE& v) { IntSSE h = Max(Shuffle<1,0,3,2>(v),v); return Max(Shuffle<2,3,0,1>(h),h); }
		__forceinline const IntSSE vReduceAdd(const IntSSE& v) { IntSSE h = Shuffle<1,0,3,2>(v) + v ; return Shuffle<2,3,0,1>(h) + h ; }

		__forceinline int ReduceMin(const IntSSE& v) { return Extract<0>(VReduceMin(v)); }
		__forceinline int ReduceMax(const IntSSE& v) { return Extract<0>(VReduceMax(v)); }
		__forceinline int ReduceAdd(const IntSSE& v) { return Extract<0>(vReduceAdd(v)); }

		__forceinline size_t SelectMin(const IntSSE& v) { return __bsf(_mm_movemask_ps(v == VReduceMin(v))); }
		__forceinline size_t SelectMax(const IntSSE& v) { return __bsf(_mm_movemask_ps(v == VReduceMax(v))); }

		__forceinline size_t SelectMin(const BoolSSE& valid, const IntSSE& v) { const IntSSE tmp = Select(valid,v,IntSSE(Math::EDX_INFINITY)); return __bsf(_mm_movemask_ps(valid & (tmp == VReduceMin(tmp)))); }
		__forceinline size_t SelectMax(const BoolSSE& valid, const IntSSE& v) { const IntSSE tmp = Select(valid,v,IntSSE(Math::EDX_NEG_INFINITY)); return __bsf(_mm_movemask_ps(valid & (tmp == VReduceMax(tmp)))); }

	}
	//----------------------------------------------------------------------------------------------
	// Output Operators
	//----------------------------------------------------------------------------------------------
	inline std::ostream& operator << (std::ostream& out, const IntSSE& rhs)
	{
		return out << "<" << rhs[0] << ", " << rhs[1] << ", " << rhs[2] << ", " << rhs[3] << ">";
	}


	// 4-wide SSE float type
	class FloatSSE
	{
	public:
		typedef BoolSSE Mask;          // mask type for us
		typedef IntSSE Int ;			// int type
		enum  { size = 4 };         // number of SIMD elements
		union { __m128 m128; float v[4]; int i[4]; }; // data

	public:
		//----------------------------------------------------------------------------------------------
		// Constructors, Assignment & Cast Operators
		//----------------------------------------------------------------------------------------------
		__forceinline FloatSSE() 
			: m128(_mm_setzero_ps()) {}

		__forceinline FloatSSE(const FloatSSE& copyFrom)
			: m128(copyFrom.m128) {}

		__forceinline FloatSSE& operator = (const FloatSSE& copyFrom)
		{
			m128 = copyFrom.m128;
			return *this;
		}

		__forceinline FloatSSE(const __m128& val)
			: m128(val) {}

		__forceinline operator const __m128&(void) const { return m128; }
		__forceinline operator		 __m128&(void)		 { return m128; }

		__forceinline explicit FloatSSE(const float* const pfVal)
			: m128(_mm_loadu_ps(pfVal)) {}

		//__forceinline ssef(float rhs) : m128(_mm_set1_ps(rhs)) {}

		__forceinline FloatSSE (const float& fVal)
			: m128(_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(_mm_load_ss(&fVal)), _MM_SHUFFLE(0, 0, 0, 0)))) {}

		__forceinline FloatSSE(float a, float b)
			: m128(_mm_set_ps(b, a, b, a)) {}

		__forceinline FloatSSE(float a, float b, float c, float d)
			: m128(_mm_set_ps(d, c, b, a)) {}

		__forceinline explicit FloatSSE(const __m128i rhs)
			: m128(_mm_cvtepi32_ps(rhs)) {}

		//----------------------------------------------------------------------------------------------
		// Constants
		//----------------------------------------------------------------------------------------------
		__forceinline FloatSSE(Math::Zero)		: m128(_mm_setzero_ps()) {}
		__forceinline FloatSSE(Math::One)		: m128(_mm_set_ps(1.0f,1.0f,1.0f,1.0f)) {}
		__forceinline FloatSSE(Math::PosInf)	: m128(_mm_set_ps(Math::EDX_INFINITY,Math::EDX_INFINITY,Math::EDX_INFINITY,Math::EDX_INFINITY)) {}
		__forceinline FloatSSE(Math::NegInf)	: m128(_mm_set_ps(Math::EDX_NEG_INFINITY,Math::EDX_NEG_INFINITY,Math::EDX_NEG_INFINITY,Math::EDX_NEG_INFINITY)) {}
		__forceinline FloatSSE(Math::Step)		: m128(_mm_set_ps(3.0f, 2.0f, 1.0f, 0.0f)) {}
		__forceinline FloatSSE(Math::NaN)		: m128(_mm_set1_ps(Math::EDX_NAN)) {}

		//----------------------------------------------------------------------------------------------
		// Array Access
		//----------------------------------------------------------------------------------------------
		__forceinline const float& operator [] (const size_t i) const { assert(i < 4); return v[i]; }
		__forceinline		float& operator [] (const size_t i)		  { assert(i < 4); return v[i]; }

		//----------------------------------------------------------------------------------------------
		// Unary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const FloatSSE operator + () const { return *this; }
		__forceinline const FloatSSE operator - () const { return _mm_xor_ps(m128, _mm_castsi128_ps(_mm_set1_epi32(0x80000000))); }

		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const FloatSSE operator + (const FloatSSE& rhs) const { return _mm_add_ps(m128, rhs.m128); }
		__forceinline const FloatSSE operator + (const float& rhs) const { return *this + FloatSSE(rhs); }
		__forceinline const FloatSSE operator - (const FloatSSE& rhs) const { return _mm_sub_ps(m128, rhs.m128); }
		__forceinline const FloatSSE operator - (const float& rhs) const { return *this - FloatSSE(rhs); }
		__forceinline const FloatSSE operator * (const FloatSSE& rhs) const { return _mm_mul_ps(m128, rhs.m128); }
		__forceinline const FloatSSE operator * (const float& rhs) const { return *this * FloatSSE(rhs); }
		__forceinline const FloatSSE operator / (const FloatSSE& rhs) const { return *this * _mm_rcp_ps(rhs); }
		__forceinline const FloatSSE operator / (const float& rhs) const { return *this * (1.0f / rhs); }
		__forceinline const FloatSSE operator ^ (const FloatSSE& rhs) const { return _mm_xor_ps(m128, rhs.m128); }
		__forceinline const FloatSSE operator ^ (const IntSSE& rhs) const { return _mm_xor_ps(m128, _mm_castsi128_ps(rhs.m128)); }

		//----------------------------------------------------------------------------------------------
		// Assignment Operators
		//----------------------------------------------------------------------------------------------
		__forceinline FloatSSE& operator += (const FloatSSE& rhs) { return *this = *this + rhs; }
		__forceinline FloatSSE& operator += (const float& rhs) { return *this = *this + rhs; }

		__forceinline FloatSSE& operator -= (const FloatSSE& rhs) { return *this = *this - rhs; }
		__forceinline FloatSSE& operator -= (const float& rhs) { return *this = *this - rhs; }

		__forceinline FloatSSE& operator *= (const FloatSSE& rhs) { return *this = *this * rhs; }
		__forceinline FloatSSE& operator *= (const float& rhs) { return *this = *this * rhs; }

		__forceinline FloatSSE& operator /= (const FloatSSE& rhs) { return *this = *this / rhs; }
		__forceinline FloatSSE& operator /= (const float& rhs) { return *this = *this / rhs; }

		//----------------------------------------------------------------------------------------------
		// Comparison Operators + Select
		//----------------------------------------------------------------------------------------------
		__forceinline const BoolSSE operator == (const FloatSSE& rhs) const { return _mm_cmpeq_ps (m128, rhs.m128); }
		__forceinline const BoolSSE operator == (const float& rhs) const { return *this == FloatSSE(rhs); }
		__forceinline const BoolSSE operator != (const FloatSSE& rhs) const { return _mm_cmpneq_ps(m128, rhs.m128); }
		__forceinline const BoolSSE operator != (const float& rhs) const { return *this != FloatSSE(rhs); }
		__forceinline const BoolSSE operator < (const FloatSSE& rhs) const { return _mm_cmplt_ps (m128, rhs.m128); }
		__forceinline const BoolSSE operator < (const float& rhs) const { return *this < FloatSSE(rhs); }
		__forceinline const BoolSSE operator >= (const FloatSSE& rhs) const { return _mm_cmpnlt_ps(m128, rhs.m128); }
		__forceinline const BoolSSE operator > (const FloatSSE& rhs) const { return _mm_cmpnle_ps(m128, rhs.m128); }
		__forceinline const BoolSSE operator > (const float& rhs) const { return *this > FloatSSE(rhs); }
		__forceinline const BoolSSE operator >= (const float& rhs) const { return *this >= FloatSSE(rhs); }
		__forceinline const BoolSSE operator <= (const FloatSSE& rhs) const { return _mm_cmple_ps (m128, rhs.m128); }
		__forceinline const BoolSSE operator <= (const float& rhs) const { return *this <= FloatSSE(rhs); }

	};

	//----------------------------------------------------------------------------------------------
	// Comparison Operators
	//----------------------------------------------------------------------------------------------
	__forceinline const BoolSSE operator == (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) == rhs; }
	__forceinline const BoolSSE operator != (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) != rhs; }
	__forceinline const BoolSSE operator < (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) < rhs; }
	__forceinline const BoolSSE operator >= (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) >= rhs; }
	__forceinline const BoolSSE operator > (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) > rhs; }
	__forceinline const BoolSSE operator <= (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) <= rhs; }

	namespace SSE
	{
		__forceinline const FloatSSE Abs(const FloatSSE& rhs) { return _mm_and_ps(rhs.m128, _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff))); }
		__forceinline const FloatSSE Sign(const FloatSSE& rhs) { return _mm_blendv_ps(FloatSSE(Math::EDX_ONE), -FloatSSE(Math::EDX_ONE), _mm_cmplt_ps (rhs,FloatSSE(Math::EDX_ZERO))); }
		__forceinline const FloatSSE SignMask(const FloatSSE& rhs) { return _mm_and_ps(rhs.m128,_mm_castsi128_ps(_mm_set1_epi32(0x80000000))); }

		__forceinline const FloatSSE Rcp(const FloatSSE& rhs)
		{
			const FloatSSE r = _mm_rcp_ps(rhs.m128);
			return _mm_sub_ps(_mm_add_ps(r, r), _mm_mul_ps(_mm_mul_ps(r, r), rhs));
		}
		__forceinline const FloatSSE Sqr(const FloatSSE& rhs) { return _mm_mul_ps(rhs,rhs); }
		__forceinline const FloatSSE Sqrt(const FloatSSE& rhs) { return _mm_sqrt_ps(rhs.m128); }
		__forceinline const FloatSSE Rsqrt(const FloatSSE& rhs)
		{
			const FloatSSE r = _mm_rsqrt_ps(rhs.m128);
			return _mm_add_ps(_mm_mul_ps(_mm_set_ps(1.5f, 1.5f, 1.5f, 1.5f), r),
				_mm_mul_ps(_mm_mul_ps(_mm_mul_ps(rhs, _mm_set_ps(-0.5f, -0.5f, -0.5f, -0.5f)), r), _mm_mul_ps(r, r)));}


		//----------------------------------------------------------------------------------------------
		// Binary Operators
		//----------------------------------------------------------------------------------------------
		__forceinline const FloatSSE operator + (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) + rhs; }
		__forceinline const FloatSSE operator - (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) - rhs; }
		__forceinline const FloatSSE operator * (const float& lhs, const FloatSSE& rhs) { return FloatSSE(lhs) * rhs; }
		__forceinline const FloatSSE operator / (const float& lhs, const FloatSSE& rhs) { return lhs * Rcp(rhs); }

		__forceinline const FloatSSE Min(const FloatSSE& lhs, const FloatSSE& rhs) { return _mm_min_ps(lhs.m128, rhs.m128); }
		__forceinline const FloatSSE Min(const FloatSSE& lhs, const float& rhs) { return _mm_min_ps(lhs.m128, FloatSSE(rhs)); }
		__forceinline const FloatSSE Min(const float& lhs, const FloatSSE& rhs) { return _mm_min_ps(FloatSSE(lhs), rhs.m128); }

		__forceinline const FloatSSE Max(const FloatSSE& lhs, const FloatSSE& rhs) { return _mm_max_ps(lhs.m128, rhs.m128); }
		__forceinline const FloatSSE Max(const FloatSSE& lhs, const float& rhs) { return _mm_max_ps(lhs.m128, FloatSSE(rhs)); }
		__forceinline const FloatSSE Max(const float& lhs, const FloatSSE& rhs) { return _mm_max_ps(FloatSSE(lhs), rhs.m128); }


		//----------------------------------------------------------------------------------------------
		// Select
		//----------------------------------------------------------------------------------------------
		__forceinline const FloatSSE Select(const BoolSSE& mask, const FloatSSE& t, const FloatSSE& f)
		{ 
			return _mm_blendv_ps(f, t, mask); 
		}

		//----------------------------------------------------------------------------------------------
		// Rounding Functions
		//----------------------------------------------------------------------------------------------
		__forceinline const FloatSSE RoundEven(const FloatSSE& lhs) { return _mm_round_ps(lhs, _MM_FROUND_TO_NEAREST_INT); }
		__forceinline const FloatSSE RoundDown(const FloatSSE& lhs) { return _mm_round_ps(lhs, _MM_FROUND_TO_NEG_INF  ); }
		__forceinline const FloatSSE RoundUp(const FloatSSE& lhs) { return _mm_round_ps(lhs, _MM_FROUND_TO_POS_INF  ); }
		__forceinline const FloatSSE RoundZero(const FloatSSE& lhs) { return _mm_round_ps(lhs, _MM_FROUND_TO_ZERO   ); }
		__forceinline const FloatSSE Floor(const FloatSSE& lhs) { return _mm_round_ps(lhs, _MM_FROUND_TO_NEG_INF  ); }
		__forceinline const FloatSSE Ceil(const FloatSSE& lhs) { return _mm_round_ps(lhs, _MM_FROUND_TO_POS_INF  ); }

		//----------------------------------------------------------------------------------------------
		// Movement/Shifting/Shuffling Functions
		//----------------------------------------------------------------------------------------------
		__forceinline FloatSSE UnpackLow(const FloatSSE& lhs, const FloatSSE& rhs) { return _mm_unpacklo_ps(lhs.m128, rhs.m128); }
		__forceinline FloatSSE UnpackHigh(const FloatSSE& lhs, const FloatSSE& rhs) { return _mm_unpackhi_ps(lhs.m128, rhs.m128); }

		template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const FloatSSE Shuffle(const FloatSSE& rhs)
		{
			return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(rhs), _MM_SHUFFLE(i3, i2, i1, i0)));
		}

		template<size_t i0, size_t i1, size_t i2, size_t i3> __forceinline const FloatSSE Shuffle(const FloatSSE& lhs, const FloatSSE& rhs)
		{
			return _mm_shuffle_ps(lhs, rhs, _MM_SHUFFLE(i3, i2, i1, i0));
		}

		__forceinline const FloatSSE Shuffle8(const FloatSSE& lhs, const IntSSE& shuf)
		{ 
			return _mm_castsi128_ps(_mm_shuffle_epi8(_mm_castps_si128(lhs), shuf)); 
		}

		template<> __forceinline const FloatSSE Shuffle<0, 0, 2, 2>(const FloatSSE& rhs) { return _mm_moveldup_ps(rhs); }
		template<> __forceinline const FloatSSE Shuffle<1, 1, 3, 3>(const FloatSSE& rhs) { return _mm_movehdup_ps(rhs); }
		template<> __forceinline const FloatSSE Shuffle<0, 1, 0, 1>(const FloatSSE& rhs) { return _mm_castpd_ps(_mm_movedup_pd(_mm_castps_pd(rhs))); }

		template<size_t i>
		__forceinline float Extract(const FloatSSE& lhs) { return _mm_cvtss_f32(Shuffle<i,i,i,i>(lhs)); }
		template<size_t dst, size_t src, size_t clr>
		__forceinline const FloatSSE Insert(const FloatSSE& lhs, const FloatSSE& rhs) { return _mm_insert_ps(lhs, rhs, (dst << 4) | (src << 6) | clr); }
		template<size_t dst, size_t src>
		__forceinline const FloatSSE Insert(const FloatSSE& lhs, const FloatSSE& rhs) { return insert<dst, src, 0>(lhs, rhs); }
		template<size_t dst>
		__forceinline const FloatSSE Insert(const FloatSSE& lhs, const float rhs) { return insert<dst, 0>(lhs, _mm_set_ss(rhs)); }

		//----------------------------------------------------------------------------------------------
		// Transpose
		//----------------------------------------------------------------------------------------------
		__forceinline void Transpose(const FloatSSE& r0, const FloatSSE& r1, const FloatSSE& r2, const FloatSSE& r3, FloatSSE& c0, FloatSSE& c1, FloatSSE& c2, FloatSSE& c3)
		{
			FloatSSE l02 = UnpackLow(r0,r2);
			FloatSSE h02 = UnpackHigh(r0,r2);
			FloatSSE l13 = UnpackLow(r1,r3);
			FloatSSE h13 = UnpackHigh(r1,r3);
			c0 = UnpackLow(l02,l13);
			c1 = UnpackHigh(l02,l13);
			c2 = UnpackLow(h02,h13);
			c3 = UnpackHigh(h02,h13);
		}

		__forceinline void Transpose(const FloatSSE& r0, const FloatSSE& r1, const FloatSSE& r2, const FloatSSE& r3, FloatSSE& c0, FloatSSE& c1, FloatSSE& c2)
		{
			FloatSSE l02 = UnpackLow(r0,r2);
			FloatSSE h02 = UnpackHigh(r0,r2);
			FloatSSE l13 = UnpackLow(r1,r3);
			FloatSSE h13 = UnpackHigh(r1,r3);
			c0 = UnpackLow(l02,l13);
			c1 = UnpackHigh(l02,l13);
			c2 = UnpackLow(h02,h13);
		}

		//----------------------------------------------------------------------------------------------
		// Reductions
		//----------------------------------------------------------------------------------------------
		__forceinline const FloatSSE VReduceMin(const FloatSSE& v) { FloatSSE h = Min(Shuffle<1,0,3,2>(v),v); return Min(Shuffle<2,3,0,1>(h),h); }
		__forceinline const FloatSSE VReduceMax(const FloatSSE& v) { FloatSSE h = Max(Shuffle<1,0,3,2>(v),v); return Max(Shuffle<2,3,0,1>(h),h); }
		__forceinline const FloatSSE vReduceAdd(const FloatSSE& v) { FloatSSE h = Shuffle<1,0,3,2>(v) + v ; return Shuffle<2,3,0,1>(h) + h ; }

		__forceinline float ReduceMin(const FloatSSE& v) { return _mm_cvtss_f32(VReduceMin(v)); }
		__forceinline float ReduceMax(const FloatSSE& v) { return _mm_cvtss_f32(VReduceMax(v)); }
		__forceinline float ReduceAdd(const FloatSSE& v) { return _mm_cvtss_f32(vReduceAdd(v)); }

		__forceinline size_t SelectMin(const FloatSSE& v) { return __bsf(_mm_movemask_ps(v == VReduceMin(v))); }
		__forceinline size_t SelectMax(const FloatSSE& v) { return __bsf(_mm_movemask_ps(v == VReduceMax(v))); }

		__forceinline size_t SelectMin(const BoolSSE& valid, const FloatSSE& v) { const FloatSSE tmp = Select(valid, v, FloatSSE(Math::EDX_INFINITY)); return __bsf(_mm_movemask_ps(valid & (tmp == VReduceMin(tmp)))); }
		__forceinline size_t SelectMax(const BoolSSE& valid, const FloatSSE& v) { const FloatSSE tmp = Select(valid, v, FloatSSE(Math::EDX_NEG_INFINITY)); return __bsf(_mm_movemask_ps(valid & (tmp == VReduceMax(tmp)))); }

	}

	//----------------------------------------------------------------------------------------------
	// Output Operators
	//----------------------------------------------------------------------------------------------
	inline std::ostream& operator << (std::ostream& out, const FloatSSE& rhs)
	{
		return out << "<" << rhs[0] << ", " << rhs[1] << ", " << rhs[2] << ", " << rhs[3] << ">";
	}

	inline ::std::istream &operator >> (::std::istream &is, FloatSSE& rhs)
	{
		using namespace ::std;
		char c1 = 0, c2 = 0;
		is >> c1;
		if (c1 == '<') {

			is >> rhs[0] >> ws >> c2;
			for (size_t i = 1; i < 4; i++) {
				if (c2 == ',')
					is >> rhs[i] >> ws >> c2;
				else
					is.setstate(ios::failbit);
			}
		}
		else
			is.setstate(ios::failbit);

		if (c1 == '<' && c2 != '>')
			is.setstate(ios::failbit);

		return is;
	}










	typedef Vec<2, FloatSSE> Vec2f_SSE;
	typedef Vec<2, IntSSE> Vec2i_SSE;
	typedef Vec<2, BoolSSE> Vec2b_SSE;

	typedef Vec<3, FloatSSE> Vec3f_SSE;
	typedef Vec<3, IntSSE> Vec3i_SSE;
	typedef Vec<3, BoolSSE> Vec3b_SSE;

	typedef Vec<4, FloatSSE> Vec4f_SSE;
	typedef Vec<4, IntSSE> Vec4i_SSE;
	typedef Vec<4, BoolSSE> Vec4b_SSE;





	template<class T>
	inline void SafeDelete(T*& pPtr)
	{ 
		if(pPtr != NULL) 
		{
			delete pPtr;
			pPtr = NULL;
		}
	}

	template<class T>
	inline void SafeDeleteArray(T*& pPtr) 
	{
		if(pPtr != NULL)
		{
			delete [] pPtr;
			pPtr = NULL;
		}
	}

	template<class T>
	inline void SafeClear(T* pPtr, size_t iSize)
	{
		if(pPtr != NULL)
		{
			memset(pPtr, 0, sizeof(T) * iSize);
		}
	}

	template<class T> T* AllocAligned(uint uiCount, uint uiAligned = 64/*for 64-byte cache line*/)
	{
		return (T*)_aligned_malloc(uiCount * sizeof(T), uiAligned);
	}

	template<class T> void FreeAligned(T*& pPtr)
	{
		if(pPtr)
		{
			_aligned_free(pPtr);
			pPtr = NULL;
		}
	}

	class MemoryArena
	{
	private:
		uint muiBlockSize;
		uint muiCurrOffset;
		byte* mpCurrentBlock;

		vector<byte*> mvUsedBlocks, mvAvailableBlocks;

	public:
		MemoryArena(uint uiSize = 32768)
		{
			muiBlockSize = uiSize;
			muiCurrOffset = 0;
			mpCurrentBlock = AllocAligned<byte>(muiBlockSize);
		}
		~MemoryArena()
		{
			// Free all memories in the destructor
			FreeAligned(mpCurrentBlock);

			for(uint i = 0; i < mvUsedBlocks.size(); i++)
			{
				FreeAligned(mvUsedBlocks[i]);
			}

			for(uint i = 0; i < mvAvailableBlocks.size(); i++)
			{
				FreeAligned(mvAvailableBlocks[i]);
			}
		}

		template<class T>
		inline T* Alloc(uint uiCount = 1)
		{
			// Calculate the size of requested memory in terms of byte
			uint uiSize = uiCount * sizeof(T);

			// Make it aligned to 16 byte
			uiSize = (uiSize + 15) & (~15);

			// Handle situation where the current block is used up
			if(muiCurrOffset + uiSize > muiBlockSize)
			{
				// Cache the current block for future use
				mvUsedBlocks.push_back(mpCurrentBlock);

				// Use previously allocated blocks if the requested size is within block size
				if(mvAvailableBlocks.size() > 0 && uiSize < muiBlockSize)
				{
					mpCurrentBlock = mvAvailableBlocks.back();
					mvAvailableBlocks.pop_back();
				}
				else // Else allocate new block in requested size
				{
					mpCurrentBlock = AllocAligned<byte>(Math::Max(uiSize, muiBlockSize));
				}
				// Clear the current block offset
				muiCurrOffset = 0;
			}

			T* pRet = (T*)(mpCurrentBlock + muiCurrOffset);
			muiCurrOffset += uiSize;

			return pRet;
		}

		inline void FreeAll()
		{
			muiCurrOffset = 0;
			while(mvUsedBlocks.size())
			{
				mvAvailableBlocks.push_back(mvUsedBlocks.back());
				mvUsedBlocks.pop_back();
			}
		}
	};


}




