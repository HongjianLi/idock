#include <cmath>
#include <cassert>
#include "array.hpp"

size_t mr(const size_t x, const size_t y)
{
	assert(x <= y);
	return (y*(y+1)>>1) + x;
}

size_t mp(const size_t x, const size_t y)
{
	return x <= y ? mr(x, y) : mr(y, x);
}

float norm_sqr(const array<float, 3>& a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

float norm_sqr(const float a0, const float a1, const float a2)
{
	return a0 * a0 + a1 * a1 + a2 * a2;
}

float norm_sqr(const array<float, 4>& a)
{
	return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3];
}

float norm_sqr(const float a0, const float a1, const float a2, const float a3)
{
	return a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
}

float norm(const array<float, 3>& a)
{
	return sqrt(norm_sqr(a));
}

float norm(const float a0, const float a1, const float a2)
{
	return sqrt(norm_sqr(a0, a1, a2));
}

float norm(const array<float, 4>& a)
{
	return sqrt(norm_sqr(a));
}

float norm(const float a0, const float a1, const float a2, const float a3)
{
	return sqrt(norm_sqr(a0, a1, a2, a3));
}

bool normalized(const array<float, 3>& a)
{
	return fabs(norm_sqr(a) - 1.0f) < 1e-2f;
}

bool normalized(const float a0, const float a1, const float a2)
{
	return fabs(norm_sqr(a0, a1, a2) - 1.0f) < 1e-2f;
}

//! Returns true if the current quaternion is normalized.
bool normalized(const array<float, 4>& a)
{
	return fabs(norm_sqr(a) - 1.0f) < 1e-2f;
}

bool normalized(const float a0, const float a1, const float a2, const float a3)
{
	return fabs(norm_sqr(a0, a1, a2, a3) - 1.0f) < 1e-2f;
}

array<float, 3> normalize(const array<float, 3>& a)
{
	const float norm_inv = 1.0f / norm(a);
	return
	{
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv,
	};
}

void normalize(float &a0, float &a1, float &a2)
{
	const float norm_inv = 1.0f / norm(a0, a1, a2);
	a0 *= norm_inv;
	a1 *= norm_inv;
	a2 *= norm_inv;
}

array<float, 4> normalize(const array<float, 4>& a)
{
	const float norm_inv = 1.0f / norm(a);
	return
	{
		a[0] * norm_inv,
		a[1] * norm_inv,
		a[2] * norm_inv,
		a[3] * norm_inv,
	};
}

void normalize(float &a0, float &a1, float &a2, float &a3)
{
	const float norm_inv = 1.0f / norm(a0, a1, a2, a3);
	a0 *= norm_inv;
	a1 *= norm_inv;
	a2 *= norm_inv;
	a3 *= norm_inv;
}

array<float, 3> operator+(const array<float, 3>& a, const array<float, 3>& b)
{
	return
	{
		a[0] + b[0],
		a[1] + b[1],
		a[2] + b[2],
	};
}

array<float, 3> operator-(const array<float, 3>& a, const array<float, 3>& b)
{
	return
	{
		a[0] - b[0],
		a[1] - b[1],
		a[2] - b[2],
	};
}

void operator+=(array<float, 3>& a, const array<float, 3>& b)
{
	a[0] += b[0];
	a[1] += b[1];
	a[2] += b[2];
}

void operator-=(array<float, 3>& a, const array<float, 3>& b)
{
	a[0] -= b[0];
	a[1] -= b[1];
	a[2] -= b[2];
}

array<float, 3> operator*(const float s, const array<float, 3>& a)
{
	return
	{
		s * a[0],
		s * a[1],
		s * a[2],
	};
}

array<float, 3> operator*(const array<float, 3>& a, const array<float, 3>& b)
{
	return
	{
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0],
	};
}

float distance_sqr(const array<float, 3>& a, const array<float, 3>& b)
{
	const float d0 = a[0] - b[0];
	const float d1 = a[1] - b[1];
	const float d2 = a[2] - b[2];
	return d0 * d0 + d1 * d1 + d2 * d2;
}

array<float, 4> vec4_to_qtn4(const array<float, 3>& axis, const float angle)
{
	//assert(normalized(axis));
	const float h = angle * 0.5f;
	const float s = sin(h);
	const float c = cos(h);
	return
	{
		c,
		s * axis[0],
		s * axis[1],
		s * axis[2],
	};
}

array<float, 4> operator*(const array<float, 4>& a, const array<float, 4>& b)
{
	//assert(normalized(a));
	//assert(normalized(b));
	return
	{
		a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
		a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
		a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
		a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
	};
}

array<float, 9> qtn4_to_mat3(const array<float, 4>& a)
{
	assert(normalized(a));
	const float ww = a[0]*a[0];
	const float wx = a[0]*a[1];
	const float wy = a[0]*a[2];
	const float wz = a[0]*a[3];
	const float xx = a[1]*a[1];
	const float xy = a[1]*a[2];
	const float xz = a[1]*a[3];
	const float yy = a[2]*a[2];
	const float yz = a[2]*a[3];
	const float zz = a[3]*a[3];

	// http://www.boost.org/doc/libs/1_46_1/libs/math/quaternion/TQE.pdf
	// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
	return
	{
		ww+xx-yy-zz, 2*(-wz+xy), 2*(wy+xz),
		2*(wz+xy), ww-xx+yy-zz, 2*(-wx+yz),
		2*(-wy+xz), 2*(wx+yz), ww-xx-yy+zz,
	};
}

array<float, 3> operator*(const array<float, 9>& m, const array<float, 3>& v)
{
	return
	{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
		m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
		m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
	};
}
