#pragma once
#ifndef IDOCK_ARRAY_HPP
#define IDOCK_ARRAY_HPP

#include <array>
using namespace std;

//! Returns the flattened 1D index of a triangular 2D index (x, y) where x is the lowest dimension.
size_t mr(const size_t x, const size_t y);

//! Returns the flattened 1D index of a triangular 2D index (x, y) where either x or y is the lowest dimension.
size_t mp(const size_t x, const size_t y);

//! Returns the square norm of a vector.
float norm_sqr(const array<float, 3>& a);

float norm_sqr(const float a0, const float a1, const float a2);

//! Returns the square norm of a quaternion.
float norm_sqr(const array<float, 4>& a);

float norm_sqr(const float a0, const float a1, const float a2, const float a3);

//! Returns the norm of a vector.
float norm(const array<float, 3>& a);

float norm(const float a0, const float a1, const float a2);

//! Returns the norm of a quaternion.
float norm(const array<float, 4>& a);

float norm(const float a0, const float a1, const float a2, const float a3);

//! Returns true if the norm of a vector is approximately 1.
bool normalized(const array<float, 3>& a);

bool normalized(const float a0, const float a1, const float a2);

//! Returns true if the norm of a quaternion is approximately 1.
bool normalized(const array<float, 4>& a);

bool normalized(const float a0, const float a1, const float a2, const float a3);

//! Normalizes a vector.
array<float, 3> normalize(const array<float, 3>& a);

void normalize(float &a0, float &a1, float &a2);

//! Normalizes a quaternion.
array<float, 4> normalize(const array<float, 4>& a);

void normalize(float &a0, float &a1, float &a2, float &a3);

//! Elementwise adds the second vectors to the first vector.
array<float, 3> operator+(const array<float, 3>& a, const array<float, 3>& b);

//! Elementwise subtracts the second vectors from the first vector.
array<float, 3> operator-(const array<float, 3>& a, const array<float, 3>& b);

//! Elementwise adds the second vectors to the first vector.
void operator+=(array<float, 3>& a, const array<float, 3>& b);

//! Elementwise subtracts the second vectors from the first vector.
void operator-=(array<float, 3>& a, const array<float, 3>& b);

//! Multiplies a scalar to a vector.
array<float, 3> operator*(const float s, const array<float, 3>& a);

//! Returns the cross product of two vectors.
array<float, 3> operator*(const array<float, 3>& a, const array<float, 3>& b);

//! Returns the square Euclidean distance between two vectors.
float distance_sqr(const array<float, 3>& a, const array<float, 3>& b);

//! Constructs a quaternion by a normalized axis and a rotation angle.
array<float, 4> vec4_to_qtn4(const array<float, 3>& axis, const float angle);

//! Returns the product of two quaternions.
array<float, 4> operator*(const array<float, 4>& a, const array<float, 4>& b);

//! Transforms the current quaternion into a 3x3 transformation matrix, e.g. quaternion(1, 0, 0, 0) => identity matrix.
array<float, 9> qtn4_to_mat3(const array<float, 4>& a);

//! Transforms a vector by a 3x3 matrix.
array<float, 3> operator*(const array<float, 9>& m, const array<float, 3>& v);

#endif
