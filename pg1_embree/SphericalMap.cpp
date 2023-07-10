#include "stdafx.h"
#include "SphericalMap.h"
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

SphericalMap::SphericalMap()
{
	this->_texture = nullptr;
}

SphericalMap::SphericalMap(const std::string & filename)
{
	this->_texture = std::make_unique<Texture>(filename.c_str());
}

Color3f SphericalMap::Texel(const float x, const float y, const float z) const
{
	float theta = acosf(z);
	float phi = atan2f(y, x) + (y < 0 ? 2 * M_PI : 0);
	
	float u = 1.0f - (phi * 0.5f * (1/M_PI));
	float v = theta * (1/M_PI);
	
	return _texture->BilinearInterpolation(u, v);
}