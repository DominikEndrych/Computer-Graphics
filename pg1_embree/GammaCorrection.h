#include <cmath>
#include <cassert>

static float c_linear(float c_srgb, float gamma = 2.4f)
{
	if (c_srgb <= 0.0f) return 0.0f;
	else if (c_srgb >= 1.0f) return 1.0f;
	assert((c_srgb >= 0.0f) && (c_srgb <= 1.0f));
	if (c_srgb <= 0.04045f)
	{
		return c_srgb / 12.92f;
	}
	else
	{
		const float a = 0.055f;
		return powf((c_srgb + a) / (1.0f + a), gamma);
	}
}


static float c_srgb(float c_linear, float gamma = 2.4f)
{
	if (c_linear <= 0.0f) return 0.0f;
	else if (c_linear >= 1.0f) return 1.0f;
	assert((c_linear >= 0.0f) && (c_linear <= 1.0f));
	if (c_linear <= 0.0031308f)
	{
		return 12.92f * c_linear;
	}
	else
	{
		const float a = 0.055f;
		return (1.0f + a) * powf(c_linear, 1.0f / gamma) - a;
	}
}