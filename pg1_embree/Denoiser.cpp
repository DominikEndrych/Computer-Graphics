#include "stdafx.h"
#include "Denoiser.h"
#include "OpenImageDenoise/oidn.h"

void convert(Texture & texture3u, std::vector<float> & texture3f, const bool is_normal = false)
{
	for (int y = 0; y < texture3u.height(); ++y)
	{
		for (int x = 0; x < texture3u.width(); ++x)
		{
			for (int c = 0; c < 3; ++c)
			{
				const size_t offset = (y * texture3u.width() + x) * 3 + c;
				float value = texture3u.data_[offset] / 255.0f;

				if (is_normal)
				{
					value = 2.0f * value - 1.0f;
				}

				texture3f.push_back(value);
			}
		}
	}
}

Denoiser::Denoiser()
{
	return;
}

Denoiser::~Denoiser()
{
	return;
}

void Denoiser::Run()
{
	Texture texture("../../../data/test4.png");
	std::vector<float> texture3f;
	convert(texture, texture3f);

	const float * data = texture3f.data();
}
