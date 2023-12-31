#include "stdafx.h"
#include "texture.h"
#include "GammaCorrection.h"

Texture::Texture(const char* file_name)
{
	// image format
	FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
	// pointer to the image, once loaded
	FIBITMAP* dib = nullptr;
	// pointer to the image data
	BYTE* bits = nullptr;

	// check the file signature and deduce its format
	fif = FreeImage_GetFileType(file_name, 0);
	// if still unknown, try to guess the file format from the file extension
	if (fif == FIF_UNKNOWN)
	{
		fif = FreeImage_GetFIFFromFilename(file_name);
	}
	// if known
	if (fif != FIF_UNKNOWN)
	{
		// check that the plugin has reading capabilities and load the file
		if (FreeImage_FIFSupportsReading(fif))
		{
			dib = FreeImage_Load(fif, file_name);
		}
		// if the image loaded
		if (dib)
		{
			// retrieve the image data
			bits = FreeImage_GetBits(dib);
			//FreeImage_ConvertToRawBits()
			// get the image width and height
			width_ = int(FreeImage_GetWidth(dib));
			height_ = int(FreeImage_GetHeight(dib));

			// if each of these is ok
			if ((bits != 0) && (width_ != 0) && (height_ != 0))
			{
				// texture loaded
				scan_width_ = FreeImage_GetPitch(dib); // in bytes
				pixel_size_ = FreeImage_GetBPP(dib) / 8; // in bytes

				data_ = new BYTE[scan_width_ * height_]; // BGR(A) format
				FreeImage_ConvertToRawBits(data_, dib, scan_width_, pixel_size_ * 8, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, TRUE);
			}

			FreeImage_Unload(dib);
			bits = nullptr;
		}
	}
}

Texture::~Texture()
{
	if (data_)
	{
		// free FreeImage's copy of the data
		delete[] data_;
		data_ = nullptr;

		width_ = 0;
		height_ = 0;
	}
}

Color3f Texture::get_texel(const int x, const int y) const
{
	//assert( ( x >= 0 && x < width_ ) && ( y >= 0 && y < height_ ) );

	const int offset = y * scan_width_ + x * pixel_size_;

	if (pixel_size_ > 4 * 1) // HDR, EXR
	{
		const float r = ((float*)(data_ + offset))[0];
		const float g = ((float*)(data_ + offset))[1];
		const float b = ((float*)(data_ + offset))[2];

		return Color3f{ c_linear(r), c_linear(g), c_linear(b) };
	}
	else // PNG, JPG etc.
	{
		const float b = data_[offset] / 255.0f;
		const float g = data_[offset + 1] / 255.0f;
		const float r = data_[offset + 2] / 255.0f;

		return Color3f{ c_linear(r), c_linear(g), c_linear(b) };
	}
}

Color3f Texture::get_texel(const float u, const float v) const
{
	//assert( ( u >= 0.0f && u <= 1.0f ) && ( v >= 0.0f && v <= 1.0f ) );	

	// nearest neighbor interpolation
	const int x = max(0, min(width_ - 1, int(u * width_)));
	const int y = max(0, min(height_ - 1, int(v * height_)));

	return get_texel(x, y);
}

int Texture::width() const
{
	return width_;
}

int Texture::height() const
{
	return height_;
}

Color3f Texture::BilinearInterpolation(float x, float y) {
	if (x - 1 < 0 || x + 1 > width() || y - 1 < 0 || y + 1 > height()) {
		return get_texel(x, y);
	}

	Color3f v1 = get_texel(floor(x), floor(y));
	Color3f v2 = get_texel(floor(x), ceil(y));
	Color3f v3 = get_texel(ceil(x), floor(y));
	Color3f v4 = get_texel(ceil(x), ceil(y));

	float tmp;
	float q1 = (1 - modf(x, &tmp)) * (1 - modf(y, &tmp));
	float q2 = (1 - modf(x, &tmp)) * (modf(y, &tmp));
	float q3 = (modf(x, &tmp)) * (1 - modf(y, &tmp));
	float q4 = (modf(x, &tmp)) * (modf(y, &tmp));

	Color3f c1 = Color3f{ v1.r * q1, v1.g * q1, v1.b * q1 };
	Color3f c2 = Color3f{ v2.r * q2, v2.g * q2, v2.b * q2 };
	Color3f c3 = Color3f{ v3.r * q3, v3.g * q3, v3.b * q3 };
	Color3f c4 = Color3f{ v4.r * q4, v4.g * q4, v4.b * q4 };

	Color3f C = Color3f{ c1.r + c2.r + c3.r + c4.r, c1.g + c2.g + c3.g + c4.g, c1.b + c2.b + c3.b + c4.b };

	return C;
}