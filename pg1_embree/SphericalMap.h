#include "objloader.h"
#include "tutorials.h"
#include <string>

class SphericalMap
{
public:
	SphericalMap();
	SphericalMap(const std::string & filename);

	Color3f Texel(const float x, const float y, const float z) const;

private:
	std::unique_ptr<Texture> _texture;
};

