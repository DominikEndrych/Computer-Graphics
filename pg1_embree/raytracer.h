#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "Omnilight.h"
#include "SphericalMap.h"
#include <random>

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/
class Raytracer : public SimpleGuiDX11
{
public:
	Raytracer( const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3" );
	~Raytracer();

	int InitDeviceAndScene( const char * config );

	int ReleaseDeviceAndScene();

	void LoadScene( const std::string file_name );

	Color4f get_pixel( const int x, const int y, const float t = 0.0f ) override;

	Vector3 TraceRay(RTCRay ray, int currentDepth, int maxDepth = 1000);
	RTCRay CreateRay(Vector3 origin, Vector3 direction);							// Create ray
	bool ShadowRay(Vector3 origin, Vector3 direction, Vector3 lightPosition);		// Function for hard shadows
	float GetAngle(Vector3 v1, Vector3 v2);
	float clamp(float n, float lower, float upper);

	int Ui();

private:
	float GetRandom_Uniform(float from, float to);

	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;
	std::vector<Omnilight> lights_;

	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;

	SphericalMap *_shpericalMap;

	// Uniform random generator variables
	std::mt19937 uniformGenerator;
	std::uniform_real_distribution<float> uni_dist;
};
