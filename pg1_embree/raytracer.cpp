#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include <random>

#include "GammaCorrection.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>

Raytracer::Raytracer( const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config ) : SimpleGuiDX11( width, height )
{
	InitDeviceAndScene( config );

	camera_ = Camera( width, height, fov_y, view_from, view_at );

	lights_.push_back(Omnilight());

	srand((unsigned)time(NULL));

	// Inicialize uniform random generator
	uniformGenerator = std::mt19937(123);
	uni_dist = std::uniform_real_distribution<float>(0.0f, 1.0f);
}

Raytracer::~Raytracer()
{
	ReleaseDeviceAndScene();
}

int Raytracer::InitDeviceAndScene( const char * config )
{
	device_ = rtcNewDevice( config );
	error_handler( nullptr, rtcGetDeviceError( device_ ), "Unable to create a new device.\n" );
	rtcSetDeviceErrorFunction( device_, error_handler, nullptr );

	ssize_t triangle_supported = rtcGetDeviceProperty( device_, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED );

	// create a new scene bound to the specified device
	scene_ = rtcNewScene( device_ );

	return S_OK;
}

int Raytracer::ReleaseDeviceAndScene()
{
	rtcReleaseScene( scene_ );
	rtcReleaseDevice( device_ );

	return S_OK;
}

void Raytracer::LoadScene( const std::string file_name )
{
	const int no_surfaces = LoadOBJ( file_name.c_str(), surfaces_, materials_ );

	// surfaces loop
	for ( auto surface : surfaces_ )
	{
		RTCGeometry mesh = rtcNewGeometry( device_, RTC_GEOMETRY_TYPE_TRIANGLE );

		Vertex3f * vertices = ( Vertex3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof( Vertex3f ), 3 * surface->no_triangles() );

		Triangle3ui * triangles = ( Triangle3ui * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof( Triangle3ui ), surface->no_triangles() );

		rtcSetGeometryUserData( mesh, ( void* )( surface->get_material() ) );

		rtcSetGeometryVertexAttributeCount( mesh, 2 );

		Normal3f * normals = ( Normal3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof( Normal3f ), 3 * surface->no_triangles() );

		Coord2f * tex_coords = ( Coord2f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof( Coord2f ), 3 * surface->no_triangles() );		

		// triangles loop
		for ( int i = 0, k = 0; i < surface->no_triangles(); ++i )
		{
			Triangle & triangle = surface->get_triangle( i );

			// vertices loop
			for ( int j = 0; j < 3; ++j, ++k )
			{
				const Vertex & vertex = triangle.vertex( j );

				vertices[k].x = vertex.position.x;
				vertices[k].y = vertex.position.y;
				vertices[k].z = vertex.position.z;

				normals[k].x = vertex.normal.x;
				normals[k].y = vertex.normal.y;
				normals[k].z = vertex.normal.z;

				tex_coords[k].u = vertex.texture_coords[0].u;
				tex_coords[k].v = vertex.texture_coords[0].v;
			} // end of vertices loop

			triangles[i].v0 = k - 3;
			triangles[i].v1 = k - 2;
			triangles[i].v2 = k - 1;
		} // end of triangles loop

		rtcCommitGeometry( mesh );
		unsigned int geom_id = rtcAttachGeometry( scene_, mesh );
		rtcReleaseGeometry( mesh );
	} // end of surfaces loop

	// Spherical map
	_shpericalMap = new SphericalMap("../../../data/coastScene.jpg");

	rtcCommitScene( scene_ );
}

RTCRay Raytracer::CreateRay(Vector3 origin, Vector3 direction)
{
	RTCRay ray = RTCRay();
	ray.org_x = origin.x;
	ray.org_y = origin.y;
	ray.org_z = origin.z;
	ray.tnear = 0.01f; // start of ray segment

	ray.dir_x = direction.x;
	ray.dir_y = direction.y;
	ray.dir_z = direction.z;
	ray.time = 1.0f; // time of this ray for motion blur

	ray.tfar = FLT_MAX; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}

float Raytracer::clamp(float n, float lower, float upper) {
	return (std::max)(lower, (std::min)(n, upper));
}

float Raytracer::GetAngle(Vector3 v1, Vector3 v2)
{
	return (v1.DotProduct(v2)) / (v1.L2Norm() * v2.L2Norm());
}

float Raytracer::GetRandom_Uniform(float from, float to)
{
	// Uniform distribution generator
	//std::mt19937 generator(123);
	//std::uniform_real_distribution<float> uni_distt(0.0f, 1.0f);

	float ksi = (uni_dist)(uniformGenerator);
	return ((ksi * (to - from)) / 1.0f) + from;		// Convert and return number in range <from, to>
}

Vector3 Orthogonal(Vector3 normal)
{
	return (abs(normal.x) > abs(normal.z)) ? Vector3(normal.y, -normal.x, 0.0f) : Vector3(0.0f, normal.z, -normal.y);
}

Vector3 GetHemisphereSample(Vector3 normal)
{
	float xi1 = (float)rand() / RAND_MAX;
	float xi2 = (float)rand() / RAND_MAX;

	float x = 2 * cos(2 * M_PI * xi1) * sqrt(xi2 * (1 - xi2));
	float y = 2 * sin(2 * M_PI * xi1) * sqrt(xi2 * (1 - xi2));
	float z = 1 - (2 * xi2);

	Vector3 result = Vector3(x, y, z);
	
	result.Normalize();
	if (result.DotProduct(normal) < 0)
	{
		result *= -1;
	}

	return result;
}

Vector3 GetCosineHemisphereSample(Vector3 normal)
{
	float xi1 = (float)rand() / RAND_MAX;
	float xi2 = (float)rand() / RAND_MAX;

	float x = cos(2 * M_PI * xi1) * sqrt(1 - xi2);
	float y = sin(2 * M_PI * xi1) * sqrt(1 - xi2);
	float z = sqrt(xi2);

	Vector3 omega_i = Vector3(x, y, z);		// Sample on hemisphere

	//Transformation to correct coordinates
	Vector3 o2 = Orthogonal(normal);
	Vector3 o1 = o2.CrossProduct(normal);
	Matrix3 T_rs = Matrix3x3(o1, o2, normal);

	T_rs.Transpose();		// Switch rows and columns

	return T_rs * omega_i;
}

Vector3 GetCosineLobeSample(float gamma)
{
	float xi1 = (float)rand() / RAND_MAX;
	float xi2 = (float)rand() / RAND_MAX;

	float x = cos(2 * M_PI * xi1) * sqrt(1 - pow(xi2, 2 / gamma + 1));
	float y = sin(2 * M_PI * xi1) * sqrt(1 - pow(xi2, 2 / gamma + 1));
	float z = pow(xi2, 1 / gamma + 1);

	return Vector3(x, y, z);
}

Color4f Raytracer::get_pixel( const int x, const int y, const float t )
{
	int pixelDimension = 1;							// Make AxA matrix from pixel
	float cellSize = 1.0f / (float)pixelDimension;	// Size of one cell

	Vector3 colorSum = Vector3(0, 0, 0);

	/*
	// Uniform distribution generator
	std::mt19937 generator(123);
	std::uniform_real_distribution<float> uni_dist(0.0f, 1.0f);
	*/



	/*
	for (int i = 0; i < pixelDimension; i++)
	{
		for (int j = 0; j < pixelDimension; j++)
		{
			float p_x = GetRandom_Uniform(x + (i * cellSize), x + (i * cellSize) + cellSize, &generator, &uni_dist);
			float p_y = GetRandom_Uniform(y + (j * cellSize), y + (j * cellSize) + cellSize, &generator, &uni_dist);

			colorSum += TraceRay(camera_.GenerateRay(p_x, p_y), 10);
		}
	}
	*/

	//Vector3 color = TraceRay(camera_.GenerateRay(x, y), 10);
	
	//Vector3 color = Vector3(colorSum.x / cells, colorSum.y / cells, colorSum.z / cells);
	
	/*
	int samples = 800;
	for (int i = 0; i < samples; i++)
	{
		colorSum += TraceRay(camera_.GenerateRay(x, y), 20);
	}
	Vector3 color = Vector3(colorSum.x / samples, colorSum.y / samples, colorSum.z / samples);D
	*/

	float cells = (float)pixelDimension * (float)pixelDimension;

	int samples = 1;
	for (int i = 0; i < samples; i++)
	{
		// For each pixel:
		Vector3 pixelSum = Vector3(0.f, 0.f, 0.f);
		for (int i = 0; i < pixelDimension; i++)
		{
			for (int j = 0; j < pixelDimension; j++)
			{
				float p_x = GetRandom_Uniform(x + (i * cellSize), x + (i * cellSize) + cellSize);
				float p_y = GetRandom_Uniform(y + (j * cellSize), y + (j * cellSize) + cellSize);

				pixelSum += TraceRay(camera_.GenerateRay(p_x, p_y), 10);
			}
		}
		pixelSum = pixelSum / cells;	// Pixel after supersampling

		colorSum += pixelSum;
	}

	Vector3 color = Vector3(colorSum.x / samples, colorSum.y / samples, colorSum.z / samples);

	/*
	Ns = lesklost (gama)
	ks = spec reflexivity (ms)
	kd = diff refl. (md)
	*/

	return Color4f{ c_srgb(color.x), c_srgb(color.y), c_srgb(color.z), c_srgb(1.0f) };
}

bool Raytracer::ShadowRay(Vector3 origin, Vector3 direction, Vector3 lightPosition)
{
	RTCRay shadowRay = CreateRay(origin, direction);	// Create shadow ray

	// setup a hit
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f; // geometry normal
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;

	// merge ray and hit structures
	RTCRayHit ray_hit;
	ray_hit.ray = shadowRay;
	ray_hit.hit = hit;

	// intersect ray with the scene
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(this->scene_, &context, &ray_hit);

	Vector3 v_light = lightPosition - origin;	// Vector from origin to light
	float lightDistance = v_light.L2Norm();		// Length of segment from origin to light

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		// Somathing was hit
		if (ray_hit.ray.tfar < lightDistance)	// If hit was behind the light source, do not cast shadow
		{
			return true;
		}
	}

	// Nothing was hit
	return false;
}

Vector3 Raytracer::TraceRay(RTCRay ray, int currentDepth, int maxDepth)
{
	if (currentDepth <= 0)
	{
		return Vector3{ 0.f, 0.f, 0.f };
	}

	// setup a hit
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f; // geometry normal
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;

	// merge ray and hit structures
	RTCRayHit ray_hit;
	ray_hit.ray = ray;
	ray_hit.hit = hit;

	// intersect ray with the scene
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(this->scene_, &context, &ray_hit);

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		// we hit something
		RTCGeometry geometry = rtcGetGeometry(this->scene_, ray_hit.hit.geomID);
		Normal3f normal;
		// get interpolated normal
		rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);
		// and texture coordinates
		Coord2f tex_coord;
		rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v,
			RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord.u, 2);

		Material* material = (Material*)(rtcGetGeometryUserData(geometry));	// Material

		Vector3 d = Vector3(ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z);
		Vector3 o = Vector3(ray_hit.ray.org_x, ray_hit.ray.org_y, ray_hit.ray.org_z);
		//d.Normalize();
		Vector3 p = o + d * ray_hit.ray.tfar;	// hit point

		Vector3 l = lights_[0].position - p;	// l vector
		
		Vector3 n = Vector3(normal.x, normal.y, normal.z);
		l.Normalize();
		n.Normalize();
		Vector3 l_r = 2 * (n.DotProduct(l)) * n - l;	// Reflection from light to camera
		l_r.Normalize();

		Vector3 v = d * -1;
		//v.Normalize();

		// Flip normal
		if (n.DotProduct(v) < 0)
		{
			n *= -1;
		}

		Vector3 resultColor = Vector3(0.f, 0.f, 0.f);	// Color to return

		// Normal 
		if (material->shader == 0)
		{
			Vector3 converted = (n + Vector3(1.f, 1.f, 1.f)) / 2;
			resultColor = converted;
		}
		// Diffuse
		if(material->shader == 1)
		{
			resultColor = lights_[0].intensity * material->diffuse * (clamp(n.DotProduct(l), 0, 1));
		}
		// Phong
		else if (material->shader == 2)
		{

			Vector3 diffuse = material->diffuse;	// Diffuse part

			// If material has some kind of texture
			Texture* tex_diffuse = material->get_texture(Material::kDiffuseMapSlot);
			if (tex_diffuse)
			{
				Color3f diffuseColor = tex_diffuse->BilinearInterpolation(tex_coord.u, 1.0f - tex_coord.v);
				diffuse = Vector3(diffuseColor.r, diffuseColor.g, diffuseColor.b);				// Change diffuse color to color from texture
			}

			bool isShadowed = ShadowRay(p, l, lights_[0].position);		// Is point shadowed

			if (!isShadowed)
			{
				resultColor = lights_[0].intensity * material->ambient +
					lights_[0].intensity * diffuse * (clamp(n.DotProduct(l), 0, 1)) +
					lights_[0].intensity * material->specular * pow((v.DotProduct(l_r)), material->shininess);
			}
			else 
			{
				// This might not be correct
				resultColor = material->diffuse * (clamp(n.DotProduct(l), 0, 1));
				//resultColor = Vector3(0.f, 0.f, 0.f);
			}
			
		}
		// Whitted
		else if (material->shader == 3)
		{
			Vector3 r = (2 * (v.DotProduct(n))) * n - v;
			r.Normalize();

			RTCRay refraction = CreateRay(p, r);

			resultColor = lights_[0].intensity * material->ambient +
				lights_[0].intensity * material->specular * TraceRay(refraction, currentDepth-1);
		}
		// Dielectric materials
		else if (material->shader == 4)
		{
			// Reflection ray
			Vector3 r = (2 * (v.DotProduct(n))) * n - v;	// Reflection direction
			RTCRay reflection = CreateRay(p, r);			// Reflection ray

			// Refraction ray
			float refIndex_1, refIndex_2;
			refIndex_1 = ray_hit.ray.time;
			if (refIndex_1 == 1.0f) { refIndex_2 = 1.4f; }	// From air to material
			else { refIndex_2 = 1.0f; }						// From material to air

			/*TEST DATA FROM PRESENTATION*/
			/*
			refIndex_1 = 1.5f;
			refIndex_2 = 1.0f;
			d = Vector3(-0.429, 0.903, 0);
			v = -d;
			n = Vector3(0, 1, 0);
			*/

			float underSqrt = 1 - pow((refIndex_1 / refIndex_2), 2) * (1 - pow((d.DotProduct(n)), 2));		// cos of refraction angle

			if(underSqrt < 0)
			{
				resultColor = (TraceRay(reflection, currentDepth-1));
				
				if (refIndex_1 != 1)
				{
					resultColor.x *= exp((-material->diffuse.x) * ray_hit.ray.tfar);
					resultColor.y *= exp((-material->diffuse.y) * ray_hit.ray.tfar);
					resultColor.z *= exp((-material->diffuse.z) * ray_hit.ray.tfar);
				}
				
			}
			else
			{
				Vector3 refr_d = (refIndex_1 / refIndex_2) * d -
					((refIndex_1 / refIndex_2) * (d.DotProduct(n)) + sqrt(underSqrt)) * n;	// Direction of refraction ray

				//refr_d.Normalize();
				RTCRay refraction = CreateRay(p, refr_d);	// Refraction ray
				refraction.time = refIndex_2;				// Current refractive index in material

				// Fresnel equations
				float R_s = pow(((refIndex_2 * sqrt(underSqrt) - refIndex_1 * (n.DotProduct(v))) / (refIndex_2 * sqrt(underSqrt) + refIndex_1 * (n.DotProduct(v)))), 2);
				float R_p = pow(((refIndex_2 * (n.DotProduct(v)) - refIndex_1 * sqrt(underSqrt)) / (refIndex_2 * (n.DotProduct(v)) + refIndex_1 * sqrt(underSqrt))), 2);
				float R = (R_s + R_p) / 2;
				float T = 1 - R;

				resultColor = (TraceRay(reflection, currentDepth - 1) * R) + (TraceRay(refraction, currentDepth - 1) * T);

				
				if (refIndex_1 != 1)
				{
					resultColor.x *= exp((-material->diffuse.x) * ray_hit.ray.tfar);
					resultColor.y *= exp((-material->diffuse.y) * ray_hit.ray.tfar);
					resultColor.z *= exp((-material->diffuse.z) * ray_hit.ray.tfar);
				}
				
			}
		}
		// Lambert
		else if (material->shader == 5)
		{
			// Russian Roulette
			float russianAlpha = material->diffuse.data[material->diffuse.LargestComponent()];
			float russianRandom = (float)rand() / RAND_MAX;
			if (russianAlpha <= russianRandom) { return Vector3(0.f, 0.f, 0.f); }

			Vector3 L_e = material->emission;
			if (L_e.x != 0.0f || L_e.y != 0.0f || L_e.z != 0.0f) 
			{ 
				return L_e;		// Light source was hit and recursion can stop here
			}	
			
			//Vector3 omega_i = GetHemisphereSample(n);
			Vector3 omega_i = GetCosineHemisphereSample(n);
			float theta = omega_i.DotProduct(n);

			clamp(theta, 0, 1);

			//float pdf = 1 / (2 * M_PI);
			float pdf = (cos(theta)) / M_PI;

			Vector3 L_i = TraceRay(CreateRay(p, omega_i), currentDepth-1);

			Vector3 f_r = material->diffuse / M_PI;
			Vector3 L_r = (L_i * f_r * omega_i.DotProduct(n)) / (pdf * russianAlpha);

			return L_r;
		}
		// Mirror (global)
		else if (material->shader == 6)
		{
			Vector3 r = (2 * (v.DotProduct(n))) * n - v;
			r.Normalize();

			RTCRay reflection = CreateRay(p, r);

			resultColor = material->specular * TraceRay(reflection, currentDepth - 1);
		}
		// Phong global
		else if (material->shader == 7)
		{
			
			// Russian Roulette
			float russianAlpha = material->diffuse.data[material->diffuse.LargestComponent()];
			float russianRandom = (float)rand() / RAND_MAX;
			if (russianAlpha <= russianRandom) { return Vector3(0.f, 0.f, 0.f); }

			Vector3 L_e = material->emission;
			if (L_e.x != 0.0f || L_e.y != 0.0f || L_e.z != 0.0f)
			{
				return L_e;		// Light source was hit and recursion can stop here
			}

			float diffuseMax = material->diffuse.data[material->diffuse.LargestComponent()];	// Max value in diffuse component
			float specularMax = material->specular.data[material->specular.LargestComponent()];	// Max value in diffuse component

			float randMax = diffuseMax + specularMax;
			float ksi = (float)(rand()) / ((float)(RAND_MAX / randMax));	// Random number in range 0 - randMax

			Vector3 omega_i;
			float pdf = 1.f;
			Vector3 L_i;
			Vector3 f_r;

			if (ksi < diffuseMax)
			{
				// Use diffuse
				omega_i = GetCosineHemisphereSample(n);

				float theta = omega_i.DotProduct(n);
				clamp(theta, 0, 1);

				float pdf = (cos(theta)) / M_PI;
				pdf *= diffuseMax / (diffuseMax + specularMax);

				L_i = TraceRay(CreateRay(p, omega_i), currentDepth - 1);
				f_r = material->diffuse / M_PI;
			}
			else
			{
				// Use specular
				omega_i = GetCosineLobeSample(material->shininess);
				pdf = (material->shininess + 1 / 2 * M_PI) * pow(omega_i.z, material->shininess);
				pdf *= specularMax / (specularMax + diffuseMax);

				L_i = TraceRay(CreateRay(p, omega_i), currentDepth - 1);
				f_r = material->specular * ((material->shininess) + 2 / 2 * M_PI);	// Might not be correct
			}

			Vector3 L_r = (L_i * f_r * omega_i.DotProduct(n)) / (pdf * russianAlpha);

			return L_r;
		}
		
		return resultColor;
	}

	// Nothing was hit at this point
	
	// Return background color from map
	Color3f mColor = _shpericalMap->Texel(ray_hit.ray.dir_x, ray_hit.ray.dir_y, ray_hit.ray.dir_z);

	return Vector3(mColor.r, mColor.g, mColor.b);
	//return Vector3(0, 0, 0);	// Return black background

	//return Vector3( 1.0f, 1.0f, 1.0f );
}

int Raytracer::Ui()
{
	static float f = 0.0f;
	static int counter = 0;

	// Use a Begin/End pair to created a named window
	ImGui::Begin( "Ray Tracer Params" );
	
	ImGui::Text( "Surfaces = %d", surfaces_.size() );
	ImGui::Text( "Materials = %d", materials_.size() );
	ImGui::Separator();
	ImGui::Checkbox( "Vsync", &vsync_ );
	
	//ImGui::Checkbox( "Demo Window", &show_demo_window ); // Edit bools storing our window open/close state
	//ImGui::Checkbox( "Another Window", &show_another_window );

	ImGui::SliderFloat( "float", &f, 0.0f, 1.0f ); // Edit 1 float using a slider from 0.0f to 1.0f    
	//ImGui::ColorEdit3( "clear color", ( float* )&clear_color ); // Edit 3 floats representing a color

	// Buttons return true when clicked (most widgets return true when edited/activated)
	if ( ImGui::Button( "Button" ) )
		counter++;
	ImGui::SameLine();
	ImGui::Text( "counter = %d", counter );

	ImGui::Text( "Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate );
	ImGui::End();

	// 3. Show another simple window.
	/*if ( show_another_window )
	{
	ImGui::Begin( "Another Window", &show_another_window ); // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
	ImGui::Text( "Hello from another window!" );
	if ( ImGui::Button( "Close Me" ) )
	show_another_window = false;
	ImGui::End();
	}*/

	return 0;
}
