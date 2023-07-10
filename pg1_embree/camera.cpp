#include "stdafx.h"
#include "camera.h"
#include <stdlib.h>
#include <time.h>

Camera::Camera( const int width, const int height, const float fov_y,
	const Vector3 view_from, const Vector3 view_at )
{
	width_ = width;
	height_ = height;
	fov_y_ = fov_y;

	view_from_ = view_from;
	view_at_ = view_at;

	// TODO compute focal lenght based on the vertical field of view and the camera resolution
	// f_y_ = ...
	f_y_ = (float)height / (2 * tanf(fov_y / 2));

	// TODO build M_c_w_ matrix	
	// M_c_w_ = Matrix3x3( x_c, y_c, z_c );
	Vector3 z_c = view_from_ - view_at_;
	Vector3 x_c = this->up_.CrossProduct(z_c);
	Vector3 y_c = z_c.CrossProduct(x_c);

	x_c.Normalize();
	y_c.Normalize();
	z_c.Normalize();

	M_c_w_ = Matrix3x3(x_c, y_c, z_c);

	srand(time(NULL));
}

float Camera::GetRandomPositionInRadius(float radius)
{
	srand(time(NULL));
	return (rand() % 2 * radius) - radius;
}

RTCRay Camera::GenerateRay( const float x_i, const float y_i ) const
{
	Vector3 d_c;
	d_c.x = x_i - (float)width_ / 2;
	d_c.y = (float)height_ / 2 - y_i;
	d_c.z = -f_y_;

	d_c.Normalize();

	Vector3 d_w = this->M_c_w_ * d_c;
	//d_w.Normalize();

	// TODO fill in ray structure and compute ray direction
	// ray.org_x = ...	

	// setup a primary ray
	RTCRay ray = RTCRay();
	ray.org_x = this->view_from_.x;
	ray.org_y = this->view_from_.y;
	ray.org_z = this->view_from_.z;
	ray.tnear = 0.01; // start of ray segment

	ray.dir_x = d_w.x; // ray direction
	ray.dir_y = d_w.y;
	ray.dir_z = d_w.z;
	ray.time = 1.0f; // time of this ray for motion blur

	ray.tfar = FLT_MAX; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	// Depth of Field calculations
	
	bool DoF = true;	// Turn Depth of Field ON/OFF;	mainly used for debug, default should be true

	if (DoF)
	{
		Vector3 t_p = Vector3(ray.org_x, ray.org_y, ray.org_z) + Vector3(ray.dir_x, ray.dir_y, ray.dir_z) * f_y_;	// Focal point hit

		float aperture = 0.4f;
		float shift = ((float)(rand()) / ((float)(RAND_MAX / (aperture * 2)))) - aperture;	// Random shift of camera in aperture radius

		ray.org_y += shift;
		ray.org_x += shift;
		//ray.org_x += shift;   maybe uncomment this later

		Vector3 t_v = Vector3(t_p.x, t_p.y, t_p.z) - Vector3(ray.org_x, ray.org_y, ray.org_z);	// Direction vector from new camera origin when working with depth of field
		t_v.Normalize();

		ray.dir_x = t_v.x;
		ray.dir_y = t_v.y;
		ray.dir_z = t_v.z;
	}
	

	return ray;
}
