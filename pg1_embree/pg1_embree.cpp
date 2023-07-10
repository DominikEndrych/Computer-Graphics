#include "stdafx.h"
#include "tutorials.h"
#include "Denoiser.h"

int main()
{
	printf( "PG1, (c)2011-2020 Tomas Fabian\n\n" );

	_MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
	_MM_SET_DENORMALS_ZERO_MODE( _MM_DENORMALS_ZERO_ON );

	//return tutorial_1();
	//return tutorial_2();
	//return tutorial_3( "../../../data/6887_allied_avenger.obj" );		// Ship
	//return tutorial_3("../../../data/geosphere.obj");					// Sphere

	// Last time I worked here in PG1
	//return tutorial_3("../../../data/cornell_box2.obj");				// Cornell box
	
	Denoiser denoiser = Denoiser();

	denoiser.Run();

}
