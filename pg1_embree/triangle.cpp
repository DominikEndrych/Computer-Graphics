#include "stdafx.h"
#include "triangle.h"

Triangle::Triangle( const Vertex & v0, const Vertex & v1, const Vertex & v2, Surface * surface )
{
	vertices_[0] = v0;
	vertices_[1] = v1;
	vertices_[2] = v2;	

	assert( !is_degenerate() );
}

Vertex Triangle::vertex( const int i )
{
	return vertices_[i];
}

bool Triangle::is_degenerate() const
{
	return vertices_[0].position == vertices_[1].position ||
		vertices_[0].position == vertices_[2].position || 
		vertices_[1].position == vertices_[2].position;
}

Vector3 Triangle::GetCenter()
{
	Vector3 a, b, c;
	a = this->vertex(0).position;
	b= this->vertex(0).position;
	c= this->vertex(0).position;

	float x = (a.x + b.x + c.x) / 3;
	float y = (a.y + b.y + c.y) / 3;
	float z = (a.z + b.z + c.z) / 3;

	return Vector3(x, y, z);
}


