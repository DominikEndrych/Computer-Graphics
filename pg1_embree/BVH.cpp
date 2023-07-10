#include "stdafx.h"
#include "BVH.h"

// Constructor
BVH::BVH(std::vector<Triangle*>* items)
{
	_items = items;
	_maxLeafItems = 4;
	_root = nullptr;
}

void BVH::BuildTree()
{
	_root = BuildTree(0, _items->size(), 0);
}

BVHNode* BVH::BuildTree(int from, int to, int depth)
{
	BVHNode* node = new BVHNode(from, to);
	//TODO: node->bounds = GetNodeAABB(from, to);

	if (to - from + 1 > _maxLeafItems)
	{
		int split_axis = depth % 3;					// Get axis for sorting
		int pivot = (to - from + 1) / 2 + from;		// Mid index

		// Sorting triangles
		std::vector<Triangle*>::iterator begin = _items->begin();
		std::nth_element(begin + from, begin + pivot, begin + to + 1, TriangleComparator(split_axis));

		node->children[0] = BuildTree(from, pivot - 1, depth + 1);
		node->children[1] = BuildTree(pivot, to, depth + 1);
	}
	return node;
}

bool BVH::RayBoxIntersection(RTCRay& ray, AABB& bounds)
{
	float tx_0 = ((bounds.bounds[0].x - ray.org_x) / ray.dir_x);
	float ty_0 = ((bounds.bounds[0].y - ray.org_y) / ray.dir_y);
	float tz_0 = ((bounds.bounds[0].z - ray.org_z) / ray.dir_z);

	float tx_1 = ((bounds.bounds[1].x - ray.org_x) / ray.dir_x);
	float ty_1 = ((bounds.bounds[1].y - ray.org_y) / ray.dir_y);
	float tz_1 = ((bounds.bounds[1].z - ray.org_z) / ray.dir_z);

	if (tx_0 > tx_1) { swap(tx_0, tx_1); }
	if (ty_0 > ty_1) { swap(ty_0, ty_1); }
	if (tz_0 > tz_1) { swap(tz_0, tz_1); }

	float t_0 = max(tx_0, ty_0, tz_0);
	float t_1 = max(tx_1, ty_1, tz_1);

	if (t_0 < t_1 && 0 < t_1)
	{
		return true;
	}
	else return false;
}

bool BVH::RayTriangleIntersection(RTCRay& ray, Triangle& triangle)
{
	Vector3 n = triangle.vertex(0).normal;
	Vector3 d = Vector3(ray.dir_x, ray.dir_y, ray.dir_z);
	if (d.DotProduct(n) != 0)
	{
		Vector3 a = triangle.vertex(0).position;
		Vector3 b = triangle.vertex(0).position;
		Vector3 c = triangle.vertex(0).position;
		Vector3 origin = Vector3(ray.org_x, ray.org_y, ray.org_z);
		Vector3 direction = Vector3(ray.dir_x, ray.dir_y, ray.dir_z);
		float t_hit = ((triangle.GetCenter() - origin).DotProduct(n) / direction.DotProduct(n));
		Vector3 P = Vector3(ray.org_x, ray.org_y, ray.org_z) + (Vector3(ray.dir_x, ray.dir_y, ray.dir_z) * t_hit);



		// Compute vectors        
		Vector3 v0 = c - a;
		Vector3 v1 = b - a;
		Vector3 v2 = P - a;

		// Compute dot products
		float dot00 = v0.DotProduct(v0);
		float dot01 = v0.DotProduct(v1);
		float dot02 = v0.DotProduct(v2);
		float dot11 = v1.DotProduct(v1);
		float dot12 = v1.DotProduct(v2);

		// Compute barycentric coordinates
		float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
		float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		// Check if point is in triangle
		return (u >= 0) && (v >= 0) && (u + v < 1);
	}
	else return false;
}

void BVH::Traverse(RTCRay& ray)
{
	Traverse(ray, _root, ray.tnear, ray.tfar);
}

void BVH::Traverse(RTCRay& ray, BVHNode* node, float t0, float t1)
{
	// Ray intersected with bounding box
	if (RayBoxIntersection(ray, node->bounds))
	{
		// Node is leaf
		if (node->children[0] == nullptr && node->children[1] == nullptr)
		{
			for (int i = node->span[0]; i <= node->span[1]; i++)
			{
				// TODO: ray and triangle intersection
				// if(RayTriangleIntersection(ray, _items[i])) ...
			}
		}
		else
		{
			Traverse(ray, node->children[0], t0, ray.tfar);
			Traverse(ray, node->children[1], t0, ray.tfar);
		}
	}
	else
	{
		return;		// Ray missed the bounding box
	}
}

// Destructor
BVH::~BVH()
{
	delete _items;
}