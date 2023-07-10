#include "triangle.h"

using namespace std;

struct AABB
{
	Vector3 bounds[2];

	// Default constructor with zero-length bounding box
	AABB()
	{
		bounds[0] = Vector3(0, 0, 0);
		bounds[1] = Vector3(0, 0, 0);
	}

	AABB(Vector3 v0, Vector3 v1)
	{
		bounds[0] = v0;
		bounds[1] = v1;
	}
};

struct BVHNode
{
	AABB bounds;			// Bounding box of the node
	int span[2];			// Index span
	BVHNode* children[2];	// 0 - left; 1 - right

	BVHNode(int from, int to)
	{
		span[0] = from;
		span[1] = to;
		children[0] = children[1] = nullptr;
	}

	~BVHNode()
	{
		delete children[0];
		delete children[1];
		delete children;
		delete span;
	}
};

class BVH
{
public:
	BVH(std::vector<Triangle*>* items);
	~BVH();

	void BuildTree();
	void Traverse(RTCRay& ray);

private:
	BVHNode * BuildTree(int from, int to, int depth);
	void Traverse(RTCRay& ray, BVHNode* node, float t0, float t1);
	//AABB GetNodeAABB(int from, int to);
	bool RayBoxIntersection(RTCRay& ray, AABB &bounds);
	bool RayTriangleIntersection(RTCRay& ray, Triangle &triangle);

	int _maxLeafItems;					// Maximum items in one noed
	BVHNode* _root;						// Root node of the tree
	std::vector<Triangle*>* _items;		// All triangles in tree
};


