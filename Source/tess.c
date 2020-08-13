/*
** SGI FREE SOFTWARE LICENSE B (Version 2.0, Sept. 18, 2008)
** Copyright (C) [dates of first publication] Silicon Graphics, Inc.
** All Rights Reserved.
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
** of the Software, and to permit persons to whom the Software is furnished to do so,
** subject to the following conditions:
**
** The above copyright notice including the dates of first publication and either this
** permission notice or a reference to http://oss.sgi.com/projects/FreeB/ shall be
** included in all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
** INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
** PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL SILICON GRAPHICS, INC.
** BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
** TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
** OR OTHER DEALINGS IN THE SOFTWARE.
**
** Except as contained in this notice, the name of Silicon Graphics, Inc. shall not
** be used in advertising or otherwise to promote the sale, use or other dealings in
** this Software without prior written authorization from Silicon Graphics, Inc.
*/
/*
** Author: Eric Veach, July 1994.
*/

#include <stddef.h>
#include <assert.h>
#include <setjmp.h>
#include "bucketalloc.h"
#include "tess.h"
#include "mesh.h"
#include "sweep.h"
#include "geom.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

#define Dot(u,v)	(u[0]*v[0] + u[1]*v[1] + u[2]*v[2])

#if defined(FOR_TRITE_TEST_PROGRAM) || defined(TRUE_PROJECT)
static void Normalize( TESSreal v[3] )
{
	TESSreal len = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

	assert( len > 0 );
	len = sqrtf( len );
	v[0] /= len;
	v[1] /= len;
	v[2] /= len;
}
#endif

#define ABS(x)	((x) < 0 ? -(x) : (x))

static int LongAxis( TESSreal v[3] )
{
	int i = 0;

	if( ABS(v[1]) > ABS(v[0]) ) { i = 1; }
	if( ABS(v[2]) > ABS(v[i]) ) { i = 2; }
	return i;
}

static int ShortAxis( TESSreal v[3] )
{
	int i = 0;

	if( ABS(v[1]) < ABS(v[0]) ) { i = 1; }
	if( ABS(v[2]) < ABS(v[i]) ) { i = 2; }
	return i;
}

#if 0

#ifdef FOR_TRITE_TEST_PROGRAM
#include <stdlib.h>
extern int RandomSweep;
#define S_UNIT_X	(RandomSweep ? (2*drand48()-1) : 1.0)
#define S_UNIT_Y	(RandomSweep ? (2*drand48()-1) : 0.0)
#else
#if defined(SLANTED_SWEEP)
/* The "feature merging" is not intended to be complete.  There are
* special cases where edges are nearly parallel to the sweep line
* which are not implemented.  The algorithm should still behave
* robustly (ie. produce a reasonable tesselation) in the presence
* of such edges, however it may miss features which could have been
* merged.  We could minimize this effect by choosing the sweep line
* direction to be something unusual (ie. not parallel to one of the
* coordinate axes).
*/
#define S_UNIT_X	(TESSreal)0.50941539564955385	/* Pre-normalized */
#define S_UNIT_Y	(TESSreal)0.86052074622010633
#else
#define S_UNIT_X	(TESSreal)1.0
#define S_UNIT_Y	(TESSreal)0.0
#endif
#endif

/* Determine the polygon normal and project vertices onto the plane
* of the polygon.
*/
void tessProjectPolygon( TESStesselator *tess )
{
	TESSvertex *vHead = &tess->mesh->vHead;
	TESSreal norm[3];
	int computedNormal = FALSE;

	norm[0] = tess->normal[0];
	norm[1] = tess->normal[1];
	norm[2] = tess->normal[2];
	if( norm[0] == 0 && norm[1] == 0 && norm[2] == 0 ) {
		ComputeNormal( tess, norm );
		computedNormal = TRUE;
	}

	TESSreal* sUnit = tess->sUnit;
	TESSreal* tUnit = tess->tUnit;
	const int i = LongAxis( norm );

#if defined(FOR_TRITE_TEST_PROGRAM) || defined(TRUE_PROJECT)
	/* Choose the initial sUnit vector to be approximately perpendicular
	* to the normal.
	*/
	Normalize( norm );

	sUnit[i] = 0;
	sUnit[(i+1)%3] = S_UNIT_X;
	sUnit[(i+2)%3] = S_UNIT_Y;

	/* Now make it exactly perpendicular */
	w = Dot( sUnit, norm );
	sUnit[0] -= w * norm[0];
	sUnit[1] -= w * norm[1];
	sUnit[2] -= w * norm[2];
	Normalize( sUnit );

	/* Choose tUnit so that (sUnit,tUnit,norm) form a right-handed frame */
	tUnit[0] = norm[1]*sUnit[2] - norm[2]*sUnit[1];
	tUnit[1] = norm[2]*sUnit[0] - norm[0]*sUnit[2];
	tUnit[2] = norm[0]*sUnit[1] - norm[1]*sUnit[0];
	Normalize( tUnit );
#else
	/* Project perpendicular to a coordinate axis -- better numerically */
	sUnit[i] = 0;
	sUnit[(i+1)%3] = S_UNIT_X;
	sUnit[(i+2)%3] = S_UNIT_Y;

	tUnit[i] = 0;
	tUnit[(i+1)%3] = (norm[i] > 0) ? -S_UNIT_Y : S_UNIT_Y;
	tUnit[(i+2)%3] = (norm[i] > 0) ? S_UNIT_X : -S_UNIT_X;
#endif

	/* Project the vertices onto the sweep plane */
	for(TESSvertex* v = vHead->next; v != vHead; v = v->next )
	{
		v->s = Dot( v->coords, sUnit );
		v->t = Dot( v->coords, tUnit );
	}
	if( computedNormal ) {
		CheckOrientation( tess );
	}

	/* Compute ST bounds. */
	int first = 1;
	for(TESSvertex* v = vHead->next; v != vHead; v = v->next )
	{
		if (first)
		{
			tess->bmin[0] = tess->bmax[0] = v->s;
			tess->bmin[1] = tess->bmax[1] = v->t;
			first = 0;
		}
		else
		{
			if (v->s < tess->bmin[0]) tess->bmin[0] = v->s;
			if (v->s > tess->bmax[0]) tess->bmax[0] = v->s;
			if (v->t < tess->bmin[1]) tess->bmin[1] = v->t;
			if (v->t > tess->bmax[1]) tess->bmax[1] = v->t;
		}
	}
}

#endif

void jerryPrepareVertices(TESStesselator* tess)
{
	TESSvertex* vHead = &tess->mesh->vHead;
	
	const float flip = tess->normal[2] > 0 ? 1.0f : -1.0f;
	
	// Compute ST bounds.
	tess->bmin[0] = 1e30f;
	tess->bmax[0] = -1e30f;
	tess->bmin[1] = 1e30f;
	tess->bmax[1] = -1e30f;
	for (TESSvertex* v = vHead->next; v != vHead; v = v->next)
	{
		v->t = flip * v->t; // Flip for SimView

		if (v->s < tess->bmin[0]) tess->bmin[0] = v->s;
		if (v->s > tess->bmax[0]) tess->bmax[0] = v->s;
		if (v->t < tess->bmin[1]) tess->bmin[1] = v->t;
		if (v->t > tess->bmax[1]) tess->bmax[1] = v->t;
	}
}

#define AddWinding(eDst,eSrc)	(eDst->winding += eSrc->winding, \
	eDst->Sym->winding += eSrc->Sym->winding)

/* tessMeshTessellateMonoRegion( face ) tessellates a monotone region
* (what else would it do??)  The region must consist of a single
* loop of half-edges (see mesh.h) oriented CCW.  "Monotone" in this
* case means that any vertical line intersects the interior of the
* region in a single interval.
*
* Tessellation consists of adding interior edges (actually pairs of
* half-edges), to split the region into non-overlapping triangles.
*
* The basic idea is explained in Preparata and Shamos (which I don''t
* have handy right now), although their implementation is more
* complicated than this one.  The are two edge chains, an upper chain
* and a lower chain.  We process all vertices from both chains in order,
* from right to left.
*
* The algorithm ensures that the following invariant holds after each
* vertex is processed: the untessellated region consists of two
* chains, where one chain (say the upper) is a single edge, and
* the other chain is concave.  The left vertex of the single edge
* is always to the left of all vertices in the concave chain.
*
* Each step consists of adding the rightmost unprocessed vertex to one
* of the two chains, and forming a fan of triangles from the rightmost
* of two chain endpoints.  Determining whether we can add each triangle
* to the fan is a simple orientation test.  By making the fan as large
* as possible, we restore the invariant (check it yourself).
*/
void tessMeshTessellateMonoRegion( TESSmesh *mesh, TESSface *face )
{
	TESShalfEdge *up, *lo;

	/* All edges are oriented CCW around the boundary of the region.
	* First, find the half-edge whose origin vertex is rightmost.
	* Since the sweep goes from left to right, face->anEdge should
	* be close to the edge we want.
	*/
	up = face->anEdge;
	assert( up->Lnext != up && up->Lnext->Lnext != up );

	for( ; VertLeq( up->Dst, up->Org ); up = up->Lprev )
		;
	for( ; VertLeq( up->Org, up->Dst ); up = up->Lnext )
		;
	lo = up->Lprev;

	while( up->Lnext != lo ) {
		if( VertLeq( up->Dst, lo->Org )) {
			/* up->Dst is on the left.  It is safe to form triangles from lo->Org.
			* The EdgeGoesLeft test guarantees progress even when some triangles
			* are CW, given that the upper and lower chains are truly monotone.
			*/
			while( lo->Lnext != up && (EdgeGoesLeft( lo->Lnext )
				|| EdgeSign( lo->Org, lo->Dst, lo->Lnext->Dst ) <= 0 )) {
					TESShalfEdge *tempHalfEdge= tessMeshConnect( mesh, lo->Lnext, lo );
					lo = tempHalfEdge->Sym;
			}
			lo = lo->Lprev;
		} else {
			/* lo->Org is on the left.  We can make CCW triangles from up->Dst. */
			while( lo->Lnext != up && (EdgeGoesRight( up->Lprev )
				|| EdgeSign( up->Dst, up->Org, up->Lprev->Org ) >= 0 )) {
					TESShalfEdge *tempHalfEdge= tessMeshConnect( mesh, up, up->Lprev );
					up = tempHalfEdge->Sym;
			}
			up = up->Lnext;
		}
	}

	/* Now lo->Org == up->Dst == the leftmost vertex.  The remaining region
	* can be tessellated in a fan from this leftmost vertex.
	*/
	assert( lo->Lnext != up );
	while( lo->Lnext->Lnext != up ) {
		TESShalfEdge *tempHalfEdge= tessMeshConnect( mesh, lo->Lnext, lo );
		lo = tempHalfEdge->Sym;
	}
}

/* tessMeshTessellateInterior( mesh ) tessellates each region of
* the mesh which is marked "inside" the polygon.  Each such region
* must be monotone.
*/
void tessMeshTessellateInterior( TESSmesh *mesh )
{
	TESSface *f, *next;

	/*LINTED*/
	for( f = mesh->fHead.next; f != &mesh->fHead; f = next ) {
		/* Make sure we don''t try to tessellate the new triangles. */
		next = f->next;
		if( f->inside ) 
		{
			tessMeshTessellateMonoRegion( mesh, f );
		}
	}
}


typedef struct EdgeStackNode EdgeStackNode;
typedef struct EdgeStack EdgeStack;

struct EdgeStackNode {
	TESShalfEdge *edge;
	EdgeStackNode *next;
};

struct EdgeStack {
	EdgeStackNode *top;
	struct BucketAlloc *nodeBucket;
};

int stackInit( EdgeStack *stack, TESSalloc *alloc )
{
	stack->top = NULL;
	stack->nodeBucket = createBucketAlloc( alloc, "CDT nodes", sizeof(EdgeStackNode), 512 );
	return stack->nodeBucket != NULL;
}

void stackDelete( EdgeStack *stack )
{
    deleteBucketAlloc( stack->nodeBucket );
}

int stackEmpty( EdgeStack *stack )
{
	return stack->top == NULL;
}

void stackPush( EdgeStack *stack, TESShalfEdge *e )
{
	EdgeStackNode *node = (EdgeStackNode *)bucketAlloc( stack->nodeBucket );
	node->edge = e;
	node->next = stack->top;
	stack->top = node;
}

TESShalfEdge *stackPop( EdgeStack *stack )
{
	TESShalfEdge *e = NULL;
	EdgeStackNode *node = stack->top;
	if (node) {
		stack->top = node->next;
		e = node->edge;
		bucketFree( stack->nodeBucket, node );
	}
	return e;
}


//	Starting with a valid triangulation, uses the Edge Flip algorithm to
//	refine the triangulation into a Constrained Delaunay Triangulation.
void tessMeshRefineDelaunay( TESSmesh *mesh, TESSalloc *alloc )
{
	// At this point, we have a valid, but not optimal, triangulation.
	// We refine the triangulation using the Edge Flip algorithm
	//
	//  1) Find all internal edges
	//	2) Mark all dual edges
	//	3) insert all dual edges into a queue

	TESSface *f;
	EdgeStack stack;
	TESShalfEdge *e;
	int maxFaces = 0, maxIter = 0, iter = 0;

	stackInit(&stack, alloc);

	for( f = mesh->fHead.next; f != &mesh->fHead; f = f->next ) {
		if ( f->inside) {
			e = f->anEdge;
			do {
				e->mark = EdgeIsInternal(e); // Mark internal edges
				if (e->mark && !e->Sym->mark) stackPush(&stack, e); // Insert into queue
				e = e->Lnext;
			} while (e != f->anEdge);
			maxFaces++;
		}
	}

	// The algorithm should converge on O(n^2), since the predicate is not robust,
	// we'll save guard against infinite loop.
	maxIter = maxFaces * maxFaces;

	// Pop stack until we find a reversed edge
	// Flip the reversed edge, and insert any of the four opposite edges
	// which are internal and not already in the stack (!marked)
	while (!stackEmpty(&stack) && iter < maxIter) {
		e = stackPop(&stack);
		e->mark = e->Sym->mark = 0;
		if (!tesedgeIsLocallyDelaunay(e)) {
			TESShalfEdge *edges[4];
			int i;
			tessMeshFlipEdge(mesh, e);
			// for each opposite edge
			edges[0] = e->Lnext;
			edges[1] = e->Lprev;
			edges[2] = e->Sym->Lnext;
			edges[3] = e->Sym->Lprev;
			for (i = 0; i < 4; i++) {
				if (!edges[i]->mark && EdgeIsInternal(edges[i])) {
					edges[i]->mark = edges[i]->Sym->mark = 1;
					stackPush(&stack, edges[i]);
				}
			}
		}
		iter++;
	}

	stackDelete(&stack);
}


/* tessMeshDiscardExterior( mesh ) zaps (ie. sets to NULL) all faces
* which are not marked "inside" the polygon.  Since further mesh operations
* on NULL faces are not allowed, the main purpose is to clean up the
* mesh so that exterior loops are not represented in the data structure.
*/
void tessMeshDiscardExterior( TESSmesh *mesh )
{
	TESSface *f, *next;

	/*LINTED*/
	for( f = mesh->fHead.next; f != &mesh->fHead; f = next ) {
		/* Since f will be destroyed, save its next pointer. */
		next = f->next;
		if( ! f->inside ) {
			tessMeshZapFace( mesh, f );
		}
	}
}

/* tessMeshSetWindingNumber( mesh, value, keepOnlyBoundary ) resets the
* winding numbers on all edges so that regions marked "inside" the
* polygon have a winding number of "value", and regions outside
* have a winding number of 0.
*
* If keepOnlyBoundary is TRUE, it also deletes all edges which do not
* separate an interior region from an exterior one.
*/
void tessMeshSetWindingNumber( TESSmesh *mesh, int value,
							 int keepOnlyBoundary )
{
	TESShalfEdge *e, *eNext;

	for( e = mesh->eHead.next; e != &mesh->eHead; e = eNext ) {
		eNext = e->next;
		if( e->Rface->inside != e->Lface->inside ) {

			/* This is a boundary edge (one side is interior, one is exterior). */
			e->winding = (e->Lface->inside) ? value : -value;
		} else {

			/* Both regions are interior, or both are exterior. */
			if( ! keepOnlyBoundary ) {
				e->winding = 0;
			} else {
				tessMeshDelete(mesh, e);
			}
		}
	}
}

void* heapAlloc( void* userData, unsigned int size )
{
	TESS_NOTUSED( userData );
	return malloc( size );
}

void* heapRealloc( void *userData, void* ptr, unsigned int size )
{
	TESS_NOTUSED( userData );
	return realloc( ptr, size );
}

void heapFree( void* userData, void* ptr )
{
	TESS_NOTUSED( userData );
	free( ptr );
}

static TESSalloc defaulAlloc =
{
	heapAlloc,
	heapRealloc,
	heapFree,
	0,
	0,
	0,
	0,
	0,
	0,
	0,
};

TESStesselator* tessNewTess( TESSalloc* alloc )
{
	TESStesselator* tess;

	if (alloc == NULL)
		alloc = &defaulAlloc;

	/* Only initialize fields which can be changed by the api.  Other fields
	* are initialized where they are used.
	*/

	tess = (TESStesselator *)alloc->memalloc( alloc->userData, sizeof( TESStesselator ));
	if ( tess == NULL ) {
		return 0;          /* out of memory */
	}
	tess->alloc = *alloc;
	/* Check and set defaults. */
	if (tess->alloc.meshEdgeBucketSize == 0)
		tess->alloc.meshEdgeBucketSize = 512;
	if (tess->alloc.meshVertexBucketSize == 0)
		tess->alloc.meshVertexBucketSize = 512;
	if (tess->alloc.meshFaceBucketSize == 0)
		tess->alloc.meshFaceBucketSize = 256;
	if (tess->alloc.dictNodeBucketSize == 0)
		tess->alloc.dictNodeBucketSize = 512;
	if (tess->alloc.regionBucketSize == 0)
		tess->alloc.regionBucketSize = 256;

	tess->normal[0] = 0;
	tess->normal[1] = 0;
	tess->normal[2] = 0;

	tess->bmin[0] = 0;
	tess->bmin[1] = 0;
	tess->bmax[0] = 0;
	tess->bmax[1] = 0;

	tess->reverseContours = 0;
    
	tess->windingRule = TESS_WINDING_ODD;
	tess->processCDT = 0;

	if (tess->alloc.regionBucketSize < 16)
		tess->alloc.regionBucketSize = 16;
	if (tess->alloc.regionBucketSize > 4096)
		tess->alloc.regionBucketSize = 4096;
	tess->regionPool = createBucketAlloc( &tess->alloc, "Regions",
										 sizeof(ActiveRegion), tess->alloc.regionBucketSize );

	// Initialize to begin polygon.
	tess->mesh = NULL;

	tess->vertexIndexCounter = 0;

	tess->vertices = NULL;
	tess->vertexIndices = NULL;
	tess->vertexCount = 0;
	tess->elements = NULL;
	tess->elementCount = 0;

	return tess;
}

void tessDeleteTess( TESStesselator *tess )
{

	struct TESSalloc alloc = tess->alloc;

	deleteBucketAlloc( tess->regionPool );

	if( tess->mesh != NULL ) {
		tessMeshDeleteMesh( &alloc, tess->mesh );
		tess->mesh = NULL;
	}
	if (tess->vertices != NULL) {
		alloc.memfree( alloc.userData, tess->vertices );
		tess->vertices = NULL;
	}
	if (tess->vertexIndices != NULL) {
		alloc.memfree( alloc.userData, tess->vertexIndices );
		tess->vertexIndices = NULL;
	}
	if (tess->elements != NULL) {
		alloc.memfree( alloc.userData, tess->elements );
		tess->elements = NULL;
	}

	alloc.memfree( alloc.userData, tess );
}


static TESSindex GetNeighbourFace(TESShalfEdge* edge)
{
	if (!edge->Rface)
		return TESS_UNDEF;
	if (!edge->Rface->inside)
		return TESS_UNDEF;
	return edge->Rface->n;
}

void OutputPolymesh( TESStesselator *tess, TESSmesh *mesh, int elementType, int polySize, int vertexSize )
{
	TESSvertex* v = 0;
	TESSface* f = 0;
	TESShalfEdge* edge = 0;
	int maxFaceCount = 0;
	int maxVertexCount = 0;
	int faceVerts, i;
	TESSindex *elements = 0;
	TESSreal *vert;

	// Assume that the input data is triangles now.
	// Try to merge as many polygons as possible
	if (polySize > 3)
	{
		tessMeshMergeConvexFaces(mesh, polySize);
	}

	// Mark unused
	for ( v = mesh->vHead.next; v != &mesh->vHead; v = v->next )
		v->n = TESS_UNDEF;

	// Create unique IDs for all vertices and faces.
	for ( f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		f->n = TESS_UNDEF;
		if( !f->inside ) continue;

		edge = f->anEdge;
		faceVerts = 0;
		do
		{
			v = edge->Org;
			if ( v->n == TESS_UNDEF )
			{
				v->n = maxVertexCount;
				maxVertexCount++;
			}
			faceVerts++;
			edge = edge->Lnext;
		}
		while (edge != f->anEdge);

		assert( faceVerts <= polySize );

		f->n = maxFaceCount;
		++maxFaceCount;
	}

	tess->elementCount = maxFaceCount;
	if (elementType == TESS_CONNECTED_POLYGONS)
		maxFaceCount *= 2;

	tess->elements = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData, sizeof(TESSindex) * maxFaceCount * polySize );

	tess->vertexCount = maxVertexCount;
	tess->vertices = (TESSreal*)tess->alloc.memalloc( tess->alloc.userData, sizeof(TESSreal) * tess->vertexCount * vertexSize );
	
	tess->vertexIndices = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData, sizeof(TESSindex) * tess->vertexCount );
	
	const float flip = tess->normal[2] > 0 ? 1.0f : -1.0f;
	// Output vertices.
	for (v = mesh->vHead.next; v != &mesh->vHead; v = v->next)
	{
		if (v->n != TESS_UNDEF)
		{
			// Store coordinate
			vert = &tess->vertices[v->n * vertexSize];
			vert[0] = v->s;
			vert[1] = flip * v->t;

			if (flip < 0)
			{
				assert(vert[1] == -v->t);
			}
			else
			{
				assert(vert[1] == v->t);
			}

			// Store vertex index.
			tess->vertexIndices[v->n] = v->idx;
		}
	}
	
	// Output indices.
	elements = tess->elements;
	for ( f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		if ( !f->inside ) continue;

		// Store polygon
		edge = f->anEdge;
		faceVerts = 0;
		do
		{
			v = edge->Org;
			*elements++ = v->n;
			faceVerts++;
			edge = edge->Lnext;
		}
		while (edge != f->anEdge);
		// Fill unused.
		for (i = faceVerts; i < polySize; ++i)
			*elements++ = TESS_UNDEF;

		// Store polygon connectivity
		if ( elementType == TESS_CONNECTED_POLYGONS )
		{
			edge = f->anEdge;
			do
			{
				*elements++ = GetNeighbourFace( edge );
				edge = edge->Lnext;
			}
			while (edge != f->anEdge);
			// Fill unused.
			for (i = faceVerts; i < polySize; ++i)
				*elements++ = TESS_UNDEF;
		}
	}
}

void OutputContours( TESStesselator *tess, TESSmesh *mesh, int vertexSize )
{
	TESShalfEdge *edge = NULL;
	TESShalfEdge *start = NULL;
	 
	tess->vertexCount = 0;
	tess->elementCount = 0;

	for (TESSface* f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		if ( !f->inside ) continue;

		start = edge = f->anEdge;
		do
		{
			++tess->vertexCount;
			edge = edge->Lnext;
		}
		while ( edge != start );

		++tess->elementCount;
	}

	tess->elements = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData, sizeof(TESSindex) * tess->elementCount * 2 );
	tess->vertices = (TESSreal*)tess->alloc.memalloc( tess->alloc.userData, sizeof(TESSreal) * tess->vertexCount * vertexSize );
	tess->vertexIndices = (TESSindex*)tess->alloc.memalloc( tess->alloc.userData, sizeof(TESSindex) * tess->vertexCount );
	
	TESSreal* verts = tess->vertices;
	TESSindex* elements = tess->elements;
	TESSindex* vertInds = tess->vertexIndices;

	int startVert = 0;
	int vertCount = 0;

	const float flip = tess->normal[2] > 0 ? 1.0f : -1.0f;

	for (TESSface* f = mesh->fHead.next; f != &mesh->fHead; f = f->next )
	{
		if ( !f->inside ) continue;

		vertCount = 0;
		start = edge = f->anEdge;
		do
		{
			*verts++ = edge->Org->s;
			*verts++ = flip * edge->Org->t;
			*vertInds++ = edge->Org->idx;
			++vertCount;
			edge = edge->Lnext;
		}
		while ( edge != start );

		elements[0] = startVert;
		elements[1] = vertCount;
		elements += 2;

		startVert += vertCount;
	}
}

void tessAddContour( TESStesselator *tess, int size, const void* vertices,
					int stride, int numVertices )
{
	const unsigned char *src = (const unsigned char*)vertices;
	
	if ( tess->mesh == NULL )
	  	tess->mesh = tessMeshNewMesh( &tess->alloc );

	if ( size < 2 )	size = 2;
	if ( size > 3 )	size = 3;

	const int winding = tess->reverseContours ? -1 : 1;
	
	TESShalfEdge* e = NULL;
	
	for(int i = 0; i < numVertices; ++i )
	{
		const TESSreal* coords = (const TESSreal*)src;
		src += stride;

		if( e == NULL ) {
			/* Make a self-loop (one vertex, one edge). */
			e = tessMeshMakeEdge( tess->mesh );
			
			tessMeshSplice(tess->mesh, e, e->Sym);
		} 
		else {
			/* Create a new vertex and edge which immediately follow e
			* in the ordering around the left face.
			*/
			tessMeshSplitEdge(tess->mesh, e);
				
			e = e->Lnext;
		}

		/* The new vertex is now e->Org. */
		e->Org->s = coords[0];
		e->Org->t = coords[1];
		
		/* Store the insertion number so that the vertex can be later recognized. */
		e->Org->idx = tess->vertexIndexCounter++;

		/* The winding of an edge says how the winding number changes as we
		* cross from the edge''s right face to its left face.  We add the
		* vertices in such an order that a CCW contour will add +1 to
		* the winding number of the region inside the contour.
		*/
        e->winding = winding;
		e->Sym->winding = -winding;
	}
}

void tessSetOption( TESStesselator *tess, int option, int value )
{
	switch(option)
	{
	case TESS_CONSTRAINED_DELAUNAY_TRIANGULATION:
		tess->processCDT = value > 0 ? 1 : 0;
		break;
	case TESS_REVERSE_CONTOURS:
		tess->reverseContours = value > 0 ? 1 : 0;
		break;
	}
}

//#include <windows.h>
int tessTesselate( TESStesselator *tess, int windingRule, int elementType,
				  int polySize, int vertexSize, const TESSreal* normal )
{
	//LARGE_INTEGER tt[7];
	//QueryPerformanceCounter(&tt[0]);
	

	TESSmesh *mesh;

	if (tess->vertices != NULL) {
		tess->alloc.memfree( tess->alloc.userData, tess->vertices );
		tess->vertices = 0;
	}
	if (tess->elements != NULL) {
		tess->alloc.memfree( tess->alloc.userData, tess->elements );
		tess->elements = 0;
	}
	if (tess->vertexIndices != NULL) {
		tess->alloc.memfree( tess->alloc.userData, tess->vertexIndices );
		tess->vertexIndices = 0;
	}

	tess->vertexIndexCounter = 0;

	if (normal)
	{
		tess->normal[0] = normal[0];
		tess->normal[1] = normal[1];
		tess->normal[2] = normal[2];
	}

	tess->windingRule = windingRule;

	if (vertexSize < 2)
		vertexSize = 2;
	if (vertexSize > 3)
		vertexSize = 3;

	if (setjmp(tess->env) != 0) {
		/* come back here if out of memory */
		return 0;
	}

	if (!tess->mesh)
	{
		return 0;
	}

	/* Determine the polygon normal and project vertices onto the plane
	* of the polygon.
	*/
	//tessProjectPolygon( tess );
	jerryPrepareVertices(tess);

	//QueryPerformanceCounter(&tt[1]);

	/* tessComputeInterior( tess ) computes the planar arrangement specified
	* by the given contours, and further subdivides this arrangement
	* into regions.  Each region is marked "inside" if it belongs
	* to the polygon, according to the rule given by tess->windingRule.
	* Each interior region is guaranteed be monotone.
	*/
	if ( !tessComputeInterior( tess ) ) {
		longjmp(tess->env,1);  /* could've used a label */
	}

	mesh = tess->mesh;

	//QueryPerformanceCounter(&tt[2]);

	/* If the user wants only the boundary contours, we throw away all edges
	* except those which separate the interior from the exterior.
	* Otherwise we tessellate all the regions marked "inside".
	*/
	if (elementType == TESS_BOUNDARY_CONTOURS)
	{
		tessMeshSetWindingNumber( mesh, 1, TRUE );
	} 
	else 
	{
		tessMeshTessellateInterior( mesh );
		if (tess->processCDT != 0)
			tessMeshRefineDelaunay( mesh, &tess->alloc );
	}
	
	//QueryPerformanceCounter(&tt[3]);

	tessMeshCheckMesh( mesh );

	//QueryPerformanceCounter(&tt[4]);

	if (elementType == TESS_BOUNDARY_CONTOURS) {
		OutputContours( tess, mesh, vertexSize );     /* output contours */
	}
	else
	{
		OutputPolymesh( tess, mesh, elementType, polySize, vertexSize );     /* output polygons */
	}

	//QueryPerformanceCounter(&tt[5]);

	tessMeshDeleteMesh( &tess->alloc, mesh );
	tess->mesh = NULL;


	//QueryPerformanceCounter(&tt[6]);

	/*LARGE_INTEGER b;
	QueryPerformanceFrequency(&b);

	double aa[6];
	for (int i = 0; i < 6; i++)
	{
		aa[i] = ((double)(tt[i+1].QuadPart - tt[i].QuadPart)*1000) / b.QuadPart;
	}
	
	int abba = 10;*/


	return 1;
}

int tessGetVertexCount( TESStesselator *tess )
{
	return tess->vertexCount;
}

const TESSreal* tessGetVertices( TESStesselator *tess )
{
	return tess->vertices;
}

const TESSindex* tessGetVertexIndices( TESStesselator *tess )
{
	return tess->vertexIndices;
}

int tessGetElementCount( TESStesselator *tess )
{
	return tess->elementCount;
}

const int* tessGetElements( TESStesselator *tess )
{
	return tess->elements;
}
