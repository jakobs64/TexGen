/*=============================================================================
TexGen: Geometric textile modeller.
Copyright (C) 2006 Martin Sherburn

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
=============================================================================*/

#include "PrecompiledHeaders.h"
#include "Mesh.h"
#include "Plane.h"
#if BUILD_OCTREE
#include "MeshOctreeClasses.h"
#endif
#include "Textile.h"
#include "TexGen.h"

using namespace TexGen;

auto LessPairDoubleXYZ = [](std::pair<double, XYZ> x, std::pair<double, XYZ> y) { return x.first < y.first; };

extern "C"
{
#include "../Triangle/triangle.h"
}

CMesh::CMesh(void)
{
}

CMesh::~CMesh(void)
{
}

CMesh::CMesh(TiXmlElement &Element)
{
	FOR_EACH_TIXMLELEMENT(pNode, Element, "Node")
	{
		m_Nodes.push_back(valueify<XYZ>(pNode->Attribute("value")));
	}
	FOR_EACH_TIXMLELEMENT(pType, Element, "Element")
	{
		ELEMENT_TYPE iType = (ELEMENT_TYPE)valueify<int>(pType->Attribute("type"));
		FOR_EACH_TIXMLELEMENT(pIndex, *pType, "Index")
		{
			m_Indices[iType].push_back(valueify<int>(pIndex->Attribute("value")));
		}
	}
}

void CMesh::PopulateTiXmlElement(TiXmlElement &Element, OUTPUT_TYPE OutputType) const
{
	vector<XYZ>::const_iterator itNode;
	for (itNode = m_Nodes.begin(); itNode != m_Nodes.end(); ++itNode)
	{
		TiXmlElement Node("Node");
		Node.SetAttribute("value", stringify(*itNode));
		Element.InsertEndChild(Node);
	}
	int i;
	list<int>::const_iterator itIndex;
	for (i=0; i<NUM_ELEMENT_TYPES; ++i)
	{
		TiXmlElement Indices("Element");
		Indices.SetAttribute("type", i);
		for (itIndex = m_Indices[i].begin(); itIndex != m_Indices[i].end(); ++itIndex)
		{
			TiXmlElement Index("Index");
			Index.SetAttribute("value", *itIndex);
			Indices.InsertEndChild(Index);
		}
		Element.InsertEndChild(Indices);
	}
}

int CMesh::InsertNodes(const CMesh &Mesh, XYZ Offset)
{
	int iOffset = (int)m_Nodes.size();
	vector<XYZ>::const_iterator itNode;
	for (itNode = Mesh.m_Nodes.begin(); itNode != Mesh.m_Nodes.end(); ++itNode)
	{
		m_Nodes.push_back(*itNode+Offset);
	}
	return iOffset;
}

void CMesh::InsertMesh(const CMesh &Mesh, XYZ Offset)
{
	int iOffset = InsertNodes(Mesh, Offset);
	int i;
	list<int>::const_iterator itIndex;
	for (i=0; i<NUM_ELEMENT_TYPES; ++i)
	{
		for (itIndex = Mesh.m_Indices[i].begin(); itIndex != Mesh.m_Indices[i].end(); ++itIndex)
		{
			m_Indices[i].push_back(*itIndex + iOffset);
		}
	}
}

void CMesh::ChangeNodeIndices(int iChangeTo, int iChangeFrom)
{
	int i;
	list<int>::iterator itIndex;
	for (i=0; i<NUM_ELEMENT_TYPES; ++i)
	{
		for (itIndex = m_Indices[i].begin(); itIndex != m_Indices[i].end(); ++itIndex)
		{
			if (*itIndex == iChangeFrom)
			{
				*itIndex = iChangeTo;
			}
		}
	}
}

void CMesh::ChangeNodeIndices(int iChangeTo, int iChangeFrom, vector<vector<int*> > &References)
{
	vector<int*>::iterator itIntPointer;
	vector<int*> &ChangeFromPointers = References[iChangeFrom];
	vector<int*> &ChangeToPointers = References[iChangeTo];
	for (itIntPointer = ChangeFromPointers.begin(); itIntPointer != ChangeFromPointers.end(); ++itIntPointer)
	{
		**itIntPointer = iChangeTo;
		ChangeToPointers.push_back(*itIntPointer);
	}
	ChangeFromPointers.clear();
}

void CMesh::RemoveOpposingTriangles()
{
	list<int> &TriangleIndices = m_Indices[TRI];
	int i1[3];
	int i2[3];
	list<int>::iterator itIter1, itIter2;
	list<int>::iterator itTriStart1, itTriStart2;
//	int iCommonNodes;
//	int i, j;
	for (itIter1 = TriangleIndices.begin(); itIter1 != TriangleIndices.end(); )
	{
		itTriStart1 = itIter1;
		i1[0] = *(itIter1++);
		i1[1] = *(itIter1++);
		i1[2] = *(itIter1++);
		for (itIter2 = itIter1; itIter2 != TriangleIndices.end(); )
		{
			itTriStart2 = itIter2;
			i2[0] = *(itIter2++);
			i2[1] = *(itIter2++);
			i2[2] = *(itIter2++);

			if (i1[0] == i2[0] && i1[1] == i2[2] && i1[2] == i2[1] ||
				i1[0] == i2[1] && i1[1] == i2[0] && i1[2] == i2[2] ||
				i1[0] == i2[2] && i1[1] == i2[1] && i1[2] == i2[0])
			{
				TriangleIndices.erase(itTriStart1, itIter1);
				if (itIter1 == itTriStart2)
					itIter1 = TriangleIndices.erase(itTriStart2, itIter2);
				else
					TriangleIndices.erase(itTriStart2, itIter2);
				//continue;
				break;
			}
		}		
	}
}

void CMesh::RemoveOpposingQuads()
{
	list<int> &QuadIndices = m_Indices[QUAD];
	int i1[4];
	int i2[4];
	list<int>::iterator itIter1, itIter2;
	list<int>::iterator itQuadStart1, itQuadStart2;
//	int iCommonNodes;
//	int i, j;
	for (itIter1 = QuadIndices.begin(); itIter1 != QuadIndices.end(); )
	{
		itQuadStart1 = itIter1;
		i1[0] = *(itIter1++);
		i1[1] = *(itIter1++);
		i1[2] = *(itIter1++);
		i1[3] = *(itIter1++);
		for (itIter2 = itIter1; itIter2 != QuadIndices.end(); )
		{
			itQuadStart2 = itIter2;
			i2[0] = *(itIter2++);
			i2[1] = *(itIter2++);
			i2[2] = *(itIter2++);
			i2[3] = *(itIter2++);

			if (i1[0] == i2[0] && i1[1] == i2[3] && i1[2] == i2[2] && i1[3] == i2[1] ||
				i1[0] == i2[1] && i1[1] == i2[0] && i1[2] == i2[3] && i1[3] == i2[2] ||
				i1[0] == i2[2] && i1[1] == i2[1] && i1[2] == i2[0] && i1[3] == i2[3] ||
				i1[0] == i2[3] && i1[1] == i2[2] && i1[2] == i2[1] && i1[3] == i2[0])
			{
				QuadIndices.erase(itQuadStart1, itIter1);
				if (itIter1 == itQuadStart2)
					itIter1 = QuadIndices.erase(itQuadStart2, itIter2);
				continue;
			}
		}		
	}
}

void CMesh::RemoveDegenerateTriangles()
{
	list<int> &TriangleIndices = m_Indices[TRI];
	list<int>::iterator itIter, itTriStart;
	int i1;
	int i2;
	int i3;
	for (itIter = TriangleIndices.begin(); itIter != TriangleIndices.end(); )
	{
		itTriStart = itIter;
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		if (i1 == i2 || i2 == i3 || i3 == i1)
		{
			itIter = TriangleIndices.erase(itTriStart, itIter);
		}
	}
}

void CMesh::RemoveDuplicateElements(CMesh::ELEMENT_TYPE ElementType)
{
	list<int> &ElementIndices = m_Indices[ElementType];
	vector<int> Index1;
	vector<int> Index2;
	list<int>::iterator itIter1, itIter2;
	list<int>::iterator itStart1, itStart2;

	for (itIter1 = ElementIndices.begin(); itIter1 != ElementIndices.end(); )
	{
		itStart1 = itIter1;
		Index1.clear();
		for ( int i=0; i < GetNumNodes(ElementType); ++i )
		{
			Index1.push_back( *(itIter1++));
		}
		sort( Index1.begin(), Index1.end() ); 
		for (itIter2 = itIter1; itIter2 != ElementIndices.end(); )
		{
			itStart2 = itIter2;
			Index2.clear();
			for ( int j=0; j < GetNumNodes(ElementType); ++j )
			{
				Index2.push_back( *(itIter2++));
			}
			sort( Index2.begin(), Index2.end() );

			if ( equal( Index1.begin(), Index1.end(), Index2.begin() ) )
			{
				if (itIter1 == itStart2)
				{
					itIter1 = ElementIndices.erase(itStart2, itIter2);
					itIter2 = itIter1;
				}
				else
					itIter2 = ElementIndices.erase(itStart2, itIter2);
				//continue;
			}
		}		
	}
}

void CMesh::RemoveDuplicateTriangles()
{
	list<int> &TriangleIndices = m_Indices[TRI];
	int i1[3];
	int i2[3];
	list<int>::iterator itIter1, itIter2;
	list<int>::iterator itTriStart1, itTriStart2;
//	int iCommonNodes;
//	int i, j;
	for (itIter1 = TriangleIndices.begin(); itIter1 != TriangleIndices.end(); )
	{
		itTriStart1 = itIter1;
		i1[0] = *(itIter1++);
		i1[1] = *(itIter1++);
		i1[2] = *(itIter1++);
		for (itIter2 = itIter1; itIter2 != TriangleIndices.end(); )
		{
			itTriStart2 = itIter2;
			i2[0] = *(itIter2++);
			i2[1] = *(itIter2++);
			i2[2] = *(itIter2++);

			if ( i1[0] == i2[0] && i1[1] == i2[1] && i1[2] == i2[2] )
			{
				//TriangleIndices.erase(itTriStart1, itIter1);
				if (itIter1 == itTriStart2)
				{
					itIter1 = TriangleIndices.erase(itTriStart2, itIter2);
					itIter2 = itIter1;
				}
				else
					itIter2 = TriangleIndices.erase(itTriStart2, itIter2);
				//continue;
			}
		}		
	}
}

void CMesh::RemoveDuplicateSegments()
{
	list<int> &LineIndices = m_Indices[LINE];
	list<int>::iterator itIter, itLineStart;
	list<int>::iterator itCompare;
	int i1, i2;
	int j1, j2;
	for (itIter = LineIndices.begin(); itIter != LineIndices.end(); )
	{
		itLineStart = itIter;
		i1 = *(itIter++);
		i2 = *(itIter++);
		for (itCompare = itIter; itCompare != LineIndices.end(); )
		{
			j1 = *(itCompare++);
			j2 = *(itCompare++);
			if ((i1 == j1 && i2 == j2) ||
				(i1 == j2 && i2 == j1))
			{
				itIter = LineIndices.erase(itLineStart, itIter);
				break;
			}
		}
	}
}

pair<XYZ, XYZ> CMesh::GetAABB(double dGrowDistance) const
{
	pair<XYZ, XYZ> AABB;

	vector<XYZ>::const_iterator itNode;
	int iNode;
	for (itNode = m_Nodes.begin(), iNode=0; itNode != m_Nodes.end(); ++itNode, ++iNode)
	{
		if (iNode == 0)
		{
			AABB.first = AABB.second = *itNode;
		}
		else
		{
			AABB.first = Min(AABB.first, *itNode);
			AABB.second = Max(AABB.second, *itNode);
		}
	}

	AABB.first.x -= dGrowDistance;
	AABB.first.y -= dGrowDistance;
	AABB.first.z -= dGrowDistance;
	AABB.second.x += dGrowDistance;
	AABB.second.y += dGrowDistance;
	AABB.second.z += dGrowDistance;

	return AABB;
}

#if BUILD_OCTREE
int CMesh::MergeNodes(const double TOL)
{
	// Merging of nodes is optimised using an octree, that is the space is partitioned to group
	// nodes that are in proximity to each other. Which reduces the number of checks needed to be performed
	// when comparing node positions.

	// Get the dimensions of the octree
	pair<XYZ, XYZ> AABB = GetAABB(TOL);

	double dSize = AABB.second.x - AABB.first.x;
	if (dSize < AABB.second.y - AABB.first.y)
		dSize = AABB.second.y - AABB.first.y;
	if (dSize < AABB.second.z - AABB.first.z)
		dSize = AABB.second.z - AABB.first.z;

	// Prevent problems with creating an octree too small
	dSize = max(dSize, 1e-3);

	// Create the octree
	XYZ Avg = 0.5*(AABB.second + AABB.first);
	dSize *= 1.1;
	Octree<pair<int, XYZ> > Octree(Vector3f(float(Avg.x-0.5*dSize), float(Avg.y-0.5*dSize), float(Avg.z-0.5*dSize)), (float)dSize, 10, 10);
	COctreeAgentNode Agent(TOL);

	// Add the nodes to the octree, the indices of the nodes need to be stored with the nodes because
	// the elements refer to the nodes by their index.
	vector<pair<int, XYZ> > NodesWithIndices;
	NodesWithIndices.resize(m_Nodes.size());

	vector<XYZ>::iterator itNode;
	int iNode;
	for (itNode = m_Nodes.begin(), iNode=0; itNode != m_Nodes.end(); ++itNode, ++iNode)
	{
		NodesWithIndices[iNode] = pair<int, XYZ>(iNode, *itNode);
		Octree.insertItem(NodesWithIndices[iNode], Agent);
	}

	// Visit the octree with the visitor class, this actually does the node merging
	COctreeVisitorMergeNodes Visitor(*this, TOL);
	Octree.visit(Visitor);

	return Visitor.GetNumMerged();
}
#else
int CMesh::MergeNodes(const double /*TOL*/)
{
	// Octree support disabled: no-op
	return 0;
}
#endif

vector<pair<int, int> > CMesh::GetNodePairs(XYZ Vector, const double TOL) const
{
	vector<pair<int, int> > NodePairs;
	vector<XYZ>::const_iterator itNode1;
	vector<XYZ>::const_iterator itNode2;
	int iNode1, iNode2;
	double dTolSquared = TOL*TOL;
	for (itNode1 = m_Nodes.begin(), iNode1=0; itNode1 != m_Nodes.end(); ++itNode1, ++iNode1)
	{
		for (itNode2 = m_Nodes.begin(), iNode2=0; itNode2 != m_Nodes.end(); ++itNode2, ++iNode2)
		{
			if (GetLengthSquared(*itNode1, *itNode2-Vector) < dTolSquared)
			{
				NodePairs.push_back(make_pair(iNode1, iNode2));
				break;
			}
		}
	}
	return NodePairs;
}

void CMesh::GetNodePairs(XYZ Vector, vector<pair<int, int> > &NodePairs, const double TOL ) const
{
	vector<XYZ>::const_iterator itNode1;
	vector<XYZ>::const_iterator itNode2;

	vector<XYZ> OffsetNodes;
	for ( itNode2 = m_Nodes.begin(); itNode2 != m_Nodes.end(); ++itNode2 )
	{
		OffsetNodes.push_back( *itNode2 - Vector );
	}

	bool bProgress = m_Nodes.size() > 10000;
	int iProgressInc = m_Nodes.size()/10;

	int iNode1, iNode2;
	double dTolSquared = TOL*TOL;
	bool bFound;
	for (itNode1 = m_Nodes.begin(), iNode1=0; itNode1 != m_Nodes.end(); ++itNode1, ++iNode1)
	{
		if ( bProgress && !(iNode1 % iProgressInc) )
		{
			TGLOG( "GetNodePairs " << iNode1/iProgressInc*10 << "% complete");
		}
		bFound = false;
		for ( iNode2 = iNode1; iNode2 < OffsetNodes.size(); ++ iNode2 )  // Matching node normally later in mesh than itNode1 so start checking here
		{                                                                // to improve speed
			if (GetLengthSquared(*itNode1, OffsetNodes[iNode2]) < dTolSquared)
			{
				NodePairs.push_back(make_pair(iNode1, iNode2));
				bFound = true;
				break;
			}
		}
		if ( !bFound )
		{
			for ( iNode2 = 0; iNode2 < iNode1; ++iNode2 )   // Check from beginning if node earlier in mesh than one checking against
			{
				if (GetLengthSquared(*itNode1, OffsetNodes[iNode2]) < dTolSquared)
				{
					NodePairs.push_back(make_pair(iNode1, iNode2));
					break;
				}
			}
		}
	}
}

int CMesh::GetClosestNode(XYZ Position) const
{
	vector<XYZ>::const_iterator itNode;
	double dClosestDistSqrd = -1;
	double dDistSqrd;
	int iClosestNode = -1, i;
	for (itNode = m_Nodes.begin(), i=0; itNode != m_Nodes.end(); ++itNode, ++i)
	{
		dDistSqrd = GetLengthSquared(Position, *itNode);
		if ( dDistSqrd < dClosestDistSqrd || dClosestDistSqrd == -1 )
		{
			dClosestDistSqrd = dDistSqrd;
			iClosestNode = i;
		}
	}
	return iClosestNode;
}

int CMesh::GetClosestNodeDistance(XYZ Position, double dTol ) const
{
	vector<XYZ>::const_iterator itNode;
	double dClosestDistSqrd = -1;
	double dDistSqrd;
	int iClosestNode = -1, i;
	for (itNode = m_Nodes.begin(), i=0; itNode != m_Nodes.end(); ++itNode, ++i)
	{
		dDistSqrd = GetLengthSquared(Position, *itNode);
		if ( dDistSqrd < dClosestDistSqrd || dClosestDistSqrd == -1 )
		{
			dClosestDistSqrd = dDistSqrd;
			if ( dClosestDistSqrd <= dTol )
				iClosestNode = i;
		}
	}
	return iClosestNode;
}

void CMesh::ConvertToSurfaceMesh()
{
	Convert3Dto2D();
	RemoveOpposingTriangles();
	RemoveOpposingQuads();
}

void CMesh::Convert3Dto2D()
{
	ConvertHextoQuad();
	ConvertWedgeto2D();
	ConvertTettoTriangle();
	ConvertPyramidto2D();
}

void CMesh::ConvertToTetMesh()
{
	ConvertWedgetoTet();
	ConvertHextoTet();
	ConvertPyramidtoTet();
	RemoveAllElementsExcept(TET);
}

void CMesh::ConvertToTriangleMesh()
{
	Convert3Dto2D();
	ConvertQuadstoTriangles();
	RemoveAllElementsExcept(TRI);
}

void CMesh::ConvertToSegmentMesh()
{
	ConvertToTriangleMesh();
	ConvertTrianglestoSegments();
	RemoveDuplicateSegments();
}

void CMesh::ConvertTriToQuad( double Tolerance )
{
	list<int> &TriangleIndices = m_Indices[TRI];
	vector<int> i1;
	int i2[3];
	list<int>::iterator itIter1, itIter2;
	list<int>::iterator itTriStart1, itTriStart2;
	vector<int>::iterator iti1;

	int iNumNodes = GetNumNodes(TRI);
	for (itIter1 = TriangleIndices.begin(); itIter1 != TriangleIndices.end(); )
	{
		itTriStart1 = itIter1;
		i1.clear();
		i1.push_back(*(itIter1++));
		i1.push_back(*(itIter1++));
		i1.push_back(*(itIter1++));
		for ( itIter2 = itIter1; itIter2 != TriangleIndices.end(); )
		{
			itTriStart2 = itIter2;
			list<int> RemInd;
			i2[0] = *(itIter2++);
			RemInd.push_back(i2[0]);
			i2[1] = *(itIter2++);
			RemInd.push_back(i2[1]);
			i2[2] = *(itIter2++);
			RemInd.push_back(i2[2]);
			vector<int> CommonInd;
			
			int i;
			for ( iti1 = i1.begin(), i = 0; iti1 != i1.end(); iti1++, ++i )
			{
				
				for ( int j = 0; j < iNumNodes; ++j )
				{
					if ( *iti1 == i2[j] )
					{
						CommonInd.push_back(i);
						RemInd.remove(i2[j]);
						break;
					}
				}
			}
			if ( CommonInd.size() == 2 )
			{
				// Check normals
				XYZ N1 = m_Nodes[i1[0]];
				XYZ N2 = m_Nodes[i1[1]];
				XYZ N3 = m_Nodes[i1[2]];
				XYZ Normal1 = CrossProduct(N2-N1, N3-N1);
				Normal1 /= GetLength(Normal1);  // Normalise

				N1 = m_Nodes[i2[0]];
				N2 = m_Nodes[i2[1]];
				N3 = m_Nodes[i2[2]];
				XYZ Normal2 = CrossProduct(N2-N1, N3-N1);
				Normal2 /= GetLength(Normal2);  // Normalise

				if ( Normal1 == Normal2 || Normal1 == -Normal2 )
				{
					if ( CommonInd[0] == 0 && CommonInd[1] == 2 )
						i1.push_back( *(RemInd.begin()) );
					else
						i1.insert( i1.begin()+CommonInd[1],*(RemInd.begin()));
					AddElement(	QUAD, i1 );
					if ( itIter1 == itTriStart2 )
						itIter1 = TriangleIndices.erase( itTriStart1, itIter2 );
					else
					{
						TriangleIndices.erase( itTriStart2, itIter2 );
						itIter1 = TriangleIndices.erase( itTriStart1, itTriStart1 );
					}
					
					break;
				}
				
			}
		}
	}
}

int CMesh::RemoveUnreferencedNodes()
{
	set<int> UnreferencedNodes;
	int i;
	for (i=0; i<(int)m_Nodes.size(); ++i)
	{
		UnreferencedNodes.insert(i);
	}

	list<int>::iterator itIndex;
	for (i = 0; i < NUM_ELEMENT_TYPES; ++i)
	{
		for (itIndex = m_Indices[i].begin(); itIndex != m_Indices[i].end(); ++itIndex)
		{
			UnreferencedNodes.erase(*itIndex);   // Get rid of the nodes which are referenced..
		}
	}

	return DeleteNodes(UnreferencedNodes);   // and delete the unused ones remaining
}

void CMesh::RemoveAllElementsExcept(ELEMENT_TYPE Type)
{
	int i;
	for (i = 0; i < NUM_ELEMENT_TYPES; ++i)
	{
		if (i == Type)
			continue;
		else
			m_Indices[i].clear();
	}
}

void CMesh::RemoveElementType( ELEMENT_TYPE Type )
{
	m_Indices[Type].clear();
}

int CMesh::DeleteNodes(const set<int> &Nodes)
{
	vector<int> NewIndices;
	NewIndices.resize(m_Nodes.size(), 0);
	int i, j;
	int iNumNodesDeleted = 0;
	for (i=0, j=0; i<(int)NewIndices.size(); ++i)
	{
		NewIndices[i] = j;
		if (Nodes.count(i))
			++iNumNodesDeleted;
		else
			++j;
	}

	list<int>::iterator itIndex;
	for (i = 0; i < NUM_ELEMENT_TYPES; ++i)
	{
		for (itIndex = m_Indices[i].begin(); itIndex != m_Indices[i].end(); ++itIndex)
		{
			*itIndex = NewIndices[*itIndex];
		}
	}

	vector<XYZ>::iterator itNode;
	vector<XYZ> NewNodes;
	int iNodeIndex;
	for (itNode = m_Nodes.begin(), iNodeIndex=0; itNode != m_Nodes.end(); ++itNode, ++iNodeIndex)
	{
		if (!Nodes.count(iNodeIndex))
			NewNodes.push_back(*itNode);
	}
	m_Nodes = NewNodes;

	return iNumNodesDeleted;
}

void CMesh::ConvertHextoQuad()
{
	list<int>::iterator itIter;
	int i1, i2, i3, i4, i5, i6, i7, i8;
	for (itIter = m_Indices[HEX].begin(); itIter != m_Indices[HEX].end(); )
	{
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);
		i5 = *(itIter++);
		i6 = *(itIter++);
		i7 = *(itIter++);
		i8 = *(itIter++);

		m_Indices[QUAD].push_back(i4);
		m_Indices[QUAD].push_back(i3);
		m_Indices[QUAD].push_back(i2);
		m_Indices[QUAD].push_back(i1);

		m_Indices[QUAD].push_back(i1);
		m_Indices[QUAD].push_back(i2);
		m_Indices[QUAD].push_back(i6);
		m_Indices[QUAD].push_back(i5);

		m_Indices[QUAD].push_back(i2);
		m_Indices[QUAD].push_back(i3);
		m_Indices[QUAD].push_back(i7);
		m_Indices[QUAD].push_back(i6);

		m_Indices[QUAD].push_back(i3);
		m_Indices[QUAD].push_back(i4);
		m_Indices[QUAD].push_back(i8);
		m_Indices[QUAD].push_back(i7);

		m_Indices[QUAD].push_back(i4);
		m_Indices[QUAD].push_back(i1);
		m_Indices[QUAD].push_back(i5);
		m_Indices[QUAD].push_back(i8);

		m_Indices[QUAD].push_back(i5);
		m_Indices[QUAD].push_back(i6);
		m_Indices[QUAD].push_back(i7);
		m_Indices[QUAD].push_back(i8);
	}
	m_Indices[HEX].clear();
}

void CMesh::ConvertWedgeto2D()
{
	list<int>::iterator itIter;
	int i1, i2, i3, i4, i5, i6;
	for (itIter = m_Indices[WEDGE].begin(); itIter != m_Indices[WEDGE].end(); )
	{
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);
		i5 = *(itIter++);
		i6 = *(itIter++);

		m_Indices[TRI].push_back(i1);
		m_Indices[TRI].push_back(i2);
		m_Indices[TRI].push_back(i3);

		m_Indices[QUAD].push_back(i2);
		m_Indices[QUAD].push_back(i1);
		m_Indices[QUAD].push_back(i4);
		m_Indices[QUAD].push_back(i5);

		m_Indices[QUAD].push_back(i3);
		m_Indices[QUAD].push_back(i2);
		m_Indices[QUAD].push_back(i5);
		m_Indices[QUAD].push_back(i6);

		m_Indices[QUAD].push_back(i1);
		m_Indices[QUAD].push_back(i3);
		m_Indices[QUAD].push_back(i6);
		m_Indices[QUAD].push_back(i4);

		m_Indices[TRI].push_back(i6);
		m_Indices[TRI].push_back(i5);
		m_Indices[TRI].push_back(i4);
	}
	m_Indices[WEDGE].clear();
}

void CMesh::ConvertTettoTriangle()
{
	list<int>::iterator itIter;
	int i1, i2, i3, i4;
	for (itIter = m_Indices[TET].begin(); itIter != m_Indices[TET].end(); )
	{
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);

		m_Indices[TRI].push_back(i3);
		m_Indices[TRI].push_back(i2);
		m_Indices[TRI].push_back(i1);

		m_Indices[TRI].push_back(i1);
		m_Indices[TRI].push_back(i2);
		m_Indices[TRI].push_back(i4);

		m_Indices[TRI].push_back(i2);
		m_Indices[TRI].push_back(i3);
		m_Indices[TRI].push_back(i4);

		m_Indices[TRI].push_back(i3);
		m_Indices[TRI].push_back(i1);
		m_Indices[TRI].push_back(i4);
	}
	m_Indices[TET].clear();
}

void CMesh::ConvertHextoWedge(bool bQuality)
{
	list<int>::iterator itIter;
	int i0, i1, i2, i3, i4, i5, i6, i7;
	for (itIter = m_Indices[HEX].begin(); itIter != m_Indices[HEX].end(); )
	{
		i0 = *(itIter++);
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);
		i5 = *(itIter++);
		i6 = *(itIter++);
		i7 = *(itIter++);

		// Split it up randomly without regard to quality
		m_Indices[WEDGE].push_back(i4);
		m_Indices[WEDGE].push_back(i5);
		m_Indices[WEDGE].push_back(i6);
		m_Indices[WEDGE].push_back(i0);
		m_Indices[WEDGE].push_back(i1);
		m_Indices[WEDGE].push_back(i2);

		m_Indices[WEDGE].push_back(i6);
		m_Indices[WEDGE].push_back(i7);
		m_Indices[WEDGE].push_back(i4);
		m_Indices[WEDGE].push_back(i2);
		m_Indices[WEDGE].push_back(i3);
		m_Indices[WEDGE].push_back(i0);

		// Todo: If the bQuality flag is set then use a more intelligent
		// way to decide how to split the element up
	}
	m_Indices[HEX].clear();
}

void CMesh::ConvertWedgetoTetandPyramid(bool bQuality)
{
	list<int>::iterator itIter;
	int i0, i1, i2, i3, i4, i5;
	for (itIter = m_Indices[WEDGE].begin(); itIter != m_Indices[WEDGE].end(); )
	{
		i0 = *(itIter++);
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);
		i5 = *(itIter++);

		// Split it up randomly without regard to quality
		m_Indices[PYRAMID].push_back(i0);
		m_Indices[PYRAMID].push_back(i1);
		m_Indices[PYRAMID].push_back(i4);
		m_Indices[PYRAMID].push_back(i3);
		m_Indices[PYRAMID].push_back(i2);

		m_Indices[TET].push_back(i3);
		m_Indices[TET].push_back(i4);
		m_Indices[TET].push_back(i5);
		m_Indices[TET].push_back(i2);

		// Todo: If the bQuality flag is set then use a more intelligent
		// way to decide how to split the element up
	}
	m_Indices[WEDGE].clear();
}

void CMesh::ConvertPyramidtoTet(bool bQuality)
{
	list<int>::iterator itIter;
	int i0, i1, i2, i3, i4;
	for (itIter = m_Indices[PYRAMID].begin(); itIter != m_Indices[PYRAMID].end(); )
	{
		i0 = *(itIter++);
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);

		double d1 = 0, d2 = 0;
		if (bQuality)
		{
			// If we don't care about the quality then just skip this calculation
			d1 = GetLengthSquared(m_Nodes[i0], m_Nodes[i2]);
			d2 = GetLengthSquared(m_Nodes[i1], m_Nodes[i3]);
		}
		if (d1<d2)
		{
			m_Indices[TET].push_back(i0);
			m_Indices[TET].push_back(i1);
			m_Indices[TET].push_back(i2);
			m_Indices[TET].push_back(i4);

			m_Indices[TET].push_back(i0);
			m_Indices[TET].push_back(i2);
			m_Indices[TET].push_back(i3);
			m_Indices[TET].push_back(i4);
		}
		else
		{
			m_Indices[TET].push_back(i1);
			m_Indices[TET].push_back(i2);
			m_Indices[TET].push_back(i3);
			m_Indices[TET].push_back(i4);

			m_Indices[TET].push_back(i3);
			m_Indices[TET].push_back(i0);
			m_Indices[TET].push_back(i1);
			m_Indices[TET].push_back(i4);
		}
	}
	m_Indices[PYRAMID].clear();
}

void CMesh::ConvertWedgetoTet(bool bQuality)
{
	ConvertWedgetoTetandPyramid(bQuality);
	ConvertPyramidtoTet(bQuality);
}

void CMesh::ConvertHextoTet(bool bQuality)
{
	ConvertHextoWedge(bQuality);
	ConvertWedgetoTet(bQuality);
}

void CMesh::ConvertPyramidto2D()
{
	list<int>::iterator itIter;
	int i0, i1, i2, i3, i4;
	for (itIter = m_Indices[PYRAMID].begin(); itIter != m_Indices[PYRAMID].end(); )
	{
		i0 = *(itIter++);
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);
		i4 = *(itIter++);

		m_Indices[QUAD].push_back(i3);
		m_Indices[QUAD].push_back(i2);
		m_Indices[QUAD].push_back(i1);
		m_Indices[QUAD].push_back(i0);

		m_Indices[TRI].push_back(i0);
		m_Indices[TRI].push_back(i1);
		m_Indices[TRI].push_back(i4);

		m_Indices[TRI].push_back(i1);
		m_Indices[TRI].push_back(i2);
		m_Indices[TRI].push_back(i4);

		m_Indices[TRI].push_back(i2);
		m_Indices[TRI].push_back(i3);
		m_Indices[TRI].push_back(i4);

		m_Indices[TRI].push_back(i3);
		m_Indices[TRI].push_back(i0);
		m_Indices[TRI].push_back(i4);
	}
	m_Indices[PYRAMID].clear();
}

list<int>::iterator CMesh::ConvertQuadtoTriangles(list<int>::iterator itQuad)
{
	int i0, i1, i2, i3;

	list<int>::iterator itEnd = itQuad;

	i0 = *(itEnd++);
	i1 = *(itEnd++);
	i2 = *(itEnd++);
	i3 = *(itEnd++);

	m_Indices[TRI].push_back(i0);
	m_Indices[TRI].push_back(i1);
	m_Indices[TRI].push_back(i2);

	m_Indices[TRI].push_back(i0);
	m_Indices[TRI].push_back(i2);
	m_Indices[TRI].push_back(i3);

	return m_Indices[QUAD].erase(itQuad, itEnd);
}

void CMesh::ConvertTrianglestoSegments()
{
	list<int>::iterator itIter;
	int i0, i1, i2;
	for (itIter = m_Indices[TRI].begin(); itIter != m_Indices[TRI].end(); )
	{
		i0 = *(itIter++);
		i1 = *(itIter++);
		i2 = *(itIter++);

		m_Indices[LINE].push_back(i0);
		m_Indices[LINE].push_back(i1);

		m_Indices[LINE].push_back(i1);
		m_Indices[LINE].push_back(i2);

		m_Indices[LINE].push_back(i2);
		m_Indices[LINE].push_back(i0);
	}
	m_Indices[TRI].clear();
}

void CMesh::ConvertQuadstoTriangles(bool bQuality)
{
	list<int>::iterator itIter;
	int i0, i1, i2, i3;
	for (itIter = m_Indices[QUAD].begin(); itIter != m_Indices[QUAD].end(); )
	{
		i0 = *(itIter++);
		i1 = *(itIter++);
		i2 = *(itIter++);
		i3 = *(itIter++);

		if (!bQuality)
		{
			m_Indices[TRI].push_back(i0);
			m_Indices[TRI].push_back(i1);
			m_Indices[TRI].push_back(i2);

			m_Indices[TRI].push_back(i0);
			m_Indices[TRI].push_back(i2);
			m_Indices[TRI].push_back(i3);
		}
		else
		{
			double d1, d2;
			d1 = GetLengthSquared(m_Nodes[i0], m_Nodes[i2]);
			d2 = GetLengthSquared(m_Nodes[i1], m_Nodes[i3]);
			if (d1<d2)
			{
				m_Indices[TRI].push_back(i0);
				m_Indices[TRI].push_back(i1);
				m_Indices[TRI].push_back(i2);

				m_Indices[TRI].push_back(i0);
				m_Indices[TRI].push_back(i2);
				m_Indices[TRI].push_back(i3);
			}
			else
			{
				m_Indices[TRI].push_back(i1);
				m_Indices[TRI].push_back(i2);
				m_Indices[TRI].push_back(i3);

				m_Indices[TRI].push_back(i3);
				m_Indices[TRI].push_back(i0);
				m_Indices[TRI].push_back(i1);
			}
		}
	}
	m_Indices[QUAD].clear();
}

void CMesh::Translate(XYZ Vector)
{
	vector<XYZ>::iterator itNode;
	for (itNode = m_Nodes.begin(); itNode != m_Nodes.end(); ++itNode)
	{
		*itNode += Vector;
	}
}

void CMesh::Rotate(WXYZ Rotation, XYZ Origin)
{
	vector<XYZ>::iterator itNode;
	for (itNode = m_Nodes.begin(); itNode != m_Nodes.end(); ++itNode)
	{
		(*itNode) = Rotation * (*itNode-Origin) + Origin;
	}	
}

void CMesh::FlipNormals()
{
	list<int>::iterator itIter;
	list<int>::iterator it0, it1, it2, it3;
	int i0, i1, i2, i3;
	for (itIter = m_Indices[QUAD].begin(); itIter != m_Indices[QUAD].end(); )
	{
		it0 = itIter++;
		it1 = itIter++;
		it2 = itIter++;
		it3 = itIter++;
		i0 = *it0;
		i1 = *it1;
		i2 = *it2;
		i3 = *it3;
		*it0 = i3;
		*it1 = i2;
		*it2 = i1;
		*it3 = i0;
	}
	for (itIter = m_Indices[TRI].begin(); itIter != m_Indices[TRI].end(); )
	{
		it0 = itIter++;
		it1 = itIter++;
		it2 = itIter++;
		i0 = *it0;
		i1 = *it1;
		i2 = *it2;
		*it0 = i2;
		*it1 = i1;
		*it2 = i0;
	}
}

void CMesh::MeshClosedLoop(const XYZ &Normal, const vector<int> &ClosedLoopVector, bool bQuality)
{
	// A more efficient algorithm is described in http://www.cs.umd.edu/~mount/754/Lects/754lects.pdf
	const double TOL = 1e-10;

	list<pair<int, int> > EdgeStack;
	list<pair<int, int> >::iterator itEdge;

	MergeNodes();

	int i;
	for (i = 0; i < NUM_ELEMENT_TYPES; ++i)
	{
		m_Indices[i].clear();
	}

	if (m_Nodes.size() < 3)
		return;

	list<int> &TriangleIndices = m_Indices[TRI];
	// Generate a list of triangle normals
	list<PLANE> TrianglePlanes;
	PLANE Plane;
	XYZ P1, P2, P3, P;

	// Find the first set of 3 points that are not collinear to create a triangle
	for (i=0; i+2<(int)m_Nodes.size(); ++i)
	{
		P1 = m_Nodes[i+0];
		P2 = m_Nodes[i+1];
		P3 = m_Nodes[i+2];
		Plane.Normal = CrossProduct(P2 - P1, P3 - P1);
		if (GetLengthSquared(Plane.Normal) > TOL*TOL)
		{
			break;
		}
	}
	if (i+2 == (int)m_Nodes.size())
	{
		TGERROR("Unable to create convex hull, all points are collinear");
		assert(false);
		return;
	}

	TriangleIndices.push_back(i+0);
	TriangleIndices.push_back(i+1);
	TriangleIndices.push_back(i+2);
	Normalise(Plane.Normal);
	Plane.d = DotProduct(Plane.Normal, P1);
	TrianglePlanes.push_back(Plane);

	TriangleIndices.push_back(i+0);
	TriangleIndices.push_back(i+2);
	TriangleIndices.push_back(i+1);
	Plane.Normal = -Plane.Normal;
	Plane.d = -Plane.d;
	TrianglePlanes.push_back(Plane);

	list<int>::iterator itIter2, itTriStart;
	list<PLANE>::iterator itPlane;
	int i1, i2, i3;
	bool bFacesAdded;
	// This loop will be iterated until no more faces are added, this is necessary when
	// we have more than 3 nodes that are coplanar. If not some nodes will be missed.
	do
	{
		bFacesAdded = false;
		for (i = 0; i < (int)m_Nodes.size(); ++i)
		{
			P = m_Nodes[i];
			// Delete faces that this vertex can see
			for (itIter2 = TriangleIndices.begin(), itPlane = TrianglePlanes.begin(); itIter2 != TriangleIndices.end(); )
			{
				itTriStart = itIter2;
				i1 = *(itIter2++);
				i2 = *(itIter2++);
				i3 = *(itIter2++);
				// If the vertex can see the plane (with a tolerance) we don't want to remove
				// triangles that are in the same place as the vertex or the algorithm will fail.
				if (DotProduct(itPlane->Normal, P) > itPlane->d+TOL)
				{
					// If the edge already exist in the edge stack then it should cancel with it
					// (remove the edge rather than adding it) otherwise add the edge.
					AddOrCancel(EdgeStack, pair<int, int>(i1, i2));
					AddOrCancel(EdgeStack, pair<int, int>(i2, i3));
					AddOrCancel(EdgeStack, pair<int, int>(i3, i1));
					// Delete the triangle
					itIter2 = TriangleIndices.erase(itTriStart, itIter2);
					itPlane = TrianglePlanes.erase(itPlane);
				}
				else
					++itPlane;
			}
			// Create new triangles and calculate the planes of the new triangle.
			// Not only is it an optimisation to calculate the triangle planes only once
			// it is also for consitency.
			for (itEdge = EdgeStack.begin(); itEdge != EdgeStack.end(); ++itEdge)
			{
				i1 = itEdge->first;
				i2 = itEdge->second;
				i3 = i;
				TriangleIndices.push_back(i1);
				TriangleIndices.push_back(i2);
				TriangleIndices.push_back(i3);
				P1 = m_Nodes[i1];
				P2 = m_Nodes[i2];
				P3 = m_Nodes[i3];
				Plane.Normal = CrossProduct(P2 - P1, P3 - P1);
				Normalise(Plane.Normal);
				Plane.d = DotProduct(Plane.Normal, P1);
				TrianglePlanes.push_back(Plane);
				bFacesAdded = true;
			}
			EdgeStack.clear();
		}
	} while (bFacesAdded);
}

void CMesh::BuildGrid(XYZ Min, XYZ Max, int iNumX, int iNumY, int iNumZ)
{
	XYZ P;
	int i, j, k;
	for (i=0; i<iNumX; ++i)
	{
		for (j=0; j<iNumY; ++j)
		{
			for (k=0; k<iNumZ; ++k)
			{
				P.x = i/double(iNumX-1);
				P.y = j/double(iNumY-1);
				P.z = k/double(iNumZ-1);
				P = Min + (Max-Min)*P;
				m_Nodes.push_back(P);
				// Create hex elements out of the grid
				if (i < iNumX-1 && j < iNumY-1 && k < iNumZ-1)
				{
					m_Indices[HEX].push_back((k+0) + (j+0)*iNumZ + (i+0)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+0) + (j+0)*iNumZ + (i+1)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+0) + (j+1)*iNumZ + (i+1)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+0) + (j+1)*iNumZ + (i+0)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+1) + (j+0)*iNumZ + (i+0)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+1) + (j+0)*iNumZ + (i+1)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+1) + (j+1)*iNumZ + (i+1)*iNumZ*iNumY);
					m_Indices[HEX].push_back((k+1) + (j+1)*iNumZ + (i+0)*iNumZ*iNumY);
				}
			}
		}
	}
}

void CMesh::BuildGrid(XYZ Min, XYZ Max, double dPointsPerUnit)
{
	int iNumX, iNumY, iNumZ;
	iNumX = Round((Max.x-Min.x)*dPointsPerUnit);
	iNumY = Round((Max.y-Min.y)*dPointsPerUnit);
	iNumZ = Round((Max.z-Min.z)*dPointsPerUnit);
	BuildGrid(Min, Max, iNumX, iNumY, iNumZ);
}

void CMesh::WriteBinaryXYZ(ostream &Output, XYZ Vector)
{
	float val;
	int i;
	for (i=0; i<3; ++i)
	{
		val = (float)Vector[i];
		Output.write((char*)&val, 4);
	}
}

bool CMesh::SaveToSTL(string Filename, bool bBinary) const
{
	AddExtensionIfMissing(Filename, ".stl");

	int i;
	for (i = 0; i < NUM_ELEMENT_TYPES; ++i)
	{
		if (i != TRI && !m_Indices[i].empty())
		{
			CMesh TriMesh(*this);
			TriMesh.ConvertToTriangleMesh();
			return TriMesh.SaveToSTL(Filename, bBinary);
		}
	}

	ofstream Output;
	if (bBinary)
		Output.open(Filename.c_str(), ios::out|ios::binary);
	else
		Output.open(Filename.c_str(), ios::out);
	if (!Output)
		return false;

	const list<int> &TriangleIndices = m_Indices[TRI];

	if (bBinary)
	{
		char szHeader[80];
		int iNumFacets = TriangleIndices.size()/3;
		strncpy(szHeader, Filename.c_str(), 80);
		Output.write(szHeader, 80);
		Output.write((char*)&iNumFacets, 4);
	}
	else
		Output << "solid " << Filename << endl;

	XYZ T1, T2, T3, Normal;
	list<int>::const_iterator itIndex;
	short int Padding = 0;
	for (itIndex = TriangleIndices.begin(); itIndex != TriangleIndices.end(); )
	{
		T1 = m_Nodes[*(itIndex++)];
		T2 = m_Nodes[*(itIndex++)];
		T3 = m_Nodes[*(itIndex++)];

		Normal = CrossProduct(T2-T1, T3-T1);
		Normalise(Normal);

		if (bBinary)
		{
			WriteBinaryXYZ(Output, Normal);
			WriteBinaryXYZ(Output, T1);
			WriteBinaryXYZ(Output, T2);
			WriteBinaryXYZ(Output, T3);
			Output.write((char*)&Padding, 2);
		}
		else
		{
			Output << " facet normal " << Normal.x << " " << Normal.y << " " << Normal.z << endl;
			Output << "  outer loop" << endl;
			Output << "   vertex " << T1.x << " " << T1.y << " " << T1.z << endl;
			Output << "   vertex " << T2.x << " " << T2.y << " " << T2.z << endl;
			Output << "   vertex " << T3.x << " " << T3.y << " " << T3.z << endl;
			Output << "  endloop" << endl;
			Output << " endfacet" << endl;
		}
	}

	if (!bBinary)
		Output << "endsolid " << Filename << endl;

	Output.close();

	return true;
}

bool CMesh::SaveToSMESH(string Filename) const
{
	AddExtensionIfMissing(Filename, ".smesh");

	ofstream Output(Filename.c_str());
	if (!Output)
		return false;
	Output << "# node count, 3 dim, no attribute, no boundary marker" << endl;
	Output << m_Nodes.size() << " 3 0 0" << endl;

	Output << "# node index, node coordinates" << endl;
	int iNodeIndex;
	vector<XYZ>::const_iterator itNode;
	for (itNode = m_Nodes.begin(), iNodeIndex=0; itNode != m_Nodes.end(); ++itNode, ++iNodeIndex)
	{
		Output << iNodeIndex << " " << itNode->x << " " << itNode->y << " " << itNode->z << endl;
	}

	int iNumTriangles = m_Indices[TRI].size()/3;
	int iNumQuads = m_Indices[QUAD].size()/4;

	Output << "# facet count, no boundary marker" << endl;
	Output << iNumTriangles+iNumQuads << " 0" << endl;
	Output << "# facets" << endl;
	list<int>::const_iterator itIndex;
	int iNumNodesPerElement;
	int i, j;
	ELEMENT_TYPE ElemType;
	for (j=0; j<2; ++j)
	{
		if (j==0)
			ElemType = QUAD;
		else
			ElemType = TRI;
		iNumNodesPerElement = GetNumNodes((ELEMENT_TYPE)j);
		int iNumElements = m_Indices[ElemType].size()/iNumNodesPerElement;
		for (itIndex = m_Indices[ElemType].begin(), i=0; itIndex != m_Indices[ElemType].end(); ++i)
		{
			Output << iNumNodesPerElement << " ";
			for (int k=0; k<iNumNodesPerElement; ++k, ++itIndex)
			{
				if (k>0)
					Output << " ";
				Output << (*itIndex);
			}
			Output << endl;
		}
	}
	return true;
}