#pragma once

#include "ILMBase.h"

class Mesh;
class Halfedge;

using namespace std;

namespace FloorplanSpace
{	
	Vec2i ToVec2i(pair<int, int> &p);
	pair<int, int> ToPair(Vec2i &v);	

	//a tile template
	const int TemplateDim = 5;
	class Template
	{
	public:
		bool v[TemplateDim][TemplateDim];  //the mask
		Template()
		{
			memset(v, false, TemplateDim * TemplateDim);
		}
		Template(bool value)
		{
			memset(v, value, TemplateDim * TemplateDim);
		}

		void Reset()
		{
			memset(v, false, TemplateDim * TemplateDim);
		}

		//align the mask to the lower-left corner
		void Align();

		//return the BB (min-max of x and y) of mask, inclusive
		Vec4i BB();
	};

	//a solved tile 
	class Tile
	{
	public:
		vector<Vec2i> cells;  //active cells (x,y)		
	};
	
	//templates and sub-tiles enumeration
	//all_templates: all ground-truth templates
	//root and leafs: the 2-level sub-tiles
	bool SubTileEnumeration(bool addition);  //way#1	
	bool SubTileEnumeration2();  //way#2

	//decompose g_all_templates into roots and leafs computationally
	bool ComputeRootLeafs();	

	//verify results of ComputeRootLeafs()
	bool VerifyComputeRootLeafs();

	//a floorplan manager
	class Floorplan
	{
	public:
		bool *m_domain;
		int m_domain_W;  //width and height of domain
		int m_domain_H;  

		vector<Tile> m_tiles;  //solved tiles
		int m_draw_tile;  //test: draw only a specific tile

		//test:
		vector<Tile> m_template_tiles;
		int m_draw_template_tile;

		vector<Tile> m_root_tiles;
		int m_draw_root_tile;
		vector<Tile> m_leaf_tiles;
		int m_draw_leaf_tile;

		Floorplan()
		{
			m_domain = NULL;
			m_domain_W = 0;
			m_domain_H = 0;		
			m_draw_template_tile = -1;
			m_draw_tile = -1;
			m_draw_root_tile = -1;
			m_draw_leaf_tile = -1;
		}
		~Floorplan()
		{
			Reset();
		}

		void Reset()
		{
			if (m_domain)
			{
				delete m_domain;
				m_domain = NULL;
			}
			m_tiles.clear();
		}

		//initialize the floorplan problem domain
		//width and height: specify a BB of a domain
		//off_cells: disabled cells in the domain
		bool Init(int width, int height, vector<Vec2i> &off_cells);

		//solve the floorplan IP problem (by forest)
		bool Solve();
		bool Solve2();  //use multi-root tile sets
		//solve by baseline placement-based method
		bool SolvePlacements();

		//rectangle-tiles solve tests:
		bool SolveRectangles();
		
		//test negative leafs
		bool SolveNegativeLeafs();

		//test20200614:
		bool SolveNonConvex();

		//main draw function
		void Draw();

	private:

		//pick_x,pick_y,w,h: parameter for a pickup matrix
		void SetCamera(bool set_pickup_matrix, double *pick_x = NULL, double *pick_y = NULL, double *w = NULL, double *h = NULL);
	};

	//ToFloorplan() solution
	class FloorplanSolution
	{
	public:
		vector<vector<bool>> edge_state_vars;  //raw Booleans of edge state vars
		vector<Vec2f> vertex_poss;  //solved vertex positions for vertices
	};

	//Floorplan IP manager

	class FloorplanIP
	{
	public:

		vector<Vec3f> m_ori_poss;  //original vertex poss (for the latest mesh's vertices)
		vector<FloorplanSolution> m_solutions;  //solved solutions

		vector<tuple<Vec2f, int, bool>> m_grid_points;  //pos, #cluster (-1 for "merged"), "boundary" flag
		vector<vector<Vec2f>> m_solved_faces;  //solved polygons

		FloorplanIP()
		{
			Reset();
		}
		void Reset()
		{
			m_ori_poss.clear();
			m_solutions.clear();
			m_solved_faces.clear();
		}

		//solve quad mesh to floor plan IP 
		bool SolveOld(Mesh *mesh, int solutions_to_get);  //line equation based
		bool Solve(Mesh* mesh, int solutions_to_get);  //pure edge direction assignments
		bool Solve2(Mesh* mesh, int solutions_to_get);  //edge direction assignment + poss and line equations
		
		//solve vertex poss given edge angles
		//edge_constraints: edges, preferred 2D angles (0~2PI angle from postive-x CCW), and "hard"-flags
		bool SolvePoss(Mesh* mesh, vector<tuple<Halfedge*, float, bool>>& edge_constraints);
			
		//apply a solution 
		//sol_index: which one in m_solutions?
		//ratio: interpolation of vertex poss? (1 = solution's vertex poss, 0 = ori poss)
		bool Apply(Mesh *mesh, int sol_index, float ratio);

		void Draw();

		////test direct meshing IP

		//solve the regular grids given a mesh's boundary
		bool SolveGrids(Mesh* mesh, float step_size);

		//solve connectivity given grid_points
		//grid_points:  {pos, #cluster, "boundary" flag}
		bool SolveMesh(vector<tuple<Vec2f, int, bool>>& grid_points, float ideal_edge_len);		
	};	
}