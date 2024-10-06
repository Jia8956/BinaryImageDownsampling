#include <windows.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <omp.h>  //OpenMP
#include "glut.h"
#include "gurobi_c++.h"
#include "Basic.h"
#include "Floorplan.h"
#include "Camera.h"
#include "Mesh.h"
#include "Draw.h"

#include "ceres/ceres.h"
//using ceres::AutoDiffCostFunction;
//using ceres::CostFunction;
//using ceres::Problem;
//using ceres::Solver;
//using ceres::Solve;

using namespace FloorplanSpace;

void AddText(char* fmt, ...);
void AddTextImportant(char* fmt, ...);
void AddTextCritical(char* fmt, ...);

extern Camera *g_camera;
extern Vec4i g_viewport;
extern Floorplan g_floorplan;
extern float g_floorplan_boundary_vertex_dist;
extern int g_Gurobi_time;
extern float g_draw_vertex_size;
extern float g_draw_edge_size;

vector<Template> g_all_templates;

//test: single-tree forest
Template g_root;
vector<Template> g_leafs;

//a multi-tree forest
class Forest
{
public:
	vector<Template> roots;
	vector<Template> leafs;
	vector<vector<int>> links;  //for every leaf, a list of indices of its roots

	void Reset()
	{
		roots.clear();
		leafs.clear();
		links.clear();
	}
};
Forest g_forest;

Vec2i FloorplanSpace::ToVec2i(pair<int, int> &p)
{
	return Vec2i(p.first, p.second);
}
pair<int, int> FloorplanSpace::ToPair(Vec2i &v)
{
	return make_pair(v.x, v.y);
}

void FloorplanSpace::Template::Align()
{
	Vec4i bb = BB();
	Vec2i offset(bb[0], bb[2]);

	bool V[TemplateDim][TemplateDim] = { false };  //new mask
	for (int x = 0; x < TemplateDim; x++)
	{
		for (int y = 0; y < TemplateDim; y++)
		{
			if (v[x][y])
			{
				V[x - offset.x][y - offset.y] = true;
			}
		}
	}

	memcpy(v, V, sizeof(bool) * TemplateDim * TemplateDim);
}

Vec4i FloorplanSpace::Template::BB()
{
	Vec4i bb(INT_MAX, -INT_MAX, INT_MAX, -INT_MAX);

	for (int x = 0; x < TemplateDim; x++)
	{
		for (int y = 0; y < TemplateDim; y++)
		{
			if (v[x][y])
			{
				if (x < bb[0])
					bb[0] = x;
				if (x > bb[1])
					bb[1] = x;
				if (y < bb[2])
					bb[2] = y;
				if (y > bb[3])
					bb[3] = y;
			}
		}
	}

	return bb;
}

//utility function to enuemrate all possible combinations of leaf candidates to form a tiling mask, ex-root
bool EnumerateCombinations(Template &mask, Template &root, vector<Template> &leafs, vector<vector<int>> &combinations, bool just_one)
{
	if (leafs.size() == 0)
		return true;  //nope

	GRBEnv env = GRBEnv();

	//exhaustive enumeration of all possible combinations
	vector< vector<bool> > existing_solutions;
	while (true)
	{
		GRBModel model = GRBModel(env);

		vector<GRBVar> vars;
		for (int i = 0; i < leafs.size(); i++)
		{
			vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		model.update();

		////constraints:

		//size constraint?
		/*{
			GRBLinExpr sum;
			for (int i = 0; i < vars.size(); i++)
			{
				sum += vars[i];
			}
			model.addConstr(sum <= 5);
		}*/

		//for every ex-root cell, it is covered by exactly one template
		for (int x = 0; x < TemplateDim; x++)
		{
			for (int y = 0; y < TemplateDim; y++)
			{
				if (root.v[x][y])
					continue;  //only ex-root cells

				GRBLinExpr eq;
				for (int i = 0; i < leafs.size(); i++)
				{
					if (leafs[i].v[x][y])
					{
						eq += vars[i];
					}
				}

				if (mask.v[x][y])
				{					
					model.addConstr(eq == 1);
				}
				else
				{
					model.addConstr(eq == 0);
				}
			}
		}		

		//forbid existing combinations to appear again
		for (int i = 0; i < existing_solutions.size(); i++)
		{
			vector<bool> &existing_solution = existing_solutions[i];

			GRBLinExpr eq;
			for (int j = 0; j < vars.size(); j++)
			{
				if (existing_solution[j])
					eq += vars[j];
				else
					eq += 1 - vars[j];
			}
			model.addConstr(eq <= existing_solution.size() - 1);
		}
		
		model.getEnv().set(GRB_IntParam_OutputFlag, false); //slient
		model.getEnv().set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF);  //avoid strange "Max constraint violation exceeds tolerance" problem
		model.optimize();
		int status = model.get(GRB_IntAttr_Status);
		if (status == 3)
		{
			//infeasible			
			break;
		}
		else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
		{
			//some other failure
			AddTextCritical("[EnumerateCombinations] optimize() failed! status:%d", status);
			return false;
		}

		//get results
		vector<bool> solution;
		for (int i = 0; i < vars.size(); i++)
		{
			bool present = round(vars[i].get(GRB_DoubleAttr_X));
			solution.push_back(present);			
		}

		existing_solutions.push_back(solution);	

		//just one solution?
		if (just_one)
			break;
	}	

	//turn to output format
	combinations.clear();
	for (int i = 0; i < existing_solutions.size(); i++)
	{
		vector<int> combination;
		for (int j = 0; j < existing_solutions[i].size(); j++)
		{
			if (existing_solutions[i][j])
				combination.push_back(j);
		}
		combinations.push_back(combination);
	}

	return true;
}

//utility function to enuemrate all possible combinations of leaf candidate that led to tiling masks not of forbidden ones, ex-root
bool EnumerateCombinations2(Template &root, vector<Template> &leafs, vector<Template> &forbiddens, vector<vector<int>> &combinations, bool just_one)
{
	if (leafs.size() == 0)
		return true;  //nope

	GRBEnv env = GRBEnv();

	//exhaustive enumeration of all possible combinations (as leaf indices)
	vector< vector<bool> > existing_solutions;
	while (true)
	{
		GRBModel model = GRBModel(env);

		vector<GRBVar> leaf_vars;
		for (int i = 0; i < leafs.size(); i++)
		{
			leaf_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		//a Boolean var for every ex-root cell
		vector<GRBVar> cell_vars;
		for (int x = 0; x < TemplateDim; x++)
		{
			for (int y = 0; y < TemplateDim; y++)
			{
				if (!root.v[x][y])
				{
					cell_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
				}
			}
		}

		model.update();

		////constraints:

		//for every ex-root cell, it is covered by up one template, and denote its coverage by the corresponding cell_var
		int cell_index = 0;
		for (int x = 0; x < TemplateDim; x++)
		{
			for (int y = 0; y < TemplateDim; y++)
			{
				if (root.v[x][y])
					continue;  //only ex-root cells

				GRBLinExpr eq;
				for (int i = 0; i < leafs.size(); i++)
				{
					if (leafs[i].v[x][y])
					{
						eq += leaf_vars[i];
					}
				}

				GRBVar cell_var = cell_vars[cell_index];
				cell_index++;

				model.addConstr(eq == cell_var);  //so eq naturally <= 1
			}
		}

		//forbid forbidden combinations of cell vars
		for (int i = 0; i < forbiddens.size(); i++)
		{
			Template &f = forbiddens[i];

			GRBLinExpr eq;
			int cell_index = 0;
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					if (root.v[x][y])
						continue;  //only ex-root cells

					GRBVar cell_var = cell_vars[cell_index];
					cell_index++;

					if (f.v[x][y])
					{
						eq += cell_var;
					}
					else
					{
						eq += (1 - cell_var);
					}
				}
			}
			model.addConstr(eq <= cell_index - 1);
		}

		//forbid existing combinations to appear again
		for (int i = 0; i < existing_solutions.size(); i++)
		{
			vector<bool> &existing_solution = existing_solutions[i];

			GRBLinExpr eq;
			for (int j = 0; j < leaf_vars.size(); j++)
			{
				if (existing_solution[j])
					eq += leaf_vars[j];
				else
					eq += 1 - leaf_vars[j];
			}
			model.addConstr(eq <= existing_solution.size() - 1);
		}

		model.getEnv().set(GRB_IntParam_OutputFlag, false); //slient
		model.getEnv().set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF);  //avoid strange "Max constraint violation exceeds tolerance" problem
		model.optimize();
		int status = model.get(GRB_IntAttr_Status);
		if (status == 3)
		{
			//infeasible			
			break;
		}
		else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
		{
			//some other failure
			AddTextCritical("[EnumerateCombinations2] optimize() failed! status:%d", status);
			return false;
		}

		//get results
		vector<bool> solution;
		for (int i = 0; i < leaf_vars.size(); i++)
		{
			bool present = round(leaf_vars[i].get(GRB_DoubleAttr_X));
			solution.push_back(present);
		}

		existing_solutions.push_back(solution);
		
		if (just_one)
			break;
	}

	//turn to output format
	combinations.clear();
	for (int i = 0; i < existing_solutions.size(); i++)
	{
		vector<int> combination;
		for (int j = 0; j < existing_solutions[i].size(); j++)
		{
			if (existing_solutions[i][j])
				combination.push_back(j);
		}
		combinations.push_back(combination);
	}

	return true;
}

bool FloorplanSpace::VerifyComputeRootLeafs()
{
	//assume g_templates is loaded

	//load roots from file
	vector<Template> roots;
	{
		ifstream infile("ComputeRootLeafs2_roots");
		if (!infile.fail())
		{	
			string current_line;
			getline(infile, current_line);

			vector<string> tokens;
			SplitString(current_line, " ", tokens);
			for (int jj = 0; jj < tokens.size(); jj++)
			{
				int index = 0;
				sscanf(tokens[jj].c_str(), "%d", &index);

				roots.push_back(g_all_templates[index]);
			}
		}
		infile.close();

		cout << "[VerifyComputeRootLeafs] roots:" << roots.size() << endl;
	}

	//load leaf candidates from file
	vector<Template> leaf_candidates;
	{
		ifstream infile("ComputeRootLeafs2_leafs");
		if (!infile.fail())
		{
			string current_line;
			while (getline(infile, current_line))
			{
				if (current_line.empty())
					continue;

				//load a leaf mask
				vector<bool> booleans;
				for (int i = 0; i < current_line.length(); i++)
				{					
					if (current_line[i] == '1')
						booleans.push_back(true);
					else if (current_line[i] == '0')
						booleans.push_back(false);					
				}

				Template l(false);
				int index = 0;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						l.v[x][y] = booleans[index];
						index++;
					}
				}
				leaf_candidates.push_back(l);
			}
		}
		infile.close();

		cout << "[VerifyComputeRootLeafs] leafs:" << leaf_candidates.size() << endl;
	}

	//load the computed results (leafs of every roots)
	vector<vector<int>> root_leafs_all;
	{
		ifstream infile("ComputeRootLeafs2_root_leafs");
		if (!infile.fail())
		{
			string current_line;
			while (getline(infile, current_line))
			{
				if (current_line.empty())
					continue;
				else
				{
					vector<int> root_leafs;
					vector<string> tokens;
					SplitString(current_line, " ", tokens);
					for (int i = 0; i < tokens.size(); i++)
					{
						int index = -1;
						sscanf(tokens[i].c_str(), "%d", &index);
						root_leafs.push_back(index);
					}
					root_leafs_all.push_back(root_leafs);
				}
			}
		}
		infile.close();

		cout << "[VerifyComputeRootLeafs] root_leafs:" << root_leafs_all.size() << endl;
	}

	////do some tests:

	//for every g_all_template, it should be achieable by some combination of root + leafs
	if (true)
	{
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			int ok_root = -1;
			for (int j = 0; j < roots.size(); j++)
			{
				//leafs of this root
				vector<Template> leafs;
				for (int k = 0; k < root_leafs_all[j].size(); k++)
				{
					leafs.push_back(leaf_candidates[root_leafs_all[j][k]]);
				}

				vector<vector<int>> combs;
				EnumerateCombinations(g_all_templates[i], roots[j], leafs, combs, true/*just one*/);
				if (combs.size() > 0)
				{
					ok_root = j;
					break;
				}
			}
			if (ok_root == -1)
				cout << "[Verify valid] oops! template#" << i << " cannot be achieved?" << endl;
		}
	}

	//no non-template mask can be achieved
	if (true)
	{
		for (int i = 0; i < roots.size(); i++)
		{
			//leafs of this root
			vector<Template> leafs;
			vector<int> leaf_indices;
			for (int j = 0; j < root_leafs_all[i].size(); j++)
			{
				leafs.push_back(leaf_candidates[root_leafs_all[i][j]]);
				leaf_indices.push_back(root_leafs_all[i][j]);
			}

			vector<vector<int>> combs;
			EnumerateCombinations2(roots[i], leafs, g_all_templates, combs, false/*just one*/);
			if (combs.size() > 0)
			{
				cout << "[Verify invalid] oops! root#" << i << " has invalid solutions?" << combs.size() << endl;
			}
		}
	}

	return true;
}

bool FloorplanSpace::ComputeRootLeafs()
{	
	//collect all possible (unique) leaf candidates first
	vector<Template> leaf_candidates;
	map<string, bool> leaf_candidates_map;
	for (int i = 0; i < g_all_templates.size(); i++)
	{
		Template &root = g_all_templates[i];

		for (int j = 0; j < g_all_templates.size(); j++)
		{
			Template &target = g_all_templates[j];

			bool outside = false;  //is the root outside the target?
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					if (root.v[x][y] && !target.v[x][y])
					{
						outside = true;
						break;
					}
				}
			}
			if (!outside)
			{
				Template residual(false);  //target - root
				int cell_count = 0;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (target.v[x][y] && !root.v[x][y])
						{
							residual.v[x][y] = true;
							cell_count++;
						}
					}
				}

				//skip all-empty leaf 
				if (cell_count == 0)
					continue;

				//to string id
				string id;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (residual.v[x][y])
							id += '1';
						else
							id += '0';
					}
				}

				//add to leaf_candidates
				if (leaf_candidates_map.count(id) == 0)
				{	
					leaf_candidates.push_back(residual);
					leaf_candidates_map[id] = true;
				}				
			}
		}
	}
	cout << "[ComputeRoofLeafs2] leaf_candidates:" << leaf_candidates.size() << endl;	

	vector< vector<int> > root_leaf_candidates_all;  //compatible leaf candidates (as indices in leaf_candidates) for every potential root
	vector< map<int, vector<vector<int>>> > valid_combinations_all;  //valid leaf-combinations of every root. key: covering template index, valud: lists of combs
	vector< vector<vector<int>> > invalid_combinations_all;  //invalid leaf-combinations of every root
	//directly load valid_combinations from file?
	ifstream infile("ComputeRootLeafs2_combinations");
	if (!infile.fail())
	{
		//place holders:
		root_leaf_candidates_all.resize(g_all_templates.size());
		valid_combinations_all.resize(g_all_templates.size());
		invalid_combinations_all.resize(g_all_templates.size());

		int cur_root = -1;
		int cur_valid_target = -1;

		string current_line;
		while (getline(infile, current_line))
		{
			if (current_line.empty())
				continue;

			if (current_line[0] == 'L')
			{
				vector<string> tokens;
				SplitString(current_line, " ", tokens);
				for (int jj = 1; jj < tokens.size(); jj++)
				{
					int l = 0;
					sscanf(tokens[jj].c_str(), "%d", &l);
					root_leaf_candidates_all[cur_root].push_back(l);
				}
			}
			else if (current_line[0] == 'R')
			{
				sscanf(current_line.c_str(), "R %d", &cur_root);
			}
			else if (current_line[0] == 'I')
			{
				//load invalid configurations
				vector<string> tokens;
				SplitString(current_line, " ", tokens);
				for (int jj = 1; jj < tokens.size(); jj++)
				{
					vector<int> comb;

					vector<string> ts;
					SplitString(tokens[jj].c_str(), ",", ts);
					for (int kk = 0; kk < ts.size(); kk++)
					{
						int i = 0;
						sscanf(ts[kk].c_str(), "%d", &i);
						comb.push_back(i);
					}

					invalid_combinations_all[cur_root].push_back(comb);
				}
			}
			else if (current_line[0] == 'V')
			{
				sscanf(current_line.c_str(), "V %d", &cur_valid_target);
			}
			else
			{
				//load a set of combinations for the current valid target
				vector< vector<int> > combs;

				vector<string> tokens;
				SplitString(current_line, " ", tokens);
				for (int jj = 0; jj < tokens.size(); jj++)
				{
					vector<int> comb;

					vector<string> ts;
					SplitString(tokens[jj].c_str(), ",", ts);
					for (int kk = 0; kk < ts.size(); kk++)
					{
						int i = 0;
						sscanf(ts[kk].c_str(), "%d", &i);
						comb.push_back(i);
					}

					valid_combinations_all[cur_root][cur_valid_target].push_back(comb);
				}
			}
		}
		infile.close();		
	}	
	else
	{
		//a map of all g_templates
		map<string, bool> templates_map;
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			Template &t = g_all_templates[i];
			string id;
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					if (t.v[x][y])
						id += '1';
					else
						id += '0';
				}
			}
			templates_map[id] = true;
		}		

		for (int i = 0; i < g_all_templates.size(); i++)
		{
			Template &root = g_all_templates[i];

			//first collect leaf candidates that jointly would lead to a valid cell mask (i.e., euqal to a template)
			vector<Template> leaf_candidates_this;
			vector<int> leaf_candidates_this_indices;  //as indices in leaf_candidates
			for (int j = 0; j < leaf_candidates.size(); j++)
			{
				Template &leaf = leaf_candidates[j];

				bool conflict = false;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (root.v[x][y] && leaf.v[x][y])
						{
							conflict = true;
							break;
						}
					}
				}
				if (!conflict)
				{
					//calculate the cell mask of root + leaf
					string id;
					for (int x = 0; x < TemplateDim; x++)
					{
						for (int y = 0; y < TemplateDim; y++)
						{
							if (root.v[x][y] || leaf.v[x][y])
								id += '1';
							else
								id += '0';
						}
					}
					if (templates_map.count(id) > 0)
					{
						leaf_candidates_this.push_back(leaf);						
						leaf_candidates_this_indices.push_back(j);
					}
				}
			}

			root_leaf_candidates_all.push_back(leaf_candidates_this_indices);

			//enumerate invalid combinations
			{
				vector < vector<int> > combs;
				EnumerateCombinations2(root, leaf_candidates_this, g_all_templates, combs, false/*juse_one?*/);
				//translate to indices in leaf_candidates
				vector < vector<int> > invalid_combinations;
				for (int j = 0; j < combs.size(); j++)
				{
					vector<int> indices;
					for (int k = 0; k < combs[j].size(); k++)
					{
						indices.push_back(leaf_candidates_this_indices[combs[j][k]]);
					}
					invalid_combinations.push_back(indices);
				}
				invalid_combinations_all.push_back(invalid_combinations);
			}

			//enumerate valid combinations:
			map<int, vector<vector<int>>> valid_combinations;
			for (int j = 0; j < g_all_templates.size(); j++)
			{
				Template &target = g_all_templates[j];

				//is the target template compatible with the root at all?
				bool compatible = true;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (root.v[x][y] && !target.v[x][y])
						{
							compatible = false;
							break;
						}
					}
				}
				if (!compatible)
					continue;

				//now, collect leafs that are further compatible with the target template
				vector<Template> compatible_leafs;
				vector<int> compatible_leaf_indices;  //as indices in leaf_candidates
				for (int k = 0; k < leaf_candidates_this.size(); k++)
				{
					Template &l = leaf_candidates_this[k];

					bool compatible = true;
					for (int x = 0; x < TemplateDim; x++)
					{
						for (int y = 0; y < TemplateDim; y++)
						{
							if (l.v[x][y] && !target.v[x][y])
							{
								compatible = false;
								break;
							}
						}
					}
					if (compatible)
					{
						compatible_leafs.push_back(l);
						compatible_leaf_indices.push_back(leaf_candidates_this_indices[k]);
					}
				}
				cout << "root#" << i << " " << leaf_candidates_this.size() << " target_template:" << j << " compatible_leafs:" << compatible_leafs.size() << endl;

				//emnumerate all possible combinations of leafs that lead to this template
				vector<vector<int>> combs;
				EnumerateCombinations(target, root, compatible_leafs, combs, false/*just one?*/);

				//now, check if each comb would not lead to invalid cell masks
				for (int k = 0; k < combs.size(); k++)
				{
					vector<Template> leafs;
					for (int q = 0; q < combs[k].size(); q++)
					{
						leafs.push_back(compatible_leafs[combs[k][q]]);
					}

					//try to solve for any cell mask that is not a template
					vector<vector<int>> cs;  //dummy
					EnumerateCombinations2(root, leafs, g_all_templates, cs, true/*just one?*/);
					if (cs.size() == 0)  //the solution should be infeasible
					{
						//NOTE: translate indices back to indices among leaf_candidates!
						vector<int> indices;
						for (int q = 0; q < combs[k].size(); q++)
						{
							indices.push_back(compatible_leaf_indices[combs[k][q]]);
						}
						valid_combinations[j].push_back(indices);
					}
				}
			}

			cout << "root#" << i << " valid_combinations:" << valid_combinations.size() << endl;
			valid_combinations_all.push_back(valid_combinations);
		}

		//save stuffs to files
		{
			ofstream outfile("ComputeRootLeafs2_combinations");
			for (int i = 0; i < g_all_templates.size(); i++)
			{
				outfile << "R " << i << endl;

				outfile << "L ";
				vector<int> &root_leaf_candidates = root_leaf_candidates_all[i];
				for (int j = 0; j < root_leaf_candidates.size(); j++)
				{
					outfile << root_leaf_candidates[j];
					if (j < root_leaf_candidates.size() - 1)
						outfile << " ";
				}
				outfile << endl;

				//invalid combinations:
				outfile << "I ";
				vector<vector<int>> &invalid_combinations = invalid_combinations_all[i];
				for (int j = 0; j < invalid_combinations.size(); j++)
				{
					vector<int> &comb = invalid_combinations[j];
					for (int k = 0; k < comb.size(); k++)
					{
						outfile << comb[k];
						if (k < comb.size() - 1)
							outfile << ",";
					}
					if (j < invalid_combinations.size() - 1)
						outfile << " ";
				}
				outfile << endl;

				//valid combinations (targets and combs)
				map<int, vector<vector<int>>> &valid_combinations = valid_combinations_all[i];
				for (map<int, vector<vector<int>>>::iterator itr = valid_combinations.begin(); itr != valid_combinations.end(); itr++)
				{
					outfile << "V " << (*itr).first << endl;

					vector<vector<int>> &combs = (*itr).second;
					for (int j = 0; j < combs.size(); j++)
					{
						for (int k = 0; k < combs[j].size(); k++)
						{
							outfile << combs[j][k];
							if (k < combs[j].size() - 1)
								outfile << ",";
						}
						if (j < combs.size() - 1)
							outfile << " ";
					}
					outfile << endl;
				}
			}
			outfile.close();
		}
	}	

	//now, solve for a IP problem of finding a  set of roots and root-leafs that together covers every g_template
	vector<int> solved_roots;
	vector<vector<int>> solved_root_leafs;
	bool loaded = true;
	do
	{
		ifstream infile1("ComputeRootLeafs2_roots");
		if (!infile1.fail())
		{
			string current_line;
			
			//the file should contain a single line of many tokens
			getline(infile1, current_line);
			vector<string> tokens;
			SplitString(current_line, " ", tokens);
			for (int i = 0; i < tokens.size(); i++)
			{
				int index = -1;
				sscanf(tokens[i].c_str(), "%d", &index);
				solved_roots.push_back(index);
			}
		}
		else
		{
			loaded = false;
			break;
		}

		ifstream infile2("ComputeRootLeafs2_root_leafs");
		if (!infile2.fail())
		{
			int cur_root = 0;
			string current_line;
			while (getline(infile2, current_line))
			{
				if (current_line.empty())
					continue;

				//'n': nothing
				if (current_line[0] == 'n')
				{
					vector<int> indices;
					solved_root_leafs.push_back(indices);
					continue;
				}

				//one line per root's leafs
				vector<int> indices;
				vector<string> tokens;
				SplitString(current_line, " ", tokens);
				for (int i = 0; i < tokens.size(); i++)
				{
					int index = -1;
					sscanf(tokens[i].c_str(), "%d", &index);
					indices.push_back(index);
				}
				solved_root_leafs.push_back(indices);
			}
		}
		else
		{
			loaded = false;
			break;
		}
	} while (false);
	if(!loaded)
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		vector<GRBVar> root_vars;
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			root_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		//per-root leaf flags
		vector< vector<GRBVar> > root_leaf_vars_all;
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			vector<GRBVar> vars;  //one Boolean flag per leaf candidate
			for (int j = 0; j < leaf_candidates.size(); j++)
			{
				vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
			}

			root_leaf_vars_all.push_back(vars);
		}

		//leaf presence flags
		vector<GRBVar> leaf_vars;
		for (int i = 0; i < leaf_candidates.size(); i++)
		{
			leaf_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		//valid combination flags
		vector< map<int, vector<GRBVar>> > valid_combination_vars_all;
		for (int i = 0; i < valid_combinations_all.size(); i++)
		{
			map<int, vector<GRBVar>> valid_combination_vars;
			for (map<int, vector<vector<int>>>::iterator itr = valid_combinations_all[i].begin(); itr != valid_combinations_all[i].end(); itr++)
			{
				vector<GRBVar> vars;
				vector<vector<int>> &combs = (*itr).second;
				for (int j = 0; j < combs.size(); j++)
				{
					vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
				}
				valid_combination_vars[(*itr).first] = vars;
			}
			valid_combination_vars_all.push_back(valid_combination_vars);
		}

		model.update();

		//obj: primal goal is minimize # of active roots
		GRBLinExpr obj;
		for (int i = 0; i < root_vars.size(); i++)
		{
			//int BIG = leaf_candidates.size();
			//obj += BIG * root_vars[i];
			
			obj += root_vars[i];
		}
		//also minimize # of leafs as a minor goal
		for (int i = 0; i < leaf_vars.size(); i++)
		{	
			obj += leaf_vars[i];
		}
		
		model.setObjective(obj);

		////constraints:	

		//immediately disable impossible root_leafs 
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			vector<int> &ls = root_leaf_candidates_all[i];
			map<int, bool> ls_map;
			for (int j = 0; j < ls.size(); j++)
			{
				ls_map[ls[j]] = true;
			}

			vector<GRBVar> &vars = root_leaf_vars_all[i];
			for (int j = 0; j < vars.size(); j++)
			{
				if (ls_map.count(j) == 0)
					model.addConstr(vars[j] == 0);
			}
		}

		//for every g_template, at least one combination var of root-leafs that cover it is active
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			GRBLinExpr eq;
			for (int j = 0; j < g_all_templates.size(); j++)
			{
				//note: a root must be able to cover itself
				if (j == i)
					eq += root_vars[j];

				map<int, vector<vector<int>>> &valid_combinations = valid_combinations_all[j];
				map<int, vector<GRBVar>> &valid_combination_vars = valid_combination_vars_all[j];

				for (map<int, vector<vector<int>>>::iterator itr = valid_combinations.begin(); itr != valid_combinations.end(); itr++)
				{
					if ((*itr).first == i)
					{
						//any one of the combination_var would do
						vector<GRBVar> &vars = valid_combination_vars[(*itr).first];
						for (int k = 0; k < vars.size(); k++)
						{
							eq += vars[k];
						}
					}
				}
			}
			model.addConstr(eq >= 1);
		}

		//any invalid combinations cannot be present
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			vector<GRBVar> &root_leaf_vars = root_leaf_vars_all[i];

			vector<vector<int>> &combs = invalid_combinations_all[i];
			for (int j = 0; j < combs.size(); j++)
			{
				GRBLinExpr sum;
				for (int k = 0; k < combs[j].size(); k++)
				{
					sum += root_leaf_vars[combs[j][k]];
				}
				model.addConstr(sum <= combs[j].size() - 1);
			}
		}

		//model the valid combination vars
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			vector<GRBVar> &root_leaf_vars = root_leaf_vars_all[i];

			map<int, vector<vector<int>>> &valid_combinations = valid_combinations_all[i];
			map<int, vector<GRBVar>> &valid_combination_vars = valid_combination_vars_all[i];

			for (map<int, vector<vector<int>>>::iterator itr = valid_combinations.begin(); itr != valid_combinations.end(); itr++)
			{
				vector<vector<int>> &combs = (*itr).second;
				vector<GRBVar> &valid_vars = valid_combination_vars[(*itr).first];

				for (int j = 0; j < combs.size(); j++)
				{
					vector<int> &comb = combs[j];
					GRBVar &valid_var = valid_vars[j];

					//the valid_var cannot be true unless all the corresponding root_leaf vars are true 
					GRBLinExpr sum;
					for (int k = 0; k < comb.size(); k++)
					{
						sum += root_leaf_vars[comb[k]];
					}
					int N = comb.size();
					model.addConstr(0 <= sum - N * valid_var);
				}
			}
		}

		//any root_leaf var cannot be true unless its root var is true
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			GRBVar &root_var = root_vars[i];
			vector<GRBVar> &root_leaf_vars = root_leaf_vars_all[i];
			for (int j = 0; j < root_leaf_vars.size(); j++)
			{
				model.addConstr(root_leaf_vars[j] <= root_var);
			}
		}

		//each leaf var is true iff any of the corresponding per-root leaf vars is true
		for (int i = 0; i < leaf_vars.size(); i++)
		{
			GRBVar &l = leaf_vars[i];

			//sum of corresponding per-root leaf vars
			GRBLinExpr sum;
			for (int j = 0; j < g_all_templates.size(); j++)
			{
				sum += root_leaf_vars_all[j][i];
			}

			int N = g_all_templates.size();
			model.addConstr(-N + 1 <= sum - N * l);
			model.addConstr(sum - N * l <= 0);
		}
		
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 200);  //some time limit
		model.optimize();
		int status = model.get(GRB_IntAttr_Status);
		if (status == 3)
		{
			//infeasible
			AddTextCritical("[find root] the problem is infeasible.");
			return false;
		}
		else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
		{
			//some other failure
			AddTextCritical("[find root] optimize() failed! status:%d", status);
			return false;
		}

		//get results

		for (int i = 0; i < root_vars.size(); i++)
		{
			bool present = round(root_vars[i].get(GRB_DoubleAttr_X));
			if (present)
			{
				cout << "[find_root] active root:" << i << "(#" << solved_roots.size() << ")" << endl;

				Template &r = g_all_templates[i];
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (r.v[x][y])
							cout << "1";
						else
							cout << "0";
					}
					cout << endl;
				}
				cout << endl;
				solved_roots.push_back(i);
			}
		}

		int num_leafs = 0;
		for (int i = 0; i < leaf_candidates.size(); i++)
		{
			bool present = round(leaf_vars[i].get(GRB_DoubleAttr_X));
			if (present)
			{
				num_leafs++;
			}
		}
		cout << "[find_root] num_leafs:" << num_leafs << endl;

		int num_root_leafs = 0;
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			vector<int> root_leafs;

			vector<GRBVar> &root_leaf_vars = root_leaf_vars_all[i];
			for (int j = 0; j < root_leaf_vars.size(); j++)
			{
				bool present = round(root_leaf_vars[j].get(GRB_DoubleAttr_X));
				if (present)
				{
					root_leafs.push_back(j);
					num_root_leafs++;
				}
			}

			solved_root_leafs.push_back(root_leafs);
		}
		cout << "[find_roots] num_root_leafs:" << num_root_leafs << endl;

		//save things to files
		
		//save g_templates
		{
			ofstream outfile("ComputeRootLeafs2_templates");
			for (int i = 0; i < g_all_templates.size(); i++)
			{
				//one line per template
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						outfile << g_all_templates[i].v[x][y];
					}
				}
				outfile << endl;
			}
			outfile.close();
		}
		//save leaf candidates
		{
			ofstream outfile("ComputeRootLeafs2_leafs");
			for (int i = 0; i < leaf_candidates.size(); i++)
			{
				//one line per leaf
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						outfile << leaf_candidates[i].v[x][y];
					}
				}
				outfile << endl;
			}
			outfile << endl;
		}
		//save solved roots
		{
			ofstream outfile("ComputeRootLeafs2_roots");
			for (int i = 0; i < solved_roots.size(); i++)
			{
				outfile << solved_roots[i];

				if (i < solved_roots.size() - 1)
					outfile << " ";
			}
			outfile.close();
		}
		//save solved root_leafs
		{
			ofstream outfile("ComputeRootLeafs2_root_leafs");

			map<int, bool> root_map;
			for (int i = 0; i < solved_roots.size(); i++)
			{
				root_map[solved_roots[i]] = true;
			}

			for (int i = 0; i < solved_root_leafs.size(); i++)
			{
				if (root_map.count(i) == 0)
					continue;

				vector<int> &root_leafs = solved_root_leafs[i];
				if (root_leafs.size() > 0)
				{
					for (int j = 0; j < root_leafs.size(); j++)
					{
						outfile << root_leafs[j];

						if (j < root_leafs.size() - 1)
							outfile << " ";
					}
					outfile << endl;
				}
				else
				{
					outfile << "n" << endl;
				}
			}
			outfile.close();
		}
	}

	//export solved_roots and solved_root_leafs to the g_forest
	{
		g_forest.Reset();

		map<int, int> root_map;  //key: root index in g_all_templates, value: index in g_forest.roots
		for (int i = 0; i < solved_roots.size(); i++)
		{
			root_map[solved_roots[i]] = g_forest.roots.size();
			g_forest.roots.push_back(g_all_templates[solved_roots[i]]);
		}
		map<int, vector<int>> root_leafs_map;
		for (int i = 0; i < solved_root_leafs.size(); i++)
		{
			int ri = root_map[solved_roots[i]];
			for (int j = 0; j < solved_root_leafs[i].size(); j++)
			{
				int l = solved_root_leafs[i][j];
				root_leafs_map[l].push_back(ri);
			}
		}
		for (map<int, vector<int>>::iterator itr = root_leafs_map.begin(); itr != root_leafs_map.end(); itr++)
		{
			g_forest.leafs.push_back(leaf_candidates[(*itr).first]);
			g_forest.links.push_back((*itr).second);
		}
	}

	//test: export somethings for show
	{
		g_floorplan.m_root_tiles.clear();
		for (int i = 0; i < solved_roots.size(); i++)
		{			
			Tile t;
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					if (g_all_templates[solved_roots[i]].v[x][y])
						t.cells.push_back(Vec2i(x, y));
				}
			}
			g_floorplan.m_root_tiles.push_back(t);
		}

		g_floorplan.m_leaf_tiles.clear();
		for (int i = 0; i < solved_root_leafs.size(); i++)
		{
			for (int j = 0; j < solved_root_leafs[i].size(); j++)
			{
				Tile t;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (leaf_candidates[solved_root_leafs[i][j]].v[x][y])
							t.cells.push_back(Vec2i(x, y));
					}
				}
				g_floorplan.m_leaf_tiles.push_back(t);
			}
			
		}

	}

	return true;
}

bool FloorplanSpace::SubTileEnumeration2()
{
	////enumerate all possible tile templates

	g_all_templates.clear();

	map<string, tuple<int,int,int,int,int>> templates_map;  //value: ii,cc0, cc1,cc2,cc3

	//6 core shapes:
	for (int ii = 0; ii < 6; ii++)
	{		
		Template core(false);
		//try to align to the center as much as possible?
		if (ii == 0)  //3x3:
		{
			for (int x = 0; x <= 2; x++)
			{
				for (int y = 0; y <= 2; y++)
				{
					core.v[x][y] = true;
				}
			}
			/*for (int x = 1; x <= 3; x++)
			{
				for (int y = 1; y <= 3; y++)
				{
					core.v[x][y] = true;
				}
			}*/
		}
		else if (ii == 1)  //4x3:
		{
			for (int x = 0; x <= 3; x++)
			{
				for (int y = 0; y <= 2; y++)
				//for (int y = 1; y <= 3; y++)
				{
					core.v[x][y] = true;
				}
			}
		}
		else if (ii == 2)  //3x4:
		{
			for (int x = 0; x <= 2; x++)
			//for (int x = 1; x <= 3; x++)
			{
				for (int y = 0; y <= 3; y++)
				{
					core.v[x][y] = true;
				}
			}
		}
		else if (ii == 3)  //4x4:
		{
			for (int x = 0; x <= 3; x++)
			{
				for (int y = 0; y <= 3; y++)
				{
					core.v[x][y] = true;
				}
			}
		}
		else if (ii == 4)  //5x4:
		{
			for (int x = 0; x <= 4; x++)
			{
				for (int y = 0; y <= 3; y++)
				{
					core.v[x][y] = true;
				}
			}
		}
		else if (ii == 5)  //4x5:
		{
			for (int x = 0; x <= 3; x++)
			{
				for (int y = 0; y <= 4; y++)
				{
					core.v[x][y] = true;
				}
			}
		}
		
		Vec4i bb = core.BB();  //min-x, max-y, min-y, max-y
		Vec2i c0(bb[0], bb[2]);  //lower-left
		Vec2i c1(bb[1], bb[2]);  //lower-right
		Vec2i c2(bb[1], bb[3]);  //upper-right
		Vec2i c3(bb[0], bb[3]);  //upper-left

		//try 5 kinds of chisels (none, 1x1, 2x1, 1x2, 2x2) at every corner
		for (int cc0 = 0; cc0 < 5; cc0++)
		{
			for (int cc1 = 0; cc1 < 5; cc1++)
			{
				for (int cc2 = 0; cc2 < 5; cc2++)
				{
					for (int cc3 = 0; cc3 < 5; cc3++)
					{						
						Template t = core;

						//the mask of chisels at the four corners
						map<pair<int, int>, bool> C0, C1, C2, C3;	
						map<pair<int, int>, bool> A0, A1, A2, A3;  //adjacent cells of chisels
						if (cc0 == 0)
						{
							A0[ToPair(c0)] = true;
						}
						else if (cc0 == 1)
						{
							C0[ToPair(c0)] = true;
							A0[ToPair(c0 + Vec2i(1, 0))] = true;
							A0[ToPair(c0 + Vec2i(1, 1))] = true;
							A0[ToPair(c0 + Vec2i(0, 1))] = true;
						}
						else if (cc0 == 2)
						{
							C0[ToPair(c0)] = true;
							C0[ToPair(c0 + Vec2i(1, 0))] = true;
							A0[ToPair(c0 + Vec2i(2, 0))] = true;
							A0[ToPair(c0 + Vec2i(2, 1))] = true;
							A0[ToPair(c0 + Vec2i(1, 1))] = true;
							A0[ToPair(c0 + Vec2i(0, 1))] = true;
						}
						else if (cc0 == 3)
						{
							C0[ToPair(c0)] = true;
							C0[ToPair(c0 + Vec2i(0, 1))] = true;
							A0[ToPair(c0 + Vec2i(1, 0))] = true;
							A0[ToPair(c0 + Vec2i(1, 1))] = true;
							A0[ToPair(c0 + Vec2i(1, 2))] = true;
							A0[ToPair(c0 + Vec2i(0, 2))] = true;
						}
						else if (cc0 == 4)
						{
							C0[ToPair(c0)] = true;
							C0[ToPair(c0 + Vec2i(0, 1))] = true;
							C0[ToPair(c0 + Vec2i(1, 1))] = true;
							C0[ToPair(c0 + Vec2i(1, 0))] = true;
							A0[ToPair(c0 + Vec2i(2, 0))] = true;
							A0[ToPair(c0 + Vec2i(2, 1))] = true;
							A0[ToPair(c0 + Vec2i(2, 2))] = true;
							A0[ToPair(c0 + Vec2i(1, 2))] = true;
							A0[ToPair(c0 + Vec2i(0, 2))] = true;
						}

						if (cc1 == 0)
						{
							A1[ToPair(c1)] = true;
						}
						else if (cc1 == 1)
						{
							C1[ToPair(c1)] = true;
							A1[ToPair(c1 + Vec2i(-1, 0))] = true;
							A1[ToPair(c1 + Vec2i(-1, 1))] = true;
							A1[ToPair(c1 + Vec2i(0, 1))] = true;
						}
						else if (cc1 == 2)
						{
							C1[ToPair(c1)] = true;
							C1[ToPair(c1 + Vec2i(-1, 0))] = true;
							A1[ToPair(c1 + Vec2i(-2, 0))] = true;
							A1[ToPair(c1 + Vec2i(-2, 1))] = true;
							A1[ToPair(c1 + Vec2i(-1, 1))] = true;
							A1[ToPair(c1 + Vec2i(0, 1))] = true;
						}
						else if (cc1 == 3)
						{
							C1[ToPair(c1)] = true;
							C1[ToPair(c1 + Vec2i(0, 1))] = true;
							A1[ToPair(c0 + Vec2i(-1, 0))] = true;
							A1[ToPair(c0 + Vec2i(-1, 1))] = true;
							A1[ToPair(c0 + Vec2i(-1, 2))] = true;
							A1[ToPair(c0 + Vec2i(0, 2))] = true;
						}
						else if (cc1 == 4)
						{
							C1[ToPair(c0)] = true;
							C1[ToPair(c0 + Vec2i(0, 1))] = true;
							C1[ToPair(c0 + Vec2i(-1, 1))] = true;
							C1[ToPair(c0 + Vec2i(-1, 0))] = true;
							A1[ToPair(c0 + Vec2i(-2, 0))] = true;
							A1[ToPair(c0 + Vec2i(-2, 1))] = true;
							A1[ToPair(c0 + Vec2i(-2, 2))] = true;
							A1[ToPair(c0 + Vec2i(-1, 2))] = true;
							A1[ToPair(c0 + Vec2i(0, 2))] = true;
						}

						if (cc2 == 0)
						{
							A2[ToPair(c2)] = true;
						}
						else if (cc2 == 1)
						{
							C2[ToPair(c2)] = true;
							A2[ToPair(c2 + Vec2i(-1, 0))] = true;
							A2[ToPair(c2 + Vec2i(-1, -1))] = true;
							A2[ToPair(c2 + Vec2i(0, -1))] = true;
						}
						else if (cc2 == 2)
						{
							C2[ToPair(c2)] = true;
							C2[ToPair(c2 + Vec2i(-1, 0))] = true;
							A2[ToPair(c2 + Vec2i(-2, 0))] = true;
							A2[ToPair(c2 + Vec2i(-2, -1))] = true;
							A2[ToPair(c2 + Vec2i(-1, -1))] = true;
							A2[ToPair(c2 + Vec2i(0, -1))] = true;
						}
						else if (cc2 == 3)
						{
							C2[ToPair(c2)] = true;
							C2[ToPair(c2 + Vec2i(0, -1))] = true;
							A2[ToPair(c2 + Vec2i(-1, 0))] = true;
							A2[ToPair(c2 + Vec2i(-1, -1))] = true;
							A2[ToPair(c2 + Vec2i(-1, -2))] = true;
							A2[ToPair(c2 + Vec2i(0, -2))] = true;
						}
						else if (cc2 == 4)
						{
							C2[ToPair(c2)] = true;
							C2[ToPair(c2 + Vec2i(0, -1))] = true;
							C2[ToPair(c2 + Vec2i(-1, -1))] = true;
							C2[ToPair(c2 + Vec2i(-1, 0))] = true;
							A2[ToPair(c2 + Vec2i(-2, 0))] = true;
							A2[ToPair(c2 + Vec2i(-2, -1))] = true;
							A2[ToPair(c2 + Vec2i(-2, -2))] = true;
							A2[ToPair(c2 + Vec2i(-1, -2))] = true;
							A2[ToPair(c2 + Vec2i(0, -2))] = true;
						}

						if (cc3 == 0)
						{
							A3[ToPair(c3)] = true;
						}
						else if (cc3 == 1)
						{
							C3[ToPair(c3)] = true;
							A3[ToPair(c3 + Vec2i(1, 0))] = true;
							A3[ToPair(c3 + Vec2i(1, -1))] = true;
							A3[ToPair(c3 + Vec2i(0, -1))] = true;
						}
						else if (cc3 == 2)
						{
							C3[ToPair(c3)] = true;
							C3[ToPair(c3 + Vec2i(1, 0))] = true;
							A3[ToPair(c3 + Vec2i(2, 0))] = true;
							A3[ToPair(c3 + Vec2i(2, -1))] = true;
							A3[ToPair(c3 + Vec2i(1, -1))] = true;
							A3[ToPair(c3 + Vec2i(0, -1))] = true;
						}
						else if (cc3 == 3)
						{
							C3[ToPair(c3)] = true;
							C3[ToPair(c3 + Vec2i(0, -1))] = true;
							A3[ToPair(c3 + Vec2i(1, 0))] = true;
							A3[ToPair(c3 + Vec2i(1, -1))] = true;
							A3[ToPair(c3 + Vec2i(1, -2))] = true;
							A3[ToPair(c3 + Vec2i(0, -2))] = true;
						}
						else if (cc3 == 4)
						{
							C3[ToPair(c3)] = true;
							C3[ToPair(c3 + Vec2i(0, -1))] = true;
							C3[ToPair(c3 + Vec2i(1, -1))] = true;
							C3[ToPair(c3 + Vec2i(1, 0))] = true;
							A3[ToPair(c3 + Vec2i(2, 0))] = true;
							A3[ToPair(c3 + Vec2i(2, -1))] = true;
							A3[ToPair(c3 + Vec2i(2, -2))] = true;
							A3[ToPair(c3 + Vec2i(1, -2))] = true;
							A3[ToPair(c3 + Vec2i(0, -2))] = true;
						}

						//chisel masks shall not overlap
						map<pair<int, int>, int> counts;  //any cell cannot have count>1
						for (int Case = 0; Case < 4; Case++)
						{
							map<pair<int, int>, bool> *C = NULL;
							if (Case == 0)
								C = &C0;
							else if (Case == 1)
								C = &C1;
							else if (Case == 2)
								C = &C2;
							else if (Case == 3)
								C = &C3;

							for (map<pair<int, int>, bool>::iterator itr = (*C).begin(); itr != (*C).end(); itr++)
							{
								counts[(*itr).first]++;
							}
						}
						bool overlap = false;
						for (map<pair<int, int>, int>::iterator itr = counts.begin(); itr != counts.end(); itr++)
						{
							if ((*itr).second > 1)
							{
								overlap = true;
								break;
							}
						}
						if (overlap)
							continue;

						//also check that chisels are not adjacent (i.e., a chisel's adjacent cells cannot overlap with any other chisels' masks)						
						for (int Case = 0; Case < 4; Case++)
						{
							map<pair<int, int>, bool> *A = NULL;
							vector<map<pair<int, int>, bool>*> Cs;
							if (Case == 0)
							{
								A = &A0;
								Cs.push_back(&C1);
								Cs.push_back(&C2);
								Cs.push_back(&C3);
							}
							else if (Case == 1)
							{
								A = &A1;
								Cs.push_back(&C0);
								Cs.push_back(&C2);
								Cs.push_back(&C3);
							}
							else if (Case == 2)
							{
								A = &A2;
								Cs.push_back(&C0);
								Cs.push_back(&C1);
								Cs.push_back(&C3);
							}
							else if (Case == 3)
							{
								A = &A3;
								Cs.push_back(&C0);
								Cs.push_back(&C1);
								Cs.push_back(&C2);
							}

							for (map<pair<int, int>, bool>::iterator itr = (*A).begin(); itr != (*A).end(); itr++)
							{
								for (int q = 0; q < 3; q++)
								{
									if ((*Cs[q]).count((*itr).first) > 0)
									{
										overlap = true;
										break;
									}									
								}
								if (overlap)
									break;
							}
							if (overlap)
								break;
						}
						if (overlap)
						{
							//cout << "overlap!" << ii << " " << cc0 << " " << cc1 << " " << cc2 << " " << cc3 << endl;
							continue;
						}

						//remove chisel cells from the core
						for (int Case = 0; Case < 4; Case++)
						{
							map<pair<int, int>, bool> *C = NULL;
							if (Case == 0)
								C = &C0;
							else if (Case == 1)
								C = &C1;
							else if (Case == 2)
								C = &C2;
							else if (Case == 3)
								C = &C3;

							for (map<pair<int, int>, bool>::iterator itr = (*C).begin(); itr != (*C).end(); itr++)
							{
								t.v[(*itr).first.first][(*itr).first.second] = false;
							}							
						}

						//area range?
						if (true)
						{
							Vec2i area_range(7, 20);  //min-max
							int area = 0;
							for (int x = 0; x < TemplateDim; x++)
							{
								for (int y = 0; y < TemplateDim; y++)
								{
									if (t.v[x][y])
										area++;
								}
							}
							if (area < area_range.x || area > area_range.y)
							{
								//kick out
								//cout << "area!" << ii << " " << cc0 << " " << cc1 << " " << cc2 << " " << cc3 << endl;
								continue;
							}
						}

						//kick out ones w/ "narrow" extrusions?
						if (true)
						{
							//look for active cells w/ off cells on two sides
							bool gotcha = false;
							for (int x = 0; x < TemplateDim; x++)
							{
								for (int y = 0; y < TemplateDim; y++)
								{
									if (t.v[x][y])
									{
										//check left-right
										if (x == 0 && !t.v[x + 1][y])
										{
											gotcha = true;
											break;
										}
										else if (x == TemplateDim - 1 && !t.v[x - 1][y])
										{
											gotcha = true;
											break;
										}
										else if (x > 0 && x < TemplateDim - 1 && !t.v[x - 1][y] && !t.v[x + 1][y])
										{
											gotcha = true;
											break;
										}

										//check up-down
										if (y == 0 && !t.v[x][y + 1])
										{
											gotcha = true;
											break;
										}
										else if (y == TemplateDim - 1 && !t.v[x][y - 1])
										{
											gotcha = true;
											break;
										}
										else if (y > 0 && y < TemplateDim - 1 && !t.v[x][y - 1] && !t.v[x][y + 1])
										{
											gotcha = true;
											break;
										}
									}
									if (gotcha)
										break;
								}
								if (gotcha)
									break;
							}

							if (gotcha)
							{
								//cout << "kick out narrow ii:" << ii << " cases:" << cc0 << "," << cc1 << "," << cc2 << "," << cc3 << endl;
								continue;
							}
						}

						//add a new template?

						string id;
						for (int x = 0; x < TemplateDim; x++)
						{
							for (int y = 0; y < TemplateDim; y++)
							{
								if (t.v[x][y])
									id += '1';
								else
									id += '0';
							}
						}
						if (templates_map.count(id) == 0)
						{
							//cout << "add new ii:" << ii << " cc0:" << cc0 << " cc1:" << cc1 << " cc2:" << cc2 << " cc3:" << cc3 << endl;

							g_all_templates.push_back(t);
							templates_map[id] = make_tuple(ii, cc0, cc1, cc2, cc3);
						}
						else
						{
							//conflict!
							cout << "hit! new ii:" << ii << " cc0:" << cc0 << " cc1:" << cc1 << " cc2:" << cc2 << " cc3:" << cc3 << endl;
							tuple<int, int, int, int, int> old = templates_map[id];
							cout << "     old ii:" << get<0>(old) << " cc0:" << get<1>(old) << " cc1:" << get<2>(old) << " cc2:" << get<3>(old) << " cc3:" << get<4>(old) << endl;

							for (int x = 0; x < TemplateDim; x++)
							{
								for (int y = 0; y < TemplateDim; y++)
								{
									cout << t.v[x][y];
								}
								cout << endl;
							}
						}
					}
				}
			}
		}
	}	

	cout << "[SubTileEnumeration2] all_templates:" << g_all_templates.size() << endl;
	return true;
}

bool FloorplanSpace::SubTileEnumeration(bool addition)
{
	////enumerate all possible tile templates
	
	g_all_templates.clear();

	//5 cases on each of the four sides of a 3x3 rectangle
	//total: 5^4 = 625 types
	for (int ii = 0; ii < 625; ii++)
	{
		int side0_case = ii % 5;  //right
		int side1_case = ii / 5 % 5;  //up
		int side2_case = ii / 25 % 5;  //left
		int side3_case = ii / 125 % 5;  //bottom
		//five configs on a side: 000, 111, 110, 011, 010

		Template t(false);
		//the 3x3 core:
		for (int x = 1; x <= 3; x++)
		{
			for (int y = 1; y <= 3; y++)
			{
				t.v[x][y] = true;
			}
		}
		//right side:
		{
			if (side0_case == 1)  //111
			{
				t.v[4][1] = true;
				t.v[4][2] = true;
				t.v[4][3] = true;
			}
			else if (side0_case == 2)  //110
			{
				t.v[4][1] = true;
				t.v[4][2] = true;
				t.v[4][3] = false;
			}
			else if (side0_case == 3)  //011
			{
				t.v[4][1] = false;
				t.v[4][2] = true;
				t.v[4][3] = true;
			}
			else if (side0_case == 4)  //010
			{
				t.v[4][1] = false;
				t.v[4][2] = true;
				t.v[4][3] = false;
			}
		}
		//up side:
		{
			if (side1_case == 1)  //111
			{
				t.v[3][4] = true;
				t.v[2][4] = true;
				t.v[1][4] = true;
			}
			else if (side1_case == 2)  //110
			{
				t.v[3][4] = true;
				t.v[2][4] = true;
				t.v[1][4] = false;
			}
			else if (side1_case == 3)  //011
			{
				t.v[3][4] = false;
				t.v[2][4] = true;
				t.v[1][4] = true;
			}
			else if (side1_case == 4)  //010
			{
				t.v[3][4] = false;
				t.v[2][4] = true;
				t.v[1][4] = false;
			}
		}
		//left side:
		{
			if (side2_case == 1)  //111
			{
				t.v[0][3] = true;
				t.v[0][2] = true;
				t.v[0][1] = true;
			}
			else if (side2_case == 2)  //110
			{
				t.v[0][3] = true;
				t.v[0][2] = true;
				t.v[0][1] = false;
			}
			else if (side2_case == 3)  //011
			{
				t.v[0][3] = false;
				t.v[0][2] = true;
				t.v[0][1] = true;
			}
			else if (side2_case == 4)  //010
			{
				t.v[0][3] = false;
				t.v[0][2] = true;
				t.v[0][1] = false;
			}
		}
		//bottom side:
		{
			if (side3_case == 1)  //111
			{
				t.v[1][0] = true;
				t.v[2][0] = true;
				t.v[3][0] = true;
			}
			else if (side3_case == 2)  //110
			{
				t.v[1][0] = true;
				t.v[2][0] = true;
				t.v[3][0] = false;
			}
			else if (side3_case == 3)  //011
			{
				t.v[1][0] = false;
				t.v[2][0] = true;
				t.v[3][0] = true;
			}
			else if (side3_case == 4)  //010
			{
				t.v[1][0] = false;
				t.v[2][0] = true;
				t.v[3][0] = false;
			}
		}

		g_all_templates.push_back(t);
	}

	//add a few more templates?
	if (addition)
	{
		//four kinds of 4x4 (with the 3x3 in the four possible corners)
		for (int Case = 0; Case < 4; Case++)
		{
			Template t(false);
			//the 3x3 core:
			for (int x = 1; x <= 3; x++)
			{
				for (int y = 1; y <= 3; y++)
				{
					t.v[x][y] = true;
				}
			}

			if (Case == 0)
			{
				t.v[0][0] = true;
				t.v[1][0] = true;
				t.v[2][0] = true;
				t.v[3][0] = true;
				t.v[0][1] = true;
				t.v[0][2] = true;
				t.v[0][3] = true;
			}
			else if (Case == 1)
			{
				t.v[1][0] = true;
				t.v[2][0] = true;
				t.v[3][0] = true;
				t.v[4][0] = true;
				t.v[4][1] = true;
				t.v[4][2] = true;
				t.v[4][3] = true;
			}
			else if (Case == 2)
			{
				t.v[4][1] = true;
				t.v[4][2] = true;
				t.v[4][3] = true;
				t.v[4][4] = true;
				t.v[4][3] = true;
				t.v[4][2] = true;
				t.v[4][1] = true;
			}
			else if (Case == 3)
			{
				t.v[3][4] = true;
				t.v[2][4] = true;
				t.v[1][4] = true;
				t.v[0][4] = true;
				t.v[0][3] = true;
				t.v[0][2] = true;
				t.v[0][1] = true;
			}

			g_all_templates.push_back(t);
		}		

		//5x5 (w/ 3x3 in the middle), w/ possibly the four corners chiseled by 1x1
		for (int Case = 0; Case < 5; Case++)
		{
			Template t(false);
			for (int x = 0; x < 5; x++)
			{
				for (int y = 0; y < 5; y++)
				{
					t.v[x][y] = true;
				}
			}

			//Case == 0: no chisel
			if (Case == 1)
			{
				t.v[0][0] = false;
			}
			else if (Case == 2)
			{
				t.v[4][0] = false;
			}
			else if (Case == 3)
			{
				t.v[4][4] = false;
			}
			else if (Case == 4)
			{
				t.v[0][4] = false;
			}

			g_all_templates.push_back(t);
		}
	}

	//align all templates?
	if (false)
	{
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			g_all_templates[i].Align();
		}
	}

	//kick out redundant ones
	if (false)
	{
		vector<Template> templates_new;
		for (int i = 0; i < g_all_templates.size(); i++)
		{
			//compare to every existing one
			bool hit= false;
			for (int j = 0; j < templates_new.size(); j++)
			{
				bool same = true;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (g_all_templates[i].v[x][j] != templates_new[j].v[x][y])
						{
							same = false;
							break;
						}
					}
				}
				if (same)
				{
					hit = true;
					break;
				}
			}
			if (!hit)
			{
				templates_new.push_back(g_all_templates[i]);
			}
		}
		g_all_templates = templates_new;
	}	

	cout << "all templates:" << g_all_templates.size() << endl;

	////sub-tiles: root tile (3x3) plus 4 leaf tiles on every side

	g_leafs.clear();
	g_root.Reset();

	g_root.v[1][1] = true;
	g_root.v[2][1] = true;
	g_root.v[3][1] = true;
	g_root.v[1][2] = true;
	g_root.v[2][2] = true;
	g_root.v[3][2] = true;
	g_root.v[1][3] = true;
	g_root.v[2][3] = true;
	g_root.v[3][3] = true;

	//right side:
	{
		Template t;
		t.v[4][1] = true;
		t.v[4][2] = true;
		t.v[4][3] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[4][1] = false;
		t.v[4][2] = true;
		t.v[4][3] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[4][1] = true;
		t.v[4][2] = true;
		t.v[4][3] = false;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[4][1] = false;
		t.v[4][2] = true;
		t.v[4][3] = false;
		g_leafs.push_back(t);
	}
	//up side:
	{
		Template t;
		t.v[3][4] = true;
		t.v[2][4] = true;
		t.v[1][4] = true; 
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[3][4] = false;
		t.v[2][4] = true;
		t.v[1][4] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[3][4] = true;
		t.v[2][4] = true;
		t.v[1][4] = false;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[3][4] = false;
		t.v[2][4] = true;
		t.v[1][4] = false;
		g_leafs.push_back(t);
	}
	//left side:
	{
		Template t;
		t.v[0][3] = true;
		t.v[0][2] = true;
		t.v[0][1] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[0][3] = false;
		t.v[0][2] = true;
		t.v[0][1] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[0][3] = true;
		t.v[0][2] = true;
		t.v[0][1] = false;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[0][3] = false;
		t.v[0][2] = true;
		t.v[0][1] = false;
		g_leafs.push_back(t);
	}
	//bottom side:
	{
		Template t;
		t.v[1][0] = true;
		t.v[2][0] = true;
		t.v[3][0] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[1][0] = false;
		t.v[2][0] = true;
		t.v[3][0] = true;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[1][0] = true;
		t.v[2][0] = true;
		t.v[3][0] = false;
		g_leafs.push_back(t);
	}
	{
		Template t;
		t.v[1][0] = false;
		t.v[2][0] = true;
		t.v[3][0] = false;
		g_leafs.push_back(t);
	}

	return true;
}

bool FloorplanSpace::Floorplan::Init(int width, int height, vector<Vec2i>& off_cells)
{
	Reset();
	
	//init the active cells of the domain
	m_domain_W = width;
	m_domain_H = height;
	m_domain = new bool[m_domain_W * m_domain_H];
	memset(m_domain, true, m_domain_W * m_domain_H);  //initially all cells are active
	for (int i = 0; i < off_cells.size(); i++)  //disable certain cells
	{
		int x = off_cells[i].x;
		int y = off_cells[i].y;
		m_domain[y * m_domain_W + x] = false;
	}	

	return true;
}

bool FloorplanSpace::Floorplan::SolvePlacements()
{
	DWORD time_begin = timeGetTime();

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	////enumerate all possible tile placements

	struct Placement
	{
		int ti;  //template index
		vector<Vec2i> cells;
	};
	vector<Placement> all_placements;
	map < pair<int, int>, vector<int>> placements_map;  //key: x-y of a cell, value: placements at the cell

#pragma omp parallel for  //OpenMP
	for (int ii = 0; ii < m_domain_W * m_domain_H; ii++)
	{
		int x = ii % m_domain_W;
		int y = ii / m_domain_W;  //offset (x,y)

		for (int t = 0; t < g_all_templates.size(); t++)
		{
			Template &T = g_all_templates[t];
			Vec4i bb = T.BB();
			int right = x + (bb[1] - bb[0]);
			int top = y + (bb[3] - bb[2]);
			if (right < m_domain_W && top < m_domain_H)
			{
				//check if all occupying domain cells are active
				bool miss = false;
				for (int xx = 0; xx < TemplateDim; xx++)
				{
					for (int yy = 0; yy < TemplateDim; yy++)
					{
						if (T.v[xx][yy] && !m_domain[(y + yy)* m_domain_W + x + xx])
						{
							miss = true;
							break;
						}
					}
				}
				if (!miss)
				{
#pragma omp critical
					{
						Placement p;
						p.ti = t;

						//enumerate occupying cells and save to map
						int index = all_placements.size();
						vector<Vec2i> cells;
						for (int xx = 0; xx < TemplateDim; xx++)
						{
							for (int yy = 0; yy < TemplateDim; yy++)
							{
								if (T.v[xx][yy])
								{
									p.cells.push_back(Vec2i(x + xx, y + yy));
									placements_map[make_pair(x + xx, y + yy)].push_back(index);
								}
							}
						}
						
						all_placements.push_back(p);
					}
				}
			}
		}
	}
	cout << "[Solve] all_placements:" << all_placements.size() << endl;

	//one Boolean for every placement

	vector<GRBVar> placement_vars;
	for (int i = 0; i < all_placements.size(); i++)
	{
		placement_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
	}

	model.update();

	//objective function: maximize # of room placements

	GRBLinExpr obj;
	for (int i = 0; i < all_placements.size(); i++)
	{
		obj += placement_vars[i];
	}
	model.setObjective(obj, GRB_MAXIMIZE);
	//model.setObjective(obj, GRB_MINIMIZE);

	//constraint: every cell in the domain is covered exactly once by a tile placement (if active) or none (if off)
	for (int x = 0; x < m_domain_W; x++)
	{
		for (int y = 0; y < m_domain_H; y++)
		{
			if (placements_map.count(make_pair(x, y)) == 0)
				continue;

			vector<int> &covers = placements_map[make_pair(x, y)];
			if (covers.size() > 0)
			{
				GRBLinExpr eq;
				for (int i = 0; i < covers.size(); i++)
				{
					eq += placement_vars[covers[i]];
				}

				if (m_domain[y*m_domain_W + x])  //active domain cell?
				{
					model.addConstr(eq == 1);
				}
				else  //off domain cell?
				{
					model.addConstr(eq == 0);
				}				
			}
		}
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, 100);  //some time limit	
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[Solve] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[Solve] optimize() failed! status:%d", status);
		return false;
	}

	//get results
	for (int i = 0; i < placement_vars.size(); i++)
	{
		bool present = round(placement_vars[i].get(GRB_DoubleAttr_X));
		if (present)
		{
			//create a "tile" record
			Tile T;
			T.cells = all_placements[i].cells;			
			m_tiles.push_back(T);
		}
	}

	return true;
}

bool FloorplanSpace::Floorplan::Solve()
{
	DWORD time_begin = timeGetTime();

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//// enumerate all possible sub-tile placements

	class Tree
	{
	public:
		Vec2i offset;  //offset of the root
		vector<int> li;  //leaf lindices
	};

	vector<Tree> all_trees;
	int num_leafs = 0;

	Template root = g_root;  //the "actual" root (TODO: fix this)
	root.Align();  //note: align to (0,0) !!
	Vec2i root_offset(root.BB()[0] - g_root.BB()[0], root.BB()[2] - g_root.BB()[2]);

#pragma omp parallel for  //OpenMP
	for (int ii = 0; ii < m_domain_W*m_domain_H; ii++)
	{
		int x = ii % m_domain_W;
		int y = ii / m_domain_W;  //offset of root

		Vec4i bb = root.BB();
		int right = x + (bb[1] - bb[0]);
		int top = y + (bb[3] - bb[2]);
		if (right < m_domain_W && top < m_domain_H)
		{
			//check if all occupying domain cells are active
			bool miss = false;
			for (int xx = 0; xx < TemplateDim; xx++)
			{
				for (int yy = 0; yy < TemplateDim; yy++)
				{
					if (root.v[xx][yy] && !m_domain[(y + yy)* m_domain_W + x + xx])
					{
						miss = true;
						break;
					}
				}
			}
			if (!miss)
			{
#pragma omp critical
				{
					Tree tree;
					tree.offset = Vec2i(x, y);

					//enumerate feasible leafs of this root
					for (int l = 0; l < g_leafs.size(); l++)
					{
						Vec4i bb = g_leafs[l].BB();
						bb[0] += x + root_offset.x;
						bb[1] += x + root_offset.x;
						bb[2] += y + root_offset.y;
						bb[3] += y + root_offset.y;
						//NOTE: assume leaf is a box (convex)
						if (bb[0] >= 0 && bb[1] < m_domain_W && bb[2] >= 0 && bb[3] < m_domain_H)
						{
							//check if all occpying cells are activ
							bool miss = false;
							for (int xx = bb[0]; xx <= bb[1]; xx++)
							{
								for (int yy = bb[2]; yy <= bb[3]; yy++)
								{
									if (!m_domain[yy* m_domain_W + xx])
									{
										miss = true;
										break;
									}
								}
							}
							if (!miss)
							{
								//gotcha!
								tree.li.push_back(l);
								num_leafs++;
							}
						}
					}

					all_trees.push_back(tree);
				}
			}
		}
	}
	cout << "[Solve] trees:" << all_trees.size() << " leafs:" << num_leafs << endl;

	vector<GRBVar> tile_vars;
	vector<int> tile_root_indices;  //-1 for root tiles, >=0 if the tile is a leaf and record the root's index
	vector<vector<Vec2i>> tile_cells;

	map<pair<int, int>, vector<int>> overlaps;  //for each domain cell record the overlapping tiles

	for (int i = 0; i < all_trees.size(); i++)
	{
		Tree &tree = all_trees[i];

		//root tile:
		int root_index = tile_vars.size();
		{
			//record to overlaps
			vector<Vec2i> cell;
			for (int xx = 0; xx < TemplateDim; xx++)
			{
				for (int yy = 0; yy < TemplateDim; yy++)
				{
					if (g_root.v[xx][yy])
					{
						int x = tree.offset.x + xx + root_offset.x;
						int y = tree.offset.y + yy + root_offset.y;

						overlaps[make_pair(x, y)].push_back(root_index);
						cell.push_back(Vec2i(x, y));
					}
				}
			}
			tile_root_indices.push_back(-1);
			tile_cells.push_back(cell);
			tile_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		//leaf tiles:
		for (int j = 0; j < tree.li.size(); j++)
		{
			int leaf_index = tile_vars.size();
			//record to overlaps
			vector<Vec2i> cell;
			for (int xx = 0; xx < TemplateDim; xx++)
			{
				for (int yy = 0; yy < TemplateDim; yy++)
				{
					if (g_leafs[j].v[xx][yy])
					{
						int x = tree.offset.x + xx + root_offset.x;
						int y = tree.offset.y + yy + root_offset.y;
						overlaps[make_pair(x, y)].push_back(leaf_index);
						cell.push_back(Vec2i(x, y));
					}
				}
			}
			tile_root_indices.push_back(root_index);
			tile_cells.push_back(cell);
			tile_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));

		}
	}

	model.update();

	//obj: maximize # of roots

	GRBLinExpr obj;
	for (int i = 0; i < tile_vars.size(); i++)
	{
		if (tile_root_indices[i]<0)
			obj += tile_vars[i];
	}
	model.setObjective(obj, GRB_MAXIMIZE);
	//model.setObjective(obj, GRB_MINIMIZE);

	////constraints:

	//coverage constraint: every domain cell is coveraged exactly once (if active) or none (if not active)

	for (int x = 0; x < m_domain_W; x++)
	{
		for (int y = 0; y < m_domain_H; y++)
		{
			if (overlaps.count(make_pair(x, y)) == 0)
				continue;
			vector<int> &covers = overlaps[make_pair(x, y)];			
			if (covers.size() > 0)
			{
				GRBLinExpr eq;
				for (int i = 0; i < covers.size(); i++)
				{
					eq += tile_vars[covers[i]];
				}

				if (m_domain[y*m_domain_W + x])  //active domain cell?
				{
					model.addConstr(eq == 1);
				}
				else  //off domain cell?
				{
					model.addConstr(eq == 0);
				}
			}
		}
	}

	//a leaf cannot appear if its root does not appear
	for (int i = 0; i < tile_vars.size(); i++)
	{
		if (tile_root_indices[i] >= 0)
		{
			model.addConstr(tile_vars[i] <= tile_vars[tile_root_indices[i]]);
		}
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, g_gurobi_time);  //some time limit	
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[Solve] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[Solve] optimize() failed! status:%d", status);
		return false;
	}

	//get results

	//retrive "tiles" - each is a root plus its active leafs

	//retrive active root tiles first to create Tile place holders
	map<int, int> tile_root_map;  //key: root index, value: index in m_tiles
	for (int i = 0; i < tile_vars.size(); i++)
	{
		bool present = round(tile_vars[i].get(GRB_DoubleAttr_X));
		if (present && tile_root_indices[i] < 0/*is root*/)
		{		
			//create a "tile" record
			Tile T;
			//the root's cells:
			for (int j = 0; j < tile_cells[i].size(); j++)
			{
				T.cells.push_back(tile_cells[i][j]);
			}
			tile_root_map[i] = m_tiles.size();
			m_tiles.push_back(T);

		}
	}

	//retrive active leaf tiles to fill tiles
	for (int i = 0; i < tile_vars.size(); i++)
	{
		bool present = round(tile_vars[i].get(GRB_DoubleAttr_X));
		if (present && tile_root_indices[i] >= 0)
		{
			//add to the corresponding Tile's cells:
			Tile &T = m_tiles[tile_root_map[tile_root_indices[i]]];
			for (int j = 0; j < tile_cells[i].size(); j++)
			{
				T.cells.push_back(tile_cells[i][j]);
			}
		}
	}

	AddText("[Solve] done. m_tiles:%d time:%d", m_tiles.size(), timeGetTime() - time_begin);
	return true;
}

bool FloorplanSpace::Floorplan::Solve2()
{
	DWORD time_begin = timeGetTime();

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//enumerate all possible forest placements and the valid roots and leafs of each forest

	class ForestPlacement
	{
	public:
		Vec2i offset;  //offset of the forest
		vector<int> rs;  //valid root indices (among all roots of g_forest)
		vector<int> ls;  //valid leaf indices (among all leafs of g_forest)
		vector<vector<int>> rs_of_leafs;  //rs of each leaf's linked roots
	};
	vector<ForestPlacement> all_forest_placements;

	int num_roots = 0, num_leafs = 0;  //for debug
	for (int ii = 0; ii < m_domain_W*m_domain_H; ii++)
	{
		int X = ii % m_domain_W;
		int Y = ii / m_domain_W;  //offset of root

		ForestPlacement fp;
		fp.offset = Vec2i(X, Y);

		//enumerate compatible roots
		map<int, int> rs_map;  //key: root index, value: index in rs
		for (int j = 0; j < g_forest.roots.size(); j++)
		{
			Template &root = g_forest.roots[j];

			//check BB first
			Vec4i bb = root.BB();
			if (bb[0] + X < 0 || bb[1] + X >= m_domain_W ||
				bb[2] + Y < 0 || bb[3] + Y >= m_domain_H)			
				continue;
			
			//any cell of the leaf that is not in the domain?
			bool miss = false;
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					int xx = X + x;
					int yy = Y + y;
					if (root.v[x][y] && !m_domain[yy * m_domain_W + xx])
					{
						miss = true;
						break;
					}
				}
			}
			if (!miss)
			{
				rs_map[j] = fp.rs.size();
				fp.rs.push_back(j);
				num_roots++;
			}
		}

		//skip forests w/o valid roots
		if (fp.rs.size() == 0)
			continue;

		//enumerate compatible leafs
		for (int j = 0; j < g_forest.leafs.size(); j++)
		{
			//prepare root_of_leafs:
			vector<int> ss;
			for (int k = 0; k < g_forest.links[j].size(); k++)
			{
				if (rs_map.count(g_forest.links[j][k]) > 0)
					ss.push_back(rs_map[g_forest.links[j][k]]);
			}
			if (ss.size() == 0)
				continue;  //no valid roots at all? skip

			Template &leafs = g_forest.leafs[j];

			//check BB first
			Vec4i bb = leafs.BB();
			if (bb[0] + X < 0 || bb[1] + X >= m_domain_W ||
				bb[2] + Y < 0 || bb[3] + Y >= m_domain_H)
				continue;
			
			//any cell of the leaf that is not in the domain?
			bool miss = false;
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					int xx = X + x;
					int yy = Y + y;
					if (leafs.v[x][y] && !m_domain[yy * m_domain_W + xx])
					{
						miss = true;
						break;
					}
				}
			}
			if (!miss)
			{
				fp.ls.push_back(j);
				fp.rs_of_leafs.push_back(ss);
				num_leafs++;				
			}
		}

		all_forest_placements.push_back(fp);
	}
	cout << "[Solve] all fps:" << all_forest_placements.size() << " roots:" << num_roots << " leafs:" << num_leafs << endl;

	//create the vars and associated indices

	vector<GRBVar> tile_vars;  //all tile vars including both root and leaf tiles
	vector<vector<Vec2i>> tile_cells;  //occupying (x,y) of each tile
	vector<int> forest_indices;  //#forest for every tile var
	vector<vector<int>> root_indices;  //for ever tile var, if non-empty, it is a leaf tile, store the indices of its connected root vars	
	vector<vector<int>> non_root_indices;  //for ever tile var, if non-empty, it is a leaf tile, store the indices of its non-connected root vars	
	map<pair<int, int>, vector<int>> overlaps;  //for each domain cell record the overlapping tiles
	vector<vector<int>> root_groups;  //indices of root tiles of each forest

	for (int i = 0; i < all_forest_placements.size(); i++)
	{
		ForestPlacement &fp = all_forest_placements[i];

		vector<int> ris_all;  //root tile indices
		for (int j = 0; j < fp.rs.size(); j++)
		{
			ris_all.push_back(tile_vars.size());

			forest_indices.push_back(i);
			
			vector<int> empty;
			root_indices.push_back(empty);  //empty list
			non_root_indices.push_back(empty);

			vector<Vec2i> cells;
			Template &root = g_forest.roots[fp.rs[j]];
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					if (root.v[x][y])
					{
						int xx = fp.offset.x + x;
						int yy = fp.offset.y + y;
						cells.push_back(Vec2i(xx, yy));
						overlaps[make_pair(xx, yy)].push_back(tile_vars.size());
					}
				}
			}			

			tile_cells.push_back(cells);
			
			tile_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}
		
		root_groups.push_back(ris_all);

		for (int j = 0; j < fp.ls.size(); j++)
		{
			forest_indices.push_back(i);

			//collect corresponding root indices
			vector<int> ris;
			for (int k = 0; k < fp.rs_of_leafs[j].size(); k++)
			{
				ris.push_back(ris_all[fp.rs_of_leafs[j][k]]);
			}			
			root_indices.push_back(ris);	
			//and non-connected root indices
			vector<int> non_ris;
			for (int k = 0; k < ris_all.size(); k++)
			{
				bool found = false;
				for (int q = 0; q < fp.rs_of_leafs[j].size(); q++)
				{
					if (fp.rs_of_leafs[j][q] == k)
					{
						found = true;
						break;
					}
				}
				if (!found)
					non_ris.push_back(ris_all[k]);				
			}
			non_root_indices.push_back(non_ris);

			vector<Vec2i> cells;
			Template &leaf = g_forest.leafs[fp.ls[j]];
			for (int x = 0; x < TemplateDim; x++)
			{
				for (int y = 0; y < TemplateDim; y++)
				{
					if (leaf.v[x][y])
					{
						int xx = fp.offset.x + x;
						int yy = fp.offset.y + y;
						cells.push_back(Vec2i(xx, yy));
						overlaps[make_pair(xx, yy)].push_back(tile_vars.size());
					}
				}
			}

			tile_cells.push_back(cells);			

			tile_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}
	}

	model.update();

	//obj: maximize # of roots (because each room has exactly one root)
	GRBLinExpr obj;
	for (int i = 0; i < tile_vars.size(); i++)
	{
		if (root_indices[i].size() == 0)
			obj += tile_vars[i];
	}
	model.setObjective(obj, GRB_MAXIMIZE);
	//model.setObjective(obj, GRB_MINIMIZE);

	////constraints:

	//coverage constraint: every domain cell is coveraged exactly once (if active) or none (if not active)

	for (int x = 0; x < m_domain_W; x++)
	{
		for (int y = 0; y < m_domain_H; y++)
		{
			if (overlaps.count(make_pair(x, y)) == 0)
				continue;

			vector<int> &covers = overlaps[make_pair(x, y)];
			if (covers.size() > 0)
			{
				GRBLinExpr eq;
				for (int i = 0; i < covers.size(); i++)
				{
					eq += tile_vars[covers[i]];
				}

				if (m_domain[y*m_domain_W + x])  //active domain cell?
				{
					model.addConstr(eq == 1);
				}
				else  //off domain cell?
				{
					model.addConstr(eq == 0);
				}
			}
		}
	}
	
	//leaf topology constraints
	for (int i = 0; i < tile_vars.size(); i++)
	{
		if (root_indices[i].size() == 0)  //must be leaf
			continue;
		
		//a leaf can appear only if any of its root(s) appear
		{
			GRBLinExpr sum;
			for (int j = 0; j < root_indices[i].size(); j++)
			{
				sum += tile_vars[root_indices[i][j]];
			}
			model.addConstr(tile_vars[i] <= sum);
		}
		//and cannot appear if any of its non-root(s) appear
		/*{			
			GRBLinExpr sum;
			for (int j = 0; j < non_root_indices[i].size(); j++)
			{
				sum += tile_vars[non_root_indices[i][j]];
			}
			model.addConstr(tile_vars[i] <= 1 - sum);
		}*/
	}	

	//for every root group, at most one of them can be true
	for (int i = 0; i < root_groups.size(); i++)
	{
		GRBLinExpr eq;
		vector<int> &group = root_groups[i];
		for (int i = 0; i < group.size(); i++)
		{
			eq += tile_vars[group[i]];
		}
		model.addConstr(eq <= 1);
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, 200);  //some time limit	
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[Solve] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[Solve] optimize() failed! status:%d", status);
		return false;
	}

	//get results

	//retrive "tiles" - each is a root plus its active leafs (of the same forest)

	//retrive active root tiles first to create Tile place holders
	map<int, int> tile_forest_map;  //key: forest index, value: index in m_tiles

	for (int i = 0; i < tile_vars.size(); i++)
	{
		bool present = round(tile_vars[i].get(GRB_DoubleAttr_X));
		if (present)
		{
			cout << "tile#" << i << " root?" << root_indices[i].size() << " cells:" << tile_cells[i].size() << endl;
			int forest_index = forest_indices[i];
			if (tile_forest_map.count(forest_index) == 0)
			{
				//create a new "tile" record
				Tile T;
				for (int j = 0; j < tile_cells[i].size(); j++)
				{
					T.cells.push_back(tile_cells[i][j]);
				}
				tile_forest_map[forest_index] = m_tiles.size();								

				m_tiles.push_back(T);				
			}
			else
			{
				//just add to the existing tile record
				Tile &T = m_tiles[tile_forest_map[forest_index]];
				for (int j = 0; j < tile_cells[i].size(); j++)
				{
					T.cells.push_back(tile_cells[i][j]);
				}
			}
		}
	}

	AddText("[Solve] done. m_tiles:%d time:%d", m_tiles.size(), timeGetTime() - time_begin);
	return true;
}

bool FloorplanSpace::Floorplan::SolveNegativeLeafs()
{
	//6 roots:
	vector<Template> roots;
	for (int Case = 0; Case < 6; Case++)
	{
		Template root;
		if (Case == 0)
		{
			for (int x = 0; x <= 2; x++)
			{
				for (int y = 0; y <= 2; y++)
				{
					root.v[x][y] = true;
				}
			}
		}
		else if (Case == 1)
		{
			for (int x = 0; x <= 3; x++)
			{
				for (int y = 0; y <= 2; y++)
				{
					root.v[x][y] = true;
				}
			}
		}
		else if (Case == 2)
		{
			for (int x = 0; x <= 2; x++)
			{
				for (int y = 0; y <= 3; y++)
				{
					root.v[x][y] = true;
				}
			}
		}
		else if (Case == 3)
		{
			for (int x = 0; x <= 3; x++)
			{
				for (int y = 0; y <= 3; y++)
				{
					root.v[x][y] = true;
				}
			}
		}
		else if (Case == 4)
		{
			for (int x = 0; x <= 4; x++)
			{
				for (int y = 0; y <= 3; y++)
				{
					root.v[x][y] = true;
				}
			}
		}
		else if (Case == 5)
		{
			for (int x = 0; x <= 3; x++)
			{
				for (int y = 0; y <= 4; y++)
				{
					root.v[x][y] = true;
				}
			}
		}
		roots.push_back(root);
	}

	//enumerate leafs:
	vector<Template> leafs;	
	vector < vector<int> > roots_of_leafs;  //roots of every leaf
	for (int Place = 0; Place < 13; Place++)
	{		
		Vec2i corner;

		if (Place == 0)
			corner = Vec2i(0, 0);
		else if (Place == 1)
			corner = Vec2i(2, 0);
		else if (Place == 2)
			corner = Vec2i(3, 0);
		else if (Place == 3)
			corner = Vec2i(4, 0);
		else if (Place == 4)
			corner = Vec2i(0, 2);
		else if (Place == 5)
			corner = Vec2i(2, 2);
		else if (Place == 6)
			corner = Vec2i(3, 2);		
		else if (Place == 7)
			corner = Vec2i(0, 3);
		else if (Place == 8)
			corner = Vec2i(2, 3);
		else if (Place == 9)
			corner = Vec2i(3, 3);
		else if (Place == 10)
			corner = Vec2i(4, 3);
		else if (Place == 11)
			corner = Vec2i(0, 4);		
		else if (Place == 12)
			corner = Vec2i(3, 4);

		//16 shape indices (some of them are the same shapes):		
		vector<int> shapes;  //each is 1x1, 2x1, 1x2, and 2x2

		if (Place == 0)
		{
			shapes.push_back(0);
			shapes.push_back(1);
			shapes.push_back(2);
			shapes.push_back(3);
		}
		else if (Place == 1)
		{
			shapes.push_back(4);
			shapes.push_back(6);
		}
		else if (Place == 2 || Place == 3)
		{
			shapes.push_back(4);
			shapes.push_back(5);
			shapes.push_back(6);
			shapes.push_back(7);
		}
		else if (Place == 4)
		{
			shapes.push_back(8);
			shapes.push_back(9);
		}
		else if (Place == 5)
		{
			shapes.push_back(12);			
		}
		else if (Place == 6)
		{
			shapes.push_back(12);
			shapes.push_back(13);
		}		
		else if (Place == 7)
		{
			shapes.push_back(8);
			shapes.push_back(9);
			shapes.push_back(10);
			shapes.push_back(11);
		}
		else if (Place == 8)
		{
			shapes.push_back(12);
			shapes.push_back(14);
		}
		else if (Place == 9 || Place == 10)
		{
			shapes.push_back(12);
			shapes.push_back(13);
			shapes.push_back(14);
			shapes.push_back(15);
		}
		else if (Place == 11)
		{
			shapes.push_back(8);
			shapes.push_back(9);
			shapes.push_back(10);
			shapes.push_back(11);
		}
		else if (Place == 12)
		{
			shapes.push_back(12);
			shapes.push_back(13);
			shapes.push_back(14);
			shapes.push_back(15);
		}

		for (int ii = 0; ii < shapes.size(); ii++)
		{
			Template leaf;

			int s = shapes[ii];
			if (s == 0 || s == 4 || s == 8 || s == 12)
			{
				//a single 1x1
				leaf.v[corner.x][corner.y] = true;				
			}
			else if (s == 1 || s == 9)
			{
				//2x1 to right
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x + 1][corner.y] = true;				
			}
			else if (s == 2 || s == 6)
			{
				//1x2 to up
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x][corner.y + 1] = true;				
			}
			else if (s == 3)
			{
				//2x2 to right up
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x + 1][corner.y] = true;
				leaf.v[corner.x + 1][corner.y + 1] = true;
				leaf.v[corner.x][corner.y + 1] = true;
			}
			else if (s == 5 || s == 13)
			{
				//2x1 to left
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x - 1][corner.y] = true;
			}			
			else if (s == 7)
			{
				//2x2 to left up
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x - 1][corner.y] = true;
				leaf.v[corner.x - 1][corner.y + 1] = true;
				leaf.v[corner.x][corner.y + 1] = true;
			}
			else if (s == 10 || s == 14)
			{
				//1x2 to down
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x][corner.y - 1] = true;
			}
			else if (s == 11)
			{
				//2x2 to right down
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x + 1][corner.y] = true;
				leaf.v[corner.x + 1][corner.y - 1] = true;
				leaf.v[corner.x][corner.y - 1] = true;
			}
			else if (s == 15)
			{
				//2x2 to left down
				leaf.v[corner.x][corner.y] = true;
				leaf.v[corner.x - 1][corner.y] = true;
				leaf.v[corner.x - 1][corner.y - 1] = true;
				leaf.v[corner.x][corner.y - 1] = true;
			}			

			leafs.push_back(leaf);			

			vector<int> rs;  //applicable roots
			if (Place == 0)
			{
				if (s == 0)
				{
					rs.push_back(0);
					rs.push_back(1);
					rs.push_back(2);
					rs.push_back(3);
					rs.push_back(4);
					rs.push_back(5);
				}
				else if (s == 1)
				{
					//rs.push_back(0);
					rs.push_back(1);
					//rs.push_back(2);
					rs.push_back(3);
					rs.push_back(4);
					rs.push_back(5);
				}
				else if (s == 2)
				{
					//rs.push_back(0);
					//rs.push_back(1);
					rs.push_back(2);
					rs.push_back(3);
					rs.push_back(4);
					rs.push_back(5);
				}
				else if (s == 3)
				{
					//rs.push_back(0);
					//rs.push_back(1);
					//rs.push_back(2);
					rs.push_back(3);
					rs.push_back(4);
					rs.push_back(5);
				}
			}
			else if (Place == 1)
			{
				if (s == 4)
				{
					rs.push_back(0);					
					rs.push_back(2);					
				}
				else if (s == 6)
				{
					//rs.push_back(0);					
					rs.push_back(2);
				}
			}
			else if (Place == 2)
			{
				if (s == 4)
				{
					rs.push_back(1);
					rs.push_back(3);
				}
				else if (s == 5)
				{
					rs.push_back(1);
					rs.push_back(3);
				}
				else if (s == 6 || s == 7)
				{
					//rs.push_back(1);
					rs.push_back(3);
				}				
			}
			else if (Place == 3)
			{
				rs.push_back(4);
				rs.push_back(5);								
			}
			else if (Place == 4)
			{
				if (s == 8)
				{
					rs.push_back(0);
					rs.push_back(1);
				}
				else if (s == 9)
				{
					//rs.push_back(0);
					rs.push_back(1);
				}				
			}
			else if (Place == 5)
			{
				rs.push_back(0);				
			}
			else if (Place == 6)
			{
				rs.push_back(1);				
			}
			else if (Place == 7)
			{
				if (s == 8 || s == 10)
				{
					rs.push_back(2);
					rs.push_back(3);
					rs.push_back(4);
				}
				else if (s == 9 || s == 11)
				{
					//rs.push_back(2);
					rs.push_back(3);
					rs.push_back(4);
				}				
			}
			else if (Place == 8)
			{
				rs.push_back(2);				
			}
			else if (Place == 9)
			{
				rs.push_back(3);
			}
			else if (Place == 10)
			{
				rs.push_back(4);
			}
			else if (Place == 11)
			{
				rs.push_back(5);				
			}
			else if (Place == 12)
			{
				rs.push_back(5);
			}

			roots_of_leafs.push_back(rs);
		}
	}
	cout << "[SolveNegativeLeafs] leafs:" << leafs.size() << endl;

	//enumerate leafs of each roots
	vector< vector<int> > leafs_of_roots(roots.size());
	for (int i = 0; i < leafs.size(); i++)
	{
		vector<int> &rs = roots_of_leafs[i];
		for (int j = 0; j < rs.size(); j++)
		{
			leafs_of_roots[rs[j]].push_back(i);
		}
	}
	
	//collect leaf-leaf conflict pairs:
	map<pair<int, int>, bool> conflict_leafs;
	for (int ii = 0; ii < roots.size(); ii++)
	{
		vector<int> &leafs_here = leafs_of_roots[ii];		
		for (int i = 0; i < leafs_here.size(); i++)
		{
			for (int j = 0; j < leafs_here.size(); j++)
			{
				if (i == j)
					continue;

				Template &leaf0 = leafs[leafs_here[i]];
				Template &leaf1 = leafs[leafs_here[j]];

				//overlapping pairs are for sure conflict
				bool overlap = false;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (leaf0.v[x][y] && leaf1.v[x][y])
						{
							overlap = true;
							break;
						}
					}
				}
				if (overlap)
				{
					if (conflict_leafs.count(make_pair(leafs_here[i], leafs_here[j])) == 0 &&
						conflict_leafs.count(make_pair(leafs_here[j], leafs_here[i])) == 0)
					{
						conflict_leafs[make_pair(leafs_here[i], leafs_here[j])] = true;
					}

					continue;
				}

				//if the two leafs are too close, conflict			
				bool gotcha = false;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (leaf0.v[x][y])
						{
							//to the left x1:
							if (x >= 1 && leaf1.v[x - 1][y])
							{
								gotcha = true;
								break;
							}
							//left x2
							if (x >= 2 && leaf1.v[x - 2][y])
							{
								gotcha = true;
								break;
							}
							//right x1
							if (x < TemplateDim - 1 && leaf1.v[x + 1][y])
							{
								gotcha = true;
								break;
							}
							//right x2
							if (x < TemplateDim - 2 && leaf1.v[x + 2][y])
							{
								gotcha = true;
								break;
							}
							//down x1
							if (y >= 1 && leaf1.v[x][y - 1])
							{
								gotcha = true;
								break;
							}
							//down x2
							if (y >= 2 && leaf1.v[x][y - 2])
							{
								gotcha = true;
								break;
							}
							//up x1
							if (y < TemplateDim - 1 && leaf1.v[x][y + 1])
							{
								gotcha = true;
								break;
							}
							//up x2
							if (y < TemplateDim - 2 && leaf1.v[x][y + 2])
							{
								gotcha = true;
								break;
							}
						}						
					}
					if (gotcha)
						break;
				}

				if (gotcha)
				{
					if (conflict_leafs.count(make_pair(leafs_here[i], leafs_here[j])) == 0 &&
						conflict_leafs.count(make_pair(leafs_here[j], leafs_here[i])) == 0)
					{
						conflict_leafs[make_pair(leafs_here[i], leafs_here[j])] = true;
					}
				}
			}
		}
	}
	cout << "[SolveNegative] conflict_leafs:" << conflict_leafs.size() << endl;

	//solve the tiling IP w/ roots and negative leafs

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	vector<GRBVar> tile_vars;  //root tiles and negative leaf tiles
	vector<vector<int>> tile_roots;  //root tile indices of leaf tiles, empty for root tiles
	map<pair<int, int>, vector<int>> overlap_tiles;  //indices of tiles (roots and leafs) at every (X,Y)	
	vector<pair<int, int>> conflict_tiles;  //pairs of conflict leaf tiles
	vector<vector<int>> root_groups;  //root tile indices of the same forest
	int num_roots = 0, num_leafs = 0;  //for debug
	for (int X = 0; X < m_domain_W; X++)
	{
		for (int Y = 0; Y < m_domain_H; Y++)
		{
			//(X,Y) is the lower-left corner

			//create root tiles:
			map<int, int> root_map;  //root tile indices. key: # in roots, value: index in tile_vars
			for (int ii = 0; ii < roots.size(); ii++)
			{
				Template &root = roots[ii];

				//check BB first
				Vec4i bb = root.BB();
				if (bb[0] + X < 0 || bb[1] + X >= m_domain_W ||
					bb[2] + Y < 0 || bb[3] + Y >= m_domain_H)
					continue;

				//any cell of the leaf that is not in the domain?
				bool miss = false;
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						int xx = X + x;
						int yy = Y + y;
						if (root.v[x][y] && !m_domain[yy * m_domain_W + xx])
						{
							miss = true;
							break;
						}
					}
				}

				if (!miss)
				{
					int index = tile_vars.size();

					root_map[ii] = index;
					
					tile_roots.push_back(vector<int>());  //empty list for root tiles
					
					//save to overlap_tiles
					for (int x = 0; x < TemplateDim; x++)
					{
						for (int y = 0; y < TemplateDim; y++)
						{
							if (root.v[x][y])
							{
								int xx = X + x;
								int yy = Y + y;
								overlap_tiles[make_pair(xx, yy)].push_back(index);
							}
						}
					}
					
					num_roots++;
					tile_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
				}
			}

			vector<int> root_group;
			for (map<int, int>::iterator itr = root_map.begin(); itr != root_map.end(); itr++)
			{
				root_group.push_back((*itr).second);
			}
			if (root_group.empty())
				continue;  //no roots at all!
			else
				root_groups.push_back(root_group);

			//create leaf tiles:
			map<int, int> leaf_map;  //key: # in leafs, value: tile indices
			for (int ii = 0; ii < leafs.size(); ii++)
			{				
				//is any of its linked roots feasible?
				vector<int> rts;  //linked root tile indices
				vector<int> &roots_here = roots_of_leafs[ii];
				for (int jj = 0; jj < roots_here.size(); jj++)
				{
					if (root_map.count(roots_here[jj]) > 0)
					{
						rts.push_back(root_map[roots_here[jj]]);
					}
				}
				if (rts.empty())
					continue;

				int index = tile_vars.size();

				tile_roots.push_back(rts);

				//save to overlap_tiles
				Template &leaf = leafs[ii];
				for (int x = 0; x < TemplateDim; x++)
				{
					for (int y = 0; y < TemplateDim; y++)
					{
						if (leaf.v[x][y])
						{
							int xx = X + x;
							int yy = Y + y;
							overlap_tiles[make_pair(xx, yy)].push_back(index);
						}
					}
				}

				leaf_map[ii] = index;

				num_leafs++;
				tile_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
			}

			//create conflict tile pairs
			for (map<pair<int, int>, bool>::iterator itr = conflict_leafs.begin(); itr != conflict_leafs.end(); itr++)
			{
				int i0 = (*itr).first.first;
				int i1 = (*itr).first.second;

				if (leaf_map.count(i0) > 0 && leaf_map.count(i1) > 0)
				{
					conflict_tiles.push_back(make_pair(leaf_map[i0], leaf_map[i1]));
				}
			}
		}
	}
	cout << "[SolveNegative] tile_vars:" << tile_vars.size() << " roots:" << num_roots << " leafs:" << num_leafs << " conflicts:" << conflict_tiles.size() << endl;

	model.update();

	//obj function: max/min # of root tiles
	GRBLinExpr obj;
	for (int i = 0; i < tile_vars.size(); i++)
	{
		if (tile_roots[i].empty())
		{
			obj += tile_vars[i];
		}
	}
	model.setObjective(obj, GRB_MAXIMIZE);

	////constraints:

	//per cell overlap counting:
	for (int X = 0; X < m_domain_W; X++)
	{
		for (int Y = 0; Y < m_domain_H; Y++)
		{
			if (!m_domain[Y * m_domain_W + X])
				continue;  //invalid cell

			if (overlap_tiles.count(make_pair(X, Y)) == 0)
			{
				cout << "nope overlap for X" << X << "Y" << Y << endl;
				return false;
			}

			GRBLinExpr eq;

			vector<int> &tis = overlap_tiles[make_pair(X, Y)];
			for (int i = 0; i < tis.size(); i++)
			{
				int ti = tis[i];
				//root or leaf?
				bool is_root = tile_roots[ti].empty();
				if (is_root)  //add root vars
				{
					eq += tile_vars[ti];
				}
				else  //subtract leaf vars
				{
					eq -= tile_vars[ti];
				}
			}

			model.addConstr(eq == 1);
		}
	}
	
	//for every leaf var, it cannot appear unless one of its root vars is present
	for (int i = 0; i < tile_vars.size(); i++)
	{
		if (tile_roots[i].empty())
			continue;

		GRBLinExpr sum;
		vector<int> &rts = tile_roots[i];
		for (int j = 0; j < rts.size(); j++)
		{
			sum += tile_vars[rts[j]];
		}

		model.addConstr(tile_vars[i] <= sum);
	}

	//each root group appeat at most once
	for (int i = 0; i < root_groups.size(); i++)
	{
		GRBLinExpr sum;
		for (int j = 0; j < root_groups[i].size(); j++)
		{
			sum += tile_vars[root_groups[i][j]];
		}
		model.addConstr(sum <= 1);
	}

	//conflict (leaf) var pairs:
	for (int i = 0; i < conflict_tiles.size(); i++)
	{
		model.addConstr(tile_vars[conflict_tiles[i].first] + tile_vars[conflict_tiles[i].second] <= 1);
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, 800);  //some time limit
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[SolveNegative] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[SolveNegative] optimize() failed! status:%d", status);
		return false;
	}

	//get results
	int count = 0;
	for (int i = 0; i < tile_vars.size(); i++)
	{
		bool present = round(tile_vars[i].get(GRB_DoubleAttr_X));
		if (present)
		{
			bool is_root = tile_roots[i].empty();

			cout << "tile#" << i << " is present root?" << is_root << endl;

		}
	}

	return true;
}

bool FloorplanSpace::Floorplan::SolveNonConvex()
{
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	try
	{
		GRBVar x0, y0, x1, y1, x2, y2;
		x0 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);
		y0 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);
		x1 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);
		y1 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);
		//x1 = model.addVar(0, 1, 0, GRB_BINARY);
		//y1 = model.addVar(0, 1, 0, GRB_BINARY);
		x2 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);
		y2 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);

		GRBVar len0, len1;
		len0 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);
		len1 = model.addVar(-1000, 1000, 0, GRB_CONTINUOUS);

		GRBVar cos0, sin0, cos1, sin1;
		cos0 = model.addVar(-1, 1, 0, GRB_CONTINUOUS);
		sin0 = model.addVar(-1, 1, 0, GRB_CONTINUOUS);
		cos1 = model.addVar(-1, 1, 0, GRB_CONTINUOUS);
		sin1 = model.addVar(-1, 1, 0, GRB_CONTINUOUS);

		model.update();

		model.setObjective(x0 * x0 + y0 * y0);

		//length constraints:
		model.addQConstr((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) == len0 * len0);
		model.addQConstr((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) == len1 * len1);

		// cos/sin constraints:
		model.addQConstr((x1 - x0) == len0 * cos0);
		model.addQConstr((y1 - y0) == len0 * sin0);
		model.addQConstr((x2 - x1) == len1 * cos1);
		model.addQConstr((y2 - y1) == len1 * sin1);

		//user constraints:
		model.addConstr(len0 == 2);
		model.addConstr(len1 == 3);

		model.set(GRB_IntParam_NonConvex, 2);

		cout << "optimize!" << endl;
		model.optimize();
		cout << "optimize done" << endl;

		int status = model.get(GRB_IntAttr_Status);
		cout << "status:" << status << endl;

		cout << "v0:" << x0.get(GRB_DoubleAttr_X) << "," << y0.get(GRB_DoubleAttr_X) << endl;
		cout << " v1:" << x1.get(GRB_DoubleAttr_X) << "," << y1.get(GRB_DoubleAttr_X) << endl;
		cout << " v2:" << x2.get(GRB_DoubleAttr_X) << "," << y2.get(GRB_DoubleAttr_X) << endl;
		cout << "cos/sin#0:" << cos0.get(GRB_DoubleAttr_X) << "," << sin0.get(GRB_DoubleAttr_X) << endl;
		cout << "cos/sin#1:" << cos1.get(GRB_DoubleAttr_X) << "," << sin1.get(GRB_DoubleAttr_X) << endl;
	}
	catch (GRBException e)
	{
		cout << "excception:" << e.getMessage() << endl;
	}

	return false;
}

bool FloorplanSpace::Floorplan::SolveRectangles()
{
	const int W_min = 6, W_max = 15;
	const int H_min = 6, H_max = 15;	

	//a potential rectangle w/ offset (x,y) and valid horizontal (W) and vertical (H) masks
	class Rectangle
	{
	public:
		Vec2i offset;  //lower-left corner	
		vector<int> valid_Ws;
		vector<int> valid_Hs;		
	};
	vector<Rectangle> rects;

	//create all potential rectangle placements
	for (int ii = 0; ii < m_domain_W*m_domain_H; ii++)
	{
		int X = ii % m_domain_W;
		int Y = ii / m_domain_W;  //offset of root

		Rectangle rect;
		rect.offset = Vec2i(X, Y);
		
		//enumerate valid horizontal-spanning masks:
		for (int W = W_min; W <= W_max; W++)
		{
			//first check BB
			if (X + W > m_domain_W)
				continue;
			//then check if covers some invalid domain cells
			bool nope = false;
			for (int xx = X; xx < X + W; xx++)
			{
				for (int yy = Y; yy < MIN2(Y + H_max, m_domain_H); yy++)
				{
					if (!m_domain[yy * m_domain_W + xx])
					{
						nope = true;
						break;
					}
				}
			}
			if (!nope)
			{
				rect.valid_Ws.push_back(W);
			}
		}

		//enumerate valid vertical-spanning masks:
		for (int H = H_min; H <= H_max; H++)
		{
			//first check BB
			if (Y + H > m_domain_H)
				continue;
			//then check if covers some invalid domain cells
			bool nope = false;
			for (int xx = X; xx < MIN2(X + W_max, m_domain_W); xx++)
			{
				for (int yy = Y; yy < Y + H; yy++)
				{
					if (!m_domain[yy * m_domain_W + xx])
					{
						nope = true;
						break;
					}
				}
			}
			if (!nope)
			{
				rect.valid_Hs.push_back(H);
			}
		}		

		if(rect.valid_Ws.size()>0 && rect.valid_Hs.size()>0)
			rects.push_back(rect);  //include rects w/ both valid Ws and Hs
	}

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//create rects' flag variables
	vector<GRBVar> rect_vars;
	for (int i = 0; i < rects.size(); i++)
	{
		rect_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
	}

	//create the mask variables
	vector<GRBVar> mask_vars;
	vector<pair< vector<int>, vector<int> >> mask_vars_of_rects;  //for every rect, store the indices of its W mask vars and H mask vars	
	map<pair<int,int>/*(x,y)*/, vector<int>> overlap_masks;  //overlapping masks at every (x,y)
	map<pair<int, int>, map<int, bool>> overlap_rects;  //overlapping rectangles at every (x,y)

	for (int i = 0; i < rects.size(); i++)
	{
		Rectangle &rect = rects[i];		

		vector<int> W_mask_vars;
		for (int j = 0; j < rect.valid_Ws.size(); j++)  //horizontal-spanning masks
		{
			int index = mask_vars.size();

			W_mask_vars.push_back(index);

			//enumerate cells of this mask and save to overlaps
			for (int xx = rect.offset.x; xx < rect.offset.x + rect.valid_Ws[j]; xx++)
			{
				for (int yy = rect.offset.y; yy < MIN2(rect.offset.y + H_max, m_domain_H); yy++)
				{
					overlap_masks[make_pair(xx, yy)].push_back(index);
					overlap_rects[make_pair(xx, yy)][i] = true;					
				}
			}
			
			mask_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		vector<int> H_mask_vars;
		for (int j = 0; j < rect.valid_Hs.size(); j++)  //vertical-spanning masks
		{
			int index = mask_vars.size();

			H_mask_vars.push_back(index);

			//enumerate cells of this mask and save to overlaps
			for (int xx = rect.offset.x; xx < MIN2(rect.offset.x + W_max, m_domain_W); xx++)
			{
				for (int yy = rect.offset.y; yy < rect.offset.y + rect.valid_Hs[j]; yy++)
				{
					overlap_masks[make_pair(xx, yy)].push_back(index);
					overlap_rects[make_pair(xx, yy)][i] = true;					
				}
			}			

			mask_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}	

		mask_vars_of_rects.push_back(make_pair(W_mask_vars, H_mask_vars));
	}

	cout << "[SolveRectnagles] rects:" << rect_vars.size() << " masks:" << mask_vars.size() << endl;

	model.update();
	
	//obj: max or min # of rectangles
	GRBLinExpr obj;
	for (int i = 0; i < rect_vars.size(); i++)
	{		
		obj += rect_vars[i];
	}
	model.setObjective(obj, GRB_MAXIMIZE);

	//constraints:

	//for every rect's Hs and Ws, sum to the rect's flag var
	for (int i = 0; i < rects.size(); i++)
	{
		GRBVar rect_var = rect_vars[i];

		for (int Case = 0; Case < 2; Case++)  //Ws or Hs
		{
			vector<int> *vars = NULL;
			if (Case == 0)
				vars = &mask_vars_of_rects[i].first;
			else if (Case == 1)
				vars = &mask_vars_of_rects[i].second;
			
			GRBLinExpr sum;
			for (int j = 0; j < (*vars).size(); j++)
			{
				sum += mask_vars[(*vars)[j]];
			}
			model.addConstr(sum == rect_var);
		}
	}

	//coverage constraint: for every cell, it is exactly covered "once" by a rectangle
	for (int x = 0; x < m_domain_W; x++)
	{
		for (int y = 0; y < m_domain_H; y++)
		{
			if (x == 7 && y == 7)
				cout << "debug";

			if (!m_domain[y * m_domain_W + x])
				continue;

			GRBLinExpr sum;			

			//sum up mask vars at this cell
			if (overlap_masks.count(make_pair(x, y)) > 0)
			{
				vector<int> &mask_vars_here = overlap_masks[make_pair(x, y)];
				for (int k = 0; k < mask_vars_here.size(); k++)
				{
					sum += mask_vars[mask_vars_here[k]];
				}
			}

			//minus rect vars at this cell
			if (overlap_rects.count(make_pair(x, y)) > 0)
			{
				map<int, bool> &rect_vars_here = overlap_rects[make_pair(x, y)];
				for (map<int, bool>::iterator itr = rect_vars_here.begin(); itr != rect_vars_here.end(); itr++)
				{
					sum -= rect_vars[(*itr).first];
				}
			}
			
			model.addConstr(sum == 1);
		}
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, 800);  //some time limit
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[SolveRectangles] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[SolveRectangles] optimize() failed! status:%d", status);
		return false;
	}

	//get results

	for (int i = 0; i < rect_vars.size(); i++)
	{
		bool present = round(rect_vars[i].get(GRB_DoubleAttr_X));
		if (present)
		{
			//learn the rect's W and H
			int W = -1;
			vector<int> &Ws = mask_vars_of_rects[i].first;
			for (int j = 0; j < Ws.size(); j++)
			{
				if (round(mask_vars[Ws[j]].get(GRB_DoubleAttr_X)))
				{
					W = rects[i].valid_Ws[j];
					break;
				}
			}
			if (W == -1)
			{
				AddTextCritical("[SolveRectangle] wrong W?");
				return false;
			}

			int H = -1;
			vector<int> &Hs = mask_vars_of_rects[i].second;
			for (int j = 0; j < Hs.size(); j++)
			{
				if (round(mask_vars[Hs[j]].get(GRB_DoubleAttr_X)))
				{
					H = rects[i].valid_Hs[j];
					break;
				}
			}
			if (H == -1)
			{
				AddTextCritical("[SolveRectangle] wrong W?");
				return false;
			}

			cout << "solved rect#" << i << "  at:" << rects[i].offset << " W:" << W << " H:" << H << endl;
		}
	}
	
	return true;
}

void FloorplanSpace::Floorplan::SetCamera(bool set_pickup_matrix, double *pick_x, double *pick_y,double *w, double *h)
{
	Vec2f BB[3] = { Vec2f(0) };  //min-max in x, y, z	

	//BB is simply the domain
	if (m_domain)
	{
		BB[0] = Vec2f(0, m_domain_W);
		BB[1] = Vec2f(0, m_domain_H);
		BB[2] = Vec2f(0, 0);
	}
	
	Vec3f BB_min(BB[0][0], BB[1][0], BB[2][0]);
	Vec3f BB_max(BB[0][1], BB[1][1], BB[2][1]);
	Vec3f BB_volume = BB_max - BB_min;

	//set g_camera according to BB
	g_camera->m_default_distance = 2.8f* MAX3(BB_volume.x, BB_volume.y, BB_volume.z);
	g_camera->m_center = (BB_min + BB_max) * 0.5;
	g_camera->m_nearPlane = 0.1f*BB_volume.length();
	g_camera->m_farPlane = g_camera->m_nearPlane + 5.0f*BB_volume.length();

	//apply model-view managed by g_camera
	g_camera->ApplyMODELView();

	//do viewport:
	glViewport(g_viewport[0], g_viewport[1], g_viewport[2], g_viewport[3]);

	//do projection matrix:
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	//set pickup matrix?
	if (set_pickup_matrix)
	{
		GLint vp[4] = { g_viewport[0], g_viewport[1], g_viewport[2], g_viewport[3] };
		gluPickMatrix(*pick_x, *pick_y, *w, *h, vp);
	}

	GLfloat aspect = (float)g_viewport[2] / (float)g_viewport[3];
	gluPerspective(g_camera->m_fieldofview / aspect, aspect, g_camera->m_nearPlane, g_camera->m_farPlane);
}

void FloorplanSpace::Floorplan::Draw()
{
	SetCamera(false/*not set pickup matrix*/);		

	//test: draw certain tles?
	if (m_draw_template_tile >= 0 || m_draw_root_tile >=0 || m_draw_leaf_tile >= 0)
	{
		//draw just TemplateDim domain
		glColor3f(0, 0, 0);
		for (int x = 0; x < TemplateDim; x++)
		{
			for (int y = 0; y < TemplateDim; y++)
			{
				glBegin(GL_LINE_LOOP);
				glVertex3f(x, y, 0);
				glVertex3f(x + 1, y, 0);
				glVertex3f(x + 1, y + 1, 0);
				glVertex3f(x, y + 1, 0);
				glEnd();
			}
		}		

		Tile *T = NULL;
		if (m_draw_template_tile >= 0)
			T = &m_template_tiles[m_draw_template_tile];
		else if (m_draw_root_tile >= 0)
			T = &m_root_tiles[m_draw_root_tile];
		else if (m_draw_leaf_tile >= 0)
			T = &m_leaf_tiles[m_draw_leaf_tile];

		for (int i = 0; i < T->cells.size(); i++)
		{
			int x = T->cells[i].x;
			int y = T->cells[i].y;
			glColor3f(0.5, 0.5, 0.5);
			glBegin(GL_POLYGON);
			glVertex3f(x, y, 0);
			glVertex3f(x + 1, y, 0);
			glVertex3f(x + 1, y + 1, 0);
			glVertex3f(x, y + 1, 0);
			glEnd();
		}
	}
	else if (m_draw_root_tile >= 0)
	{
		//draw just TemplateDim domain
		glColor3f(0, 0, 0);
		for (int x = 0; x < TemplateDim; x++)
		{
			for (int y = 0; y < TemplateDim; y++)
			{
				glBegin(GL_LINE_LOOP);
				glVertex3f(x, y, 0);
				glVertex3f(x + 1, y, 0);
				glVertex3f(x + 1, y + 1, 0);
				glVertex3f(x, y + 1, 0);
				glEnd();
			}
		}

		Tile &T = m_root_tiles[m_draw_root_tile];
		for (int i = 0; i < T.cells.size(); i++)
		{
			int x = T.cells[i].x;
			int y = T.cells[i].y;
			glBegin(GL_POLYGON);
			glVertex3f(x, y, 0);
			glVertex3f(x + 1, y, 0);
			glVertex3f(x + 1, y + 1, 0);
			glVertex3f(x, y + 1, 0);
			glEnd();
		}
	}
	else
	{
		//draw whole domain
		glColor3f(0, 0, 0);
		for (int x = 0; x < m_domain_W; x++)
		{
			for (int y = 0; y < m_domain_H; y++)
			{
				glBegin(GL_LINE_LOOP);
				glVertex3f(x, y, 0);
				glVertex3f(x + 1, y, 0);
				glVertex3f(x + 1, y + 1, 0);
				glVertex3f(x, y + 1, 0);
				glEnd();
			}
		}		

		//draw solved tiles
		srand(0);
		for (int i = 0; i < m_tiles.size(); i++)
		{
			glColor3f((float)(rand() % 1000) / 1000, (float)(rand() % 1000) / 1000, (float)(rand() % 1000) / 1000);

			if (m_draw_tile >= 0 && m_draw_tile != i)
				continue;

			for (int j = 0; j < m_tiles[i].cells.size(); j++)
			{
				int x = m_tiles[i].cells[j].x;
				int y = m_tiles[i].cells[j].y;
				glBegin(GL_POLYGON);
				glVertex3f(x, y, 0);
				glVertex3f(x + 1, y, 0);
				glVertex3f(x + 1, y + 1, 0);
				glVertex3f(x, y + 1, 0);
				glEnd();
			}
		}
	}
}

////functors for ceres solver

//Laplacian functors w/ various numbers of neighbor variables:

//dist-to-line functor 2D
struct FunctorDistToLine
{
	//the line is "ax+by=c"
	FunctorDistToLine(double a, double b, double c, double weight) : A(a), B(b), C(c), Weight(weight) {}

	template <typename T> bool operator()(const T* const v, T* residual) const
	{
		residual[0] = T(Weight) * (A * v[0] + B * v[1] - C) * (A * v[0] + B * v[1] - C);
		return true;
	}

private:
	double A, B, C;
	double Weight;  //weighting of this term	
};
 
//line gradient (e.g., y/x) functor 2D
//note: assume g is not too big
struct FunctorLineGradient
{
	FunctorLineGradient(double g, double weight) : G(g), Weight(weight) {}

	template <typename T> bool operator()(const T* const v0, const T* const v1, T* residual) const
	{
		residual[0] = T(Weight) * ((v1[1] - v0[1]) - G * (v1[0] - v0[0])) * ((v1[1] - v0[1]) - G * (v1[0] - v0[0]));
		return true;
	}

private:
	double G;
	double Weight;
};
//line gradient inverse (e.g., x/y) functor 2D
//note: assume g is not too big
struct FunctorLineGradientInverse
{
	FunctorLineGradientInverse(double i, double weight) : I(i), Weight(weight) {}

	template <typename T> bool operator()(const T* const v0, const T* const v1, T* residual) const
	{
		residual[0] = T(Weight) * ((v1[0] - v0[0]) - I * (v1[1] - v0[1])) * ((v1[0] - v0[0]) - I * (v1[1] - v0[1]));
		return true;
	}

private:
	double I;
	double Weight;
};

bool FloorplanSpace::FloorplanIP::SolveOld(Mesh * mesh, int solutions_to_get)
{
	//predefine l_min, l_max for the allowed range of edge length
	//
	//every vertex has pos vars (x,y).
	//
	//for every interior edge, it is one of the following state (each is a Boolean var):
	//one of the allowed edge dirs, D0,D1,...,
	//for each allowed edge dir, create a aux continuous var d.
	//collapsed, C.
	//gone, G.
	//It has to satisfy the following constraints:
	//for every allowed edge dir (P,Q) and its orthogonal (CCW) dir (A,B):
	//Ax0 + By0 >= d - BIG(1-Di) and Ax0 + By0 <= d + BIG(1-Di),
	//Ax1 + By1 >= d - BIG(1-Di) and Ax1 + By1 <= d + BIG(1-Di),
	//(Px1+Qy1) - (Px0+Qy0) >= l_min - BIG(1-Di) and (Px1+Qy1) - (Px0+Qy0) <= l_max + BIG(1-Di)
	//
	//also the collapsed constraint:
	//x1-x0 >= -BIG(1-C), x1-x0 <= BIG(1-C), y1-y0 >= -BIG(1-C), y1-y0 <= BIG(1-C)  
	//
	//boundary edge constraints:
	//for every boundary edge w/ line equation Ax+By = D
	//Ax0 + By0 >= D and Ax0 + By0 <= D,
	//Ax1 + By1 >= D and Ax1 + By1 <= D,
	//(also boundary vertex movement constraints along the edge dirs):

	const float BIG = mesh->m_BB_volume.length() * 2;
	const float EPSILON = 1e-4;
	const float l_min = mesh->m_avg_edge_len * 0.5;
	const float l_max = mesh->m_avg_edge_len * 1.5;
	const float corner_angle_min = D2R(45);  //min-max of allowed corner angles
	const float corner_angle_max = D2R(135);

	vector<Halfedge*> edges = mesh->GetFulledges();
	vector<Halfedge*> interior_edges;
	for (int i = 0; i < edges.size(); i++)
	{
		if (!edges[i]->Border2())
			interior_edges.push_back(edges[i]);
	}
	cout << "[Solve] edges:" << edges.size() << " interior_edges:" << interior_edges.size() << endl;
	
	vector<Vertex*> vertices = mesh->GetVertices();

	//re-save ori vertex poss
	m_ori_poss.clear();
	for (int i = 0; i < vertices.size(); i++)
	{
		m_ori_poss.push_back(vertices[i]->pos);
	}

	//collect allowed edge angles:
	//angle is 0~2PI angle from postive-x CCW
	vector<float> allowed_angles;
	{
		float threshold = D2R(5);

		for (int i = 0; i < edges.size(); i++)
		{
			Halfedge *e = edges[i];
			if (e->Border2())
			{
				for (int Case = 0; Case < 4; Case++)  //e's dir, rotated in 0,90,180,270 ways
				{
					float a = vectors_angle_signless(e->Dir2(), Vec2f(1, 0));
					a += D2R(90 * Case);
					if (a > D2R(360))
						a = a - D2R(360);

					//the same to existing allowed_dirs by the threshold?
					bool found = false;
					for (int j = 0; j < allowed_angles.size(); j++)
					{
						if (abs(allowed_angles[j] - a) <= threshold)
						{
							found = true;
							break;
						}
					}
					if (!found)
					{
						cout << "allowed_angle:" << R2D(a) << endl;
						allowed_angles.push_back(a);
					}
				}
			}
		}
	}
	cout << "[ToFloorplan] allowed_angles:" << allowed_angles.size() << endl;

	vector<vector<float>> edge_allowed_angles;
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<float> as;

		//take the closed allowed angle
		float a = vectors_angle_signless(interior_edges[i]->Dir2(), Vec2f(1, 0));

		float min_diff = FLT_MAX;
		float min_angle = -1;
		for (int j = 0; j < allowed_angles.size(); j++)
		{
			float diff = AngularDist(allowed_angles[j], a);
			if (diff < min_diff)
			{
				min_diff = diff;
				min_angle = allowed_angles[j];
			}
		}

		as.push_back(min_angle);
		
		//include more allowed_angles within a threshold?
		if(false)
		{
			float threshold = D2R(30);

			for (int j = 0; j < allowed_angles.size(); j++)
			{
				if (abs(allowed_angles[j] - min_angle) > 1e-4)  //not the closed one?
				{
					float diff = AngularDist(allowed_angles[j], a);
					if (diff < threshold)
					{
						as.push_back(allowed_angles[j]);
					}
				}
			}
		}

		edge_allowed_angles.push_back(as);
	}

	GRBEnv env = GRBEnv();

	//enumerate solutions until enough or infeasible
	m_solutions.clear();
	while (m_solutions.size() < solutions_to_get)
	{
		GRBModel model = GRBModel(env);

		//state variables (edge dirs, collapsed, gone) per edge
		vector<vector<GRBVar>> edge_state_vars;
		for (int i = 0; i < interior_edges.size(); i++)
		{
			vector<GRBVar> vars;

			//edge dirs states
			vector<float> &as = edge_allowed_angles[i];
			for (int j = 0; j < as.size(); j++)
			{
				vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
			}

			//collapsed state
			vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));

			edge_state_vars.push_back(vars);
		}

		//aux continuous "d" vars per edge dir per edge
		vector<vector<GRBVar>> edge_d_vars;
		for (int i = 0; i < interior_edges.size(); i++)
		{
			vector<GRBVar> vars;
			vector<float> &as = edge_allowed_angles[i];
			for (int j = 0; j < as.size(); j++)
			{
				vars.push_back(model.addVar(-100000, 100000, 0, GRB_CONTINUOUS));
			}
			edge_d_vars.push_back(vars);
		}

		//2D pos vars (x,y) for every vertex
		vector<vector<GRBVar>> vertex_vars;
		map<int, int> vertex_map;  //key: serial, value: index in vertices
		for (int i = 0; i < vertices.size(); i++)
		{
			vertex_map[vertices[i]->serial] = i;

			vector<GRBVar> vars;
			vars.push_back(model.addVar(-100000, 100000, 0, GRB_CONTINUOUS));
			vars.push_back(model.addVar(-100000, 100000, 0, GRB_CONTINUOUS));
			vertex_vars.push_back(vars);
		}

		model.update();

		////objective function:

		//GRBQuadExpr qobj;

		////maximize active edge states:
		//for (int i = 0; i < edge_state_vars.size(); i++)
		//{
		//	vector<GRBVar> &vars = edge_state_vars[i];
		//	for (int j = 0; j < vars.size(); j++)
		//	{
		//		//don't count "collapsed" states?
		//		if (j == vars.size() - 1)
		//			continue;

		//		qobj += vars[j];
		//	}
		//}

		////maxmize dists between vertices of edges
		//const float DIST_WEIGHT = 0.001;
		//for (int i = 0; i < interior_edges.size(); i++)
		//{
		//	Halfedge *e = interior_edges[i];
		//	GRBVar x0 = vertex_vars[vertex_map[e->o->v->serial]][0];
		//	GRBVar y0 = vertex_vars[vertex_map[e->o->v->serial]][1];
		//	GRBVar x1 = vertex_vars[vertex_map[e->v->serial]][0];
		//	GRBVar y1 = vertex_vars[vertex_map[e->v->serial]][1];

		//	qobj += - DIST_WEIGHT * (x1 - x0) * (x1 - x0);
		//	qobj += - DIST_WEIGHT * (y1 - y0) * (y1 - y0);
		//}

		//model.setObjective(qobj, GRB_MAXIMIZE);

		//maximize active edge states:
		GRBLinExpr obj;
		for (int i = 0; i < edge_state_vars.size(); i++)
		{
			vector<GRBVar> &vars = edge_state_vars[i];
			for (int j = 0; j < vars.size(); j++)
			{
				//don't count "collapsed" states?
				if (j == vars.size() - 1)
					continue;

				obj += vars[j];
			}
		}
		model.setObjective(obj, GRB_MAXIMIZE);

		////constraints

		//for every edge, exactly one of the state vars is true
		for (int i = 0; i < edge_state_vars.size(); i++)
		{
			GRBLinExpr sum;
			vector<GRBVar> &vars = edge_state_vars[i];
			for (int j = 0; j < vars.size(); j++)
			{
				sum += vars[j];
			}
			//model.addConstr(sum == 1);
			model.addConstr(sum <= 1);
		}

		//per-edge constraints:
		for (int i = 0; i < edge_state_vars.size(); i++)
		{
			Halfedge *e = interior_edges[i];
			GRBVar x0 = vertex_vars[vertex_map[e->o->v->serial]][0];
			GRBVar y0 = vertex_vars[vertex_map[e->o->v->serial]][1];
			GRBVar x1 = vertex_vars[vertex_map[e->v->serial]][0];
			GRBVar y1 = vertex_vars[vertex_map[e->v->serial]][1];

			//edge dir constraint:
			vector<float> &as = edge_allowed_angles[i];
			for (int j = 0; j < as.size(); j++)
			{
				float a = as[j];

				//the dir vector, normalized
				Vec2f PQ(cos(a), sin(a));   

				//AB is PQ 90' CCW rotated
				Vec2f AB = Vec2f(PQ.y, -PQ.x).normalize();

				GRBVar Di = edge_state_vars[i][j];
				GRBVar d = edge_d_vars[i][j];

				//(the edge is lying on the line equation)
				//Ax0 + By0 >= d - BIG(1-Di) and Ax0 + By0 <= d + BIG(1-Di),
				//Ax1 + By1 >= d - BIG(1-Di) and Ax1 + By1 <= d + BIG(1-Di),
				model.addConstr(AB.x * x0 + AB.y * y0 >= d - EPSILON - BIG*(1 - Di));
				model.addConstr(AB.x * x0 + AB.y * y0 <= d + EPSILON + BIG*(1 - Di));
				model.addConstr(AB.x * x1 + AB.y * y1 >= d - EPSILON - BIG*(1 - Di));
				model.addConstr(AB.x * x1 + AB.y * y1 <= d + EPSILON + BIG*(1 - Di));

				//the range of length along the edge dir
				//(Px1+Qy1) - (Px0+Qy0) >= l_min - BIG(1-Di) and (Px1+Qy1) - (Px0+Qy0) <= l_max + BIG(1-Di)
				model.addConstr((PQ.x*x1 + PQ.y*y1) - (PQ.x*x0 + PQ.y*y0) >= l_min - BIG*(1 - Di));
				model.addConstr((PQ.x*x1 + PQ.y*y1) - (PQ.x*x0 + PQ.y*y0) <= l_max + BIG*(1 - Di));
			}

			//edge dist/collapsed constraint:
			{
				GRBVar C = edge_state_vars[i][as.size()];
				//x1 - x0 >= -BIG(1 - C), x1 - x0 <= BIG(1 - C), y1 - y0 >= -BIG(1 - C), y1 - y0 <= BIG(1 - C)
				model.addConstr(x1 - x0 >= -EPSILON - BIG*(1 - C));
				model.addConstr(x1 - x0 <= EPSILON + BIG*(1 - C));
				model.addConstr(y1 - y0 >= -EPSILON - BIG*(1 - C));
				model.addConstr(y1 - y0 <= EPSILON + BIG*(1 - C));

				//debug: disable "collapsed" states?
				//model.addConstr(C == 0);
			}
		}

		//boundary edges constraints or fixed boundary vertex constraints:
		if (g_floorplan_boundary_vertex_dist > 0)
		{
			for (int i = 0; i < edges.size(); i++)
			{
				Halfedge *e = edges[i];
				if (e->Border2())
				{
					GRBVar x0 = vertex_vars[vertex_map[e->o->v->serial]][0];
					GRBVar y0 = vertex_vars[vertex_map[e->o->v->serial]][1];
					GRBVar x1 = vertex_vars[vertex_map[e->v->serial]][0];
					GRBVar y1 = vertex_vars[vertex_map[e->v->serial]][1];

					Vec2f Dir = (e->v->Pos2D() - e->o->v->Pos2D()).normalize();
					Vec2f AB(Dir.y, -Dir.x);
					float D = AB.dot(e->v->Pos2D());

					//the two vertices must lay on the line equation
					model.addConstr(AB.x * x0 + AB.y * y0 >= D - EPSILON);
					model.addConstr(AB.x * x0 + AB.y * y0 <= D + EPSILON);
					model.addConstr(AB.x * x1 + AB.y * y1 >= D - EPSILON);
					model.addConstr(AB.x * x1 + AB.y * y1 <= D + EPSILON);

					//the vertices cannot move too far away from their current pos 
					float d0 = Dir.dot(e->o->v->Pos2D());
					model.addConstr(Dir.x * x0 + Dir.y * y0 >= d0 - g_floorplan_boundary_vertex_dist);
					model.addConstr(Dir.x * x0 + Dir.y * y0 <= d0 + g_floorplan_boundary_vertex_dist);
					float d1 = Dir.dot(e->v->Pos2D());
					model.addConstr(Dir.x * x1 + Dir.y * y1 >= d1 - g_floorplan_boundary_vertex_dist);
					model.addConstr(Dir.x * x1 + Dir.y * y1 <= d1 + g_floorplan_boundary_vertex_dist);
				}
			}
		}
		//or simply fix boundary vertices
		else
		{
			for (int i = 0; i < vertices.size(); i++)
			{
				if (vertices[i]->Border())
				{
					Vec2f pos = vertices[i]->Pos2D();

					vector<GRBVar> &vars = vertex_vars[i];
					model.addConstr(vars[0] == pos.x);
					model.addConstr(vars[1] == pos.y);
				}
			}
		}

		//avoid existing solutions constraint:
		for (int i = 0; i < m_solutions.size(); i++)
		{
			vector<vector<bool>> &sol = m_solutions[i].edge_state_vars;

			int num = 0;  //total number of Boolean variables
			GRBLinExpr eq;
			for (int j = 0; j < edge_state_vars.size(); j++)
			{
				vector<bool> &flags = sol[j];
				vector<GRBVar> &vars = edge_state_vars[j];
				for (int k = 0; k < vars.size(); k++)
				{
					if (flags[k])
						eq += vars[k];
					else
						eq += 1 - vars[k];

					num++;
				}
			}
			model.addConstr(eq <= num - 1);
		}

		//solve!
		model.getEnv().set(GRB_DoubleParam_TimeLimit, g_Gurobi_time);  //some time limit
		model.optimize();
		int status = model.get(GRB_IntAttr_Status);
		if (status == 3)
		{
			//infeasible			
			break;
		}
		else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
		{
			//some other failure
			AddTextCritical("[ToFloorplan] optimize() failed! status:%d", status);
			return false;
		}

		//get results and save to a new_solution

		FloorplanSolution new_solution;

		//vertex positions:
		for (int i = 0; i < vertices.size(); i++)
		{
			Vec2f solved_pos(vertex_vars[i][0].get(GRB_DoubleAttr_X), vertex_vars[i][1].get(GRB_DoubleAttr_X));
			new_solution.vertex_poss.push_back(solved_pos);
		}

		for (int i = 0; i < interior_edges.size(); i++)
		{
			vector<bool> flags;
			int which_flag = -999;
			vector<GRBVar> &state_vars = edge_state_vars[i];
			for (int j = 0; j < state_vars.size(); j++)
			{
				bool present = round(state_vars[j].get(GRB_DoubleAttr_X));
				flags.push_back(present);
				if (present)
					which_flag = j;
			}
			cout << "#" << i << " edge_" << interior_edges[i]->serial << ": " << which_flag << endl;
			new_solution.edge_state_vars.push_back(flags);
		}
		m_solutions.push_back(new_solution);
	}

	if (m_solutions.size() > 0)
		AddText("Solve done. got %d solutions", m_solutions.size());
	else
		AddTextCritical("Solve failed. got %d solutions", m_solutions.size());
	return (m_solutions.size() > 0);
}


bool FloorplanSpace::FloorplanIP::Solve(Mesh* mesh, int solutions_to_get)
{
	//note: assume the faces circulate CW

	vector<Halfedge*> edges = mesh->GetFulledges();
	vector<Halfedge*> interior_edges;
	map<int, int> iedge_map;  //key: interior edge serial, value: index in interior_edges
	for (int i = 0; i < edges.size(); i++)
	{
		if (!edges[i]->Border2())
		{
			iedge_map[edges[i]->serial] = interior_edges.size();
			interior_edges.push_back(edges[i]);
		}
	}
	cout << "[Solve] edges:" << edges.size() << " interior_edges:" << interior_edges.size() << endl;

	//virtual diagonal edges:

	vector<Vertex*> vertices = mesh->GetVertices();

	//re-save ori vertex poss
	m_ori_poss.clear();
	for (int i = 0; i < vertices.size(); i++)
	{
		m_ori_poss.push_back(vertices[i]->pos);
	}

	//collect assignement edge angles:
	//(angle is 0~2PI angle from postive-x CCW)
	vector<float> assignment_angles;
	{
		float threshold = D2R(5);

		for (int i = 0; i < edges.size(); i++)
		{
			Halfedge* e = edges[i];
			if (e->Border2())
			{
				for (int Case = 0; Case < 4; Case++)  //e's dir, rotated in 0,90,180,270 ways
				{
					float a = vectors_angle_signless(e->Dir2(), Vec2f(1, 0));
					a += D2R(90 * Case);
					if (a > D2R(360))
						a = a - D2R(360);

					//the same to existing allowed_dirs by the threshold?
					bool found = false;
					for (int j = 0; j < assignment_angles.size(); j++)
					{
						if (abs(assignment_angles[j] - a) <= threshold)
						{
							found = true;
							break;
						}
					}
					if (!found)
					{
						cout << "new assignment_angles:" << R2D(a) << endl;
						assignment_angles.push_back(a);
					}
				}
			}
		}
	}
	cout << "[Solve] total assignment_angles:" << assignment_angles.size() << endl;

	vector<vector<float>> edge_assignment_angles;
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<float> as;

		//take the closed allowed angle
		float a = vectors_angle_signless(interior_edges[i]->Dir2(), Vec2f(1, 0));

		float min_diff = FLT_MAX;
		float min_angle = -1;
		for (int j = 0; j < assignment_angles.size(); j++)
		{
			float diff = AngularDist(assignment_angles[j], a);
			if (diff < min_diff)
			{
				min_diff = diff;
				min_angle = assignment_angles[j];
			}
		}

		as.push_back(min_angle);

		//include more allowed_angles within a threshold?
		if (false)
		{
			float threshold = D2R(45);

			for (int j = 0; j < assignment_angles.size(); j++)
			{
				if (abs(assignment_angles[j] - min_angle) > 1e-4)  //not the closed one?
				{
					float diff = AngularDist(assignment_angles[j], a);
					if (diff < threshold)
					{
						as.push_back(assignment_angles[j]);
					}
				}
			}
		}

		edge_assignment_angles.push_back(as);
	}

	//sampled angles for edges:
	vector<float> sample_angles;
	const int SAMPLES = 16;
	for (int i = 0; i < SAMPLES; i++)
	{
		sample_angles.push_back((float)i / (float)SAMPLES * 2 * MYPI);
	}

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//for every interior edge, it must be given an angle, which is either an assignment angle, or 
	//an unassigned (sampled) angle

	vector<vector<GRBVar>> edge_vars;

	//edge dirs states
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<GRBVar> vars;

		vector<float>& a_angles = edge_assignment_angles[i];
		for (int j = 0; j < a_angles.size(); j++)
		{
			vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		for (int j = 0; j < sample_angles.size(); j++)
		{
			vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		edge_vars.push_back(vars);
	}

	//objective function: maximize # of assignment angles

	GRBLinExpr obj;
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<GRBVar>& vars = edge_vars[i];

		vector<float>& a_angles = edge_assignment_angles[i];
		for (int j = 0; j < a_angles.size(); j++)
		{
			obj += vars[j];
		}
	}
	model.setObjective(obj, GRB_MAXIMIZE);

	////constraints:

	//for every interior edge, exactly one of the angle is true
	for (int i = 0; i < edge_vars.size(); i++)
	{
		GRBLinExpr sum;

		vector<GRBVar>& vars = edge_vars[i];
		for (int j = 0; j < vars.size(); j++)
		{
			sum += vars[j];
		}

		model.addConstr(sum == 1);
	}

	//for every corner, pairs of edge angles that led to bad corner angles are prohibited
	const float angle_threshold = D2R(45);  //deviation from 90'
	vector<Face*> faces = mesh->GetFaces();
	for (int i = 0; i < faces.size(); i++)
	{
		vector<Halfedge*> es = faces[i]->Edges();  //note: assume CW in Euclidean!
		for (int j = 0; j < es.size(); j++)
		{
			Halfedge* e0 = es[j];
			Halfedge* e1 = es[j]->next;
			int e0_type = 0;  //0: e is in interior_edges, 1: eo is in inteior_edges, 2:boundary edge
			int e0_index = -1;
			if (iedge_map.count(e0->serial) > 0)
			{
				e0_type = 0;
				e0_index = iedge_map[e0->serial];
			}
			else if (iedge_map.count(e0->o->serial) > 0)
			{
				e0_type = 1;
				e0_index = iedge_map[e0->o->serial];
			}
			else
				e0_type = 2;

			int e1_type = 0;
			int e1_index = -1;
			if (iedge_map.count(e1->serial) > 0)
			{
				e1_type = 0;
				e1_index = iedge_map[e1->serial];
			}
			else if (iedge_map.count(e1->o->serial) > 0)
			{
				e1_type = 1;
				e1_index = iedge_map[e1->o->serial];
			}
			else
				e1_type = 2;

			//boundary-adjacent corner cases:
			if (e0_type != 2 && e1_type == 2)
			{
				float angle1 = vectors_angle_signless(e1->Dir2(), Vec2f(1, 0));

				for (int x = 0; x < edge_vars[e0_index].size(); x++)
				{
					float angle0 = 0;
					if (x < edge_assignment_angles[e0_index].size())
						angle0 = edge_assignment_angles[e0_index][x];
					else
						angle0 = sample_angles[x - edge_assignment_angles[e0_index].size()];
					if (e0_type == 1)
					{
						//flip by 180'
						angle0 = angle0 + MYPI;
						if (angle0 >= 2 * MYPI)
							angle0 -= 2 * MYPI;
					}

					//corner angle = 180 - angle0->angle1 signed CCW angle
					float corner_angle = D2R(180) - vectors_angle_signed(Vec2f(cos(angle0), sin(angle0)), Vec2f(cos(angle1), sin(angle1)));
					//float corner_angle = AngularDist(angle0, angle1);
					if (corner_angle < MYPI / 2 - angle_threshold || corner_angle > MYPI / 2 + angle_threshold)
					{
						//disble this e0's angle var
						model.addConstr(edge_vars[e0_index][x] == 0);
					}
				}
			}
			else if (e0_type == 2 && e1_type != 2)
			{
				float angle0 = vectors_angle_signless(e0->Dir2(), Vec2f(1, 0));

				for (int y = 0; y < edge_vars[e1_index].size(); y++)
				{
					float angle1 = 0;
					if (y < edge_assignment_angles[e1_index].size())
						angle1 = edge_assignment_angles[e1_index][y];
					else
						angle1 = sample_angles[y - edge_assignment_angles[e1_index].size()];
					if (e1_type == 1)
					{
						//flip by 180'
						angle1 = angle1 + MYPI;
						if (angle1 >= 2 * MYPI)
							angle1 -= 2 * MYPI;
					}

					//corner angle = 180 - angle0->angle1 signed CCW angle
					float corner_angle = D2R(180) - vectors_angle_signed(Vec2f(cos(angle0), sin(angle0)), Vec2f(cos(angle1), sin(angle1)));
					//float corner_angle = AngularDist(angle0, angle1);
					if (corner_angle < MYPI / 2 - angle_threshold || corner_angle > MYPI / 2 + angle_threshold)
					{
						//disble this e1's angle var
						model.addConstr(edge_vars[e1_index][y] == 0);
					}
				}
			}
			//both e0 and e1 are interior edges:
			else if (e0_type != 2 && e1_type != 2)
			{
				//now, check every pair of edge angles
				for (int x = 0; x < edge_vars[e0_index].size(); x++)
				{
					//now, check every pair of edge angles
					for (int y = 0; y < edge_vars[e1_index].size(); y++)
					{
						float angle0 = 0;
						if (x < edge_assignment_angles[e0_index].size())
							angle0 = edge_assignment_angles[e0_index][x];
						else
							angle0 = sample_angles[x - edge_assignment_angles[e0_index].size()];
						if (e0_type == 1)
						{
							//flip by 180'
							angle0 = angle0 + MYPI;
							if (angle0 >= 2 * MYPI)
								angle0 -= 2 * MYPI;
						}

						float angle1 = 0;
						if (y < edge_assignment_angles[e1_index].size())
							angle1 = edge_assignment_angles[e1_index][y];
						else
							angle1 = sample_angles[y - edge_assignment_angles[e1_index].size()];
						if (e1_type == 1)
						{
							//flip by 180'
							angle1 = angle1 + MYPI;
							if (angle1 >= 2 * MYPI)
								angle1 -= 2 * MYPI;
						}

						//cout << Vec2f(cos(D2R(90)), sin(D2R(90))) << " " << Vec2f(cos(D2R(50)), sin(D2R(50))) << " " << 
						//	R2D(vectors_angle_signed(Vec2f(cos(D2R(90)), sin(D2R(90))), Vec2f(cos(D2R(50)), sin(D2R(50)))) ) << endl;

						//corner angle = 180 - angle0->angle1 signed CCW angle
						float corner_angle = D2R(180) - vectors_angle_signed(Vec2f(cos(angle0), sin(angle0)), Vec2f(cos(angle1), sin(angle1)));

						//float corner_angle = AngularDist(angle0, angle1);

						/*if (abs(corner_angle - corner_angle2) > D2R(2))
							cout << "!!" << R2D(angle0) << "," << R2D(angle1) << "  " << R2D(corner_angle) << "," << R2D(corner_angle2) << endl;
						else
							cout << "??" << R2D(angle0) << "," << R2D(angle1) << "  " << R2D(corner_angle) << "," << R2D(corner_angle2) << endl;*/

						if (corner_angle < MYPI / 2/*90*/ - angle_threshold || corner_angle > MYPI / 2 + angle_threshold)
						{
							//disble this angle-angle var pair
							model.addConstr(edge_vars[e0_index][x] + edge_vars[e1_index][y] <= 1);
						}
					}
				}

			}
		}
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, g_gurobi_time);  //some time limit	
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[Solve] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[Solve] optimize() failed! status:%d", status);
		return false;
	}

	//get results

	float obj_value = model.get(GRB_DoubleAttr_ObjVal);

	mesh->m_selected_edges.clear();
	vector<tuple<Halfedge*, float, bool>> solved_edge_angles;  //edge-angle-assigned_flag
	int num_assigned = 0;
	for (int i = 0; i < edge_vars.size(); i++)
	{
		float angle = -1;
		bool assigned = false;

		vector<GRBVar>& vars = edge_vars[i];
		for (int j = 0; j < vars.size(); j++)
		{
			bool present = round(vars[j].get(GRB_DoubleAttr_X));
			if (present)
			{
				//gotcha!
				if (j < edge_assignment_angles[i].size())
				{
					angle = edge_assignment_angles[i][j];
					assigned = true;
					num_assigned++;
				}
				else
				{
					angle = sample_angles[j - edge_assignment_angles[i].size()];
					assigned = false;

					//select un-assigned edges for visual inspection
					mesh->SelectEdge(interior_edges[i]);
				}

				break;
			}
		}
		if (angle < 0)
			AddTextCritical("[solve] no solved angle? e_%d", i);

		solved_edge_angles.push_back(make_tuple(interior_edges[i], angle, assigned));
	}

	//debug
	for (int i = 0; i < solved_edge_angles.size(); i++)
	{
		cout << "e:" << get<0>(solved_edge_angles[i])->serial << " a:" << R2D(get<1>(solved_edge_angles[i])) << " fixed:" << get<2>(solved_edge_angles[i]) << endl;
	}

	AddText("[Solve] solve angles ok! obj:%f assigned:%d/%d", obj_value, num_assigned, interior_edges.size());

	//now, solve vertex poss by given edge angles
	SolvePoss(mesh, solved_edge_angles);

	//save to m_solutions
	{
		FloorplanSolution fs;
		for (int i = 0; i < vertices.size(); i++)
		{
			fs.vertex_poss.push_back(vertices[i]->Pos2D());
		}
		m_solutions.push_back(fs);
	}

	mesh->Refresh();

	return true;
}

bool FloorplanSpace::FloorplanIP::Solve2(Mesh* mesh, int solutions_to_get)
{
	const float BIG = mesh->m_BB_volume.length() * 2;
	const float EPSILON = 1e-4;
	const float l_min = mesh->m_avg_edge_len * 0.5;
	const float l_max = mesh->m_avg_edge_len * 1.5;

	//note: assume the faces circulate CW!

	vector<Halfedge*> edges = mesh->GetFulledges();
	vector<Halfedge*> interior_edges;
	map<int, int> iedge_map;  //key: interior edge serial, value: index in interior_edges
	for (int i = 0; i < edges.size(); i++)
	{
		if (!edges[i]->Border2())  //interior edges only!
		{
			iedge_map[edges[i]->serial] = interior_edges.size();
			interior_edges.push_back(edges[i]);
		}
	}
	cout << "[Solve] edges:" << edges.size() << " interior_edges:" << interior_edges.size() << endl;

	vector<Vertex*> vertices = mesh->GetVertices();	

	//re-save ori vertex poss
	m_ori_poss.clear();
	for (int i = 0; i < vertices.size(); i++)
	{
		m_ori_poss.push_back(vertices[i]->pos);
	}

	//collect assignement edge angles:
	//(angle is 0~2PI angle from postive-x CCW)
	vector<float> assignment_angles;
	{
		float threshold = D2R(5);

		for (int i = 0; i < edges.size(); i++)
		{
			Halfedge* e = edges[i];
			if (e->Border2())
			{
				for (int Case = 0; Case < 4; Case++)  //e's dir, rotated in 0,90,180,270 ways
				{
					float a = vectors_angle_signless(e->Dir2(), Vec2f(1, 0));
					a += D2R(90 * Case);
					if (a > D2R(360))
						a = a - D2R(360);

					//the same to existing allowed_dirs by the threshold?
					bool found = false;
					for (int j = 0; j < assignment_angles.size(); j++)
					{
						if (abs(assignment_angles[j] - a) <= threshold)
						{
							found = true;
							break;
						}
					}
					if (!found)
					{
						cout << "new assignment_angles:" << R2D(a) << endl;
						assignment_angles.push_back(a);
					}
				}
			}
		}
	}
	cout << "[Solve] total assignment_angles:" << assignment_angles.size() << endl;

	vector<vector<float>> edge_assignment_angles;
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<float> as;

		//take the closed allowed angle
		float a = vectors_angle_signless(interior_edges[i]->Dir2(), Vec2f(1, 0));

		float min_diff = FLT_MAX;
		float min_angle = -1;
		for (int j = 0; j < assignment_angles.size(); j++)
		{
			float diff = AngularDist(assignment_angles[j], a);
			if (diff < min_diff)
			{
				min_diff = diff;
				min_angle = assignment_angles[j];
			}
		}

		as.push_back(min_angle);

		//include more allowed_angles within a threshold?
		if (false)
		{
			float threshold = D2R(45);

			for (int j = 0; j < assignment_angles.size(); j++)
			{
				if (abs(assignment_angles[j] - min_angle) > 1e-4)  //not the closed one?
				{
					float diff = AngularDist(assignment_angles[j], a);
					if (diff < threshold)
					{
						as.push_back(assignment_angles[j]);
					}
				}
			}
		}

		edge_assignment_angles.push_back(as);
	}

	//sampled angles for edges:
	vector<float> sample_angles;
	const int SAMPLES = 8;
	//const int SAMPLES = 0;
	for (int i = 0; i < SAMPLES; i++)
	{
		sample_angles.push_back((float)i / (float)SAMPLES * 2 * MYPI);
	}

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//for every interior edge, its "state" is either:
	//1) one of the 0~2PI CCW "angles" (whether an "assigned" angle or an unassigned (sampled) angle), or
	//2) "collapsed"

	vector<vector<GRBVar>> edge_vars;
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<GRBVar> vars;

		//edge "angle" states:

		vector<float>& a_angles = edge_assignment_angles[i];
		for (int j = 0; j < a_angles.size(); j++)
		{
			vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		for (int j = 0; j < sample_angles.size(); j++)
		{
			vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		//"collapsed" state:
		vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));

		edge_vars.push_back(vars);
	}

	//for every vertex, 2D pos vars:
	vector<vector<GRBVar>> vertex_vars;
	map<int, int> vertex_map;  //key: serial, value: index in vertex_vars
	for (int i = 0; i < vertices.size(); i++)
	{
		vertex_map[vertices[i]->serial] = i;

		vector<GRBVar> vars;
		vars.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));  //x
		vars.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));  //y
		vertex_vars.push_back(vars);
	}

	//objective function: maximize # of assignment angles
	GRBLinExpr obj;
	for (int i = 0; i < interior_edges.size(); i++)
	{
		vector<GRBVar> &vars = edge_vars[i];
		
		for (int j = 0; j < edge_assignment_angles[i].size(); j++)
		{
			obj += vars[j];
		}
	}
	model.setObjective(obj, GRB_MAXIMIZE);

	////constraints:

	//per-edge angle-poss constraints:
	for (int i = 0; i < edge_vars.size(); i++)
	{
		Halfedge* e = interior_edges[i];
		GRBVar x0 = vertex_vars[vertex_map[e->o->v->serial]][0];
		GRBVar y0 = vertex_vars[vertex_map[e->o->v->serial]][1];
		GRBVar x1 = vertex_vars[vertex_map[e->v->serial]][0];
		GRBVar y1 = vertex_vars[vertex_map[e->v->serial]][1];

		//edge dir (assignment and sampled) constraint:
		for (int j = 0; j < edge_vars[i].size() -1; j++)  //except the last one ("collapsed" state)
		{
			float angle = 0;
			if (j < edge_assignment_angles[i].size())
				angle = edge_assignment_angles[i][j];
			else
				angle = sample_angles[j - edge_assignment_angles[i].size()];

			//the dir vector, normalized
			Vec2f PQ(cos(angle), sin(angle));

			//AB is PQ 90' CCW rotated
			Vec2f AB = Vec2f(PQ.y, -PQ.x).normalize();

			GRBVar Di = edge_vars[i][j];

			//(the edge is lying on the line equation)
			// - BIG(1-Di) - EPSILON <= (Ax1 + By1) - (Ax0 + By0) <= EPSILON + BIG(1-Di)
			model.addConstr((AB.x * x1 + AB.y * y1) - (AB.x * x0 + AB.y * y0) >= - EPSILON - BIG * (1 - Di));
			model.addConstr((AB.x * x1 + AB.y * y1) - (AB.x * x0 + AB.y * y0) <= EPSILON + BIG * (1 - Di));

			//the range of length along the edge dir
			//(Px1+Qy1) - (Px0+Qy0) >= l_min - BIG(1-Di) and (Px1+Qy1) - (Px0+Qy0) <= l_max + BIG(1-Di)
			model.addConstr((PQ.x * x1 + PQ.y * y1) - (PQ.x * x0 + PQ.y * y0) >= l_min - BIG * (1 - Di));
			model.addConstr((PQ.x * x1 + PQ.y * y1) - (PQ.x * x0 + PQ.y * y0) <= l_max + BIG * (1 - Di));
		}

		//edge collapsed constraint:
		{
			GRBVar C = edge_vars[i].back();

			//x1 - x0 >= -BIG(1 - C), x1 - x0 <= BIG(1 - C), y1 - y0 >= -BIG(1 - C), y1 - y0 <= BIG(1 - C)
			model.addConstr(x1 - x0 >= -EPSILON - BIG * (1 - C));
			model.addConstr(x1 - x0 <= EPSILON + BIG * (1 - C));
			model.addConstr(y1 - y0 >= -EPSILON - BIG * (1 - C));
			model.addConstr(y1 - y0 <= EPSILON + BIG * (1 - C));

			//debug: disable "collapsed" states?
			//model.addConstr(C == 0);
		}
	}

	//for every interior edge, exactly one of its states is true
	for (int i = 0; i < edge_vars.size(); i++)
	{
		GRBLinExpr sum;

		vector<GRBVar>& vars = edge_vars[i];
		for (int j = 0; j < vars.size(); j++)
		{
			sum += vars[j];
		}

		model.addConstr(sum == 1);
	}

	//for every corner, pairs of edge angles that led to bad corner angles are prohibited
	if (true)
	{
		const float angle_threshold = D2R(45);  //deviation from 90'
		vector<Face*> faces = mesh->GetFaces();
		for (int i = 0; i < faces.size(); i++)
		{
			vector<Halfedge*> es = faces[i]->Edges();  //note: assume CW in Euclidean!
			for (int j = 0; j < es.size(); j++)
			{
				Halfedge* e0 = es[j];
				Halfedge* e1 = es[j]->next;
				int e0_type = 0;  //0: e is in interior_edges, 1: eo is in inteior_edges, 2:boundary edge
				int e0_index = -1;
				if (iedge_map.count(e0->serial) > 0)
				{
					e0_type = 0;
					e0_index = iedge_map[e0->serial];
				}
				else if (iedge_map.count(e0->o->serial) > 0)
				{
					e0_type = 1;
					e0_index = iedge_map[e0->o->serial];
				}
				else
					e0_type = 2;

				int e1_type = 0;
				int e1_index = -1;
				if (iedge_map.count(e1->serial) > 0)
				{
					e1_type = 0;
					e1_index = iedge_map[e1->serial];
				}
				else if (iedge_map.count(e1->o->serial) > 0)
				{
					e1_type = 1;
					e1_index = iedge_map[e1->o->serial];
				}
				else
					e1_type = 2;

				//boundary-adjacent corner cases:
				if (e0_type != 2 && e1_type == 2)
				{
					float angle1 = vectors_angle_signless(e1->Dir2(), Vec2f(1, 0));

					for (int x = 0; x < edge_vars[e0_index].size() - 1; x++)
					{
						float angle0 = 0;
						if (x < edge_assignment_angles[e0_index].size())
							angle0 = edge_assignment_angles[e0_index][x];
						else
							angle0 = sample_angles[x - edge_assignment_angles[e0_index].size()];
						if (e0_type == 1)
						{
							//flip by 180'
							angle0 = angle0 + MYPI;
							if (angle0 >= 2 * MYPI)
								angle0 -= 2 * MYPI;
						}

						//corner angle = 180 - angle0->angle1 signed CCW angle
						float corner_angle = D2R(180) - vectors_angle_signed(Vec2f(cos(angle0), sin(angle0)), Vec2f(cos(angle1), sin(angle1)));

						if (corner_angle < MYPI / 2 - angle_threshold || corner_angle > MYPI / 2 + angle_threshold)
						{
							//disble this e0's angle var
							model.addConstr(edge_vars[e0_index][x] == 0);
						}
					}
				}
				else if (e0_type == 2 && e1_type != 2)
				{
					float angle0 = vectors_angle_signless(e0->Dir2(), Vec2f(1, 0));

					for (int y = 0; y < edge_vars[e1_index].size() - 1; y++)
					{
						float angle1 = 0;
						if (y < edge_assignment_angles[e1_index].size())
							angle1 = edge_assignment_angles[e1_index][y];
						else
							angle1 = sample_angles[y - edge_assignment_angles[e1_index].size()];
						if (e1_type == 1)
						{
							//flip by 180'
							angle1 = angle1 + MYPI;
							if (angle1 >= 2 * MYPI)
								angle1 -= 2 * MYPI;
						}

						//corner angle = 180 - angle0->angle1 signed CCW angle
						float corner_angle = D2R(180) - vectors_angle_signed(Vec2f(cos(angle0), sin(angle0)), Vec2f(cos(angle1), sin(angle1)));

						if (corner_angle < MYPI / 2 - angle_threshold || corner_angle > MYPI / 2 + angle_threshold)
						{
							//disble this e1's angle var
							model.addConstr(edge_vars[e1_index][y] == 0);
						}
					}
				}
				//both e0 and e1 are interior edges:
				else if (e0_type != 2 && e1_type != 2)
				{
					//now, check every pair of edge angles
					for (int x = 0; x < edge_vars[e0_index].size() - 1; x++)
					{
						//now, check every pair of edge angles
						for (int y = 0; y < edge_vars[e1_index].size() - 1; y++)
						{
							float angle0 = 0;
							if (x < edge_assignment_angles[e0_index].size())
								angle0 = edge_assignment_angles[e0_index][x];
							else
								angle0 = sample_angles[x - edge_assignment_angles[e0_index].size()];
							if (e0_type == 1)
							{
								//flip by 180'
								angle0 = angle0 + MYPI;
								if (angle0 >= 2 * MYPI)
									angle0 -= 2 * MYPI;
							}

							float angle1 = 0;
							if (y < edge_assignment_angles[e1_index].size())
								angle1 = edge_assignment_angles[e1_index][y];
							else
								angle1 = sample_angles[y - edge_assignment_angles[e1_index].size()];
							if (e1_type == 1)
							{
								//flip by 180'
								angle1 = angle1 + MYPI;
								if (angle1 >= 2 * MYPI)
									angle1 -= 2 * MYPI;
							}

							//corner angle = 180 - angle0->angle1 signed CCW angle
							float corner_angle = D2R(180) - vectors_angle_signed(Vec2f(cos(angle0), sin(angle0)), Vec2f(cos(angle1), sin(angle1)));

							if (corner_angle < MYPI / 2/*90*/ - angle_threshold || corner_angle > MYPI / 2 + angle_threshold)
							{
								//disble this angle-angle var pair
								model.addConstr(edge_vars[e0_index][x] + edge_vars[e1_index][y] <= 1);
							}
						}
					}

				}
			}
		}
	}

	//boundary edges constraints or fixed boundary vertex constraints:
	if (g_floorplan_boundary_vertex_dist > 0)
	{
		for (int i = 0; i < edges.size(); i++)
		{
			Halfedge* e = edges[i];
			if (e->Border2())
			{
				GRBVar x0 = vertex_vars[vertex_map[e->o->v->serial]][0];
				GRBVar y0 = vertex_vars[vertex_map[e->o->v->serial]][1];
				GRBVar x1 = vertex_vars[vertex_map[e->v->serial]][0];
				GRBVar y1 = vertex_vars[vertex_map[e->v->serial]][1];

				Vec2f Dir = (e->v->Pos2D() - e->o->v->Pos2D()).normalize();
				Vec2f AB(Dir.y, -Dir.x);
				float D = AB.dot(e->v->Pos2D());

				//the two vertices must lay on the line equation
				model.addConstr(AB.x * x0 + AB.y * y0 >= D - EPSILON);
				model.addConstr(AB.x * x0 + AB.y * y0 <= D + EPSILON);
				model.addConstr(AB.x * x1 + AB.y * y1 >= D - EPSILON);
				model.addConstr(AB.x * x1 + AB.y * y1 <= D + EPSILON);

				//the vertices cannot move too far away from their current pos 
				float d0 = Dir.dot(e->o->v->Pos2D());
				model.addConstr(Dir.x * x0 + Dir.y * y0 >= d0 - g_floorplan_boundary_vertex_dist);
				model.addConstr(Dir.x * x0 + Dir.y * y0 <= d0 + g_floorplan_boundary_vertex_dist);
				float d1 = Dir.dot(e->v->Pos2D());
				model.addConstr(Dir.x * x1 + Dir.y * y1 >= d1 - g_floorplan_boundary_vertex_dist);
				model.addConstr(Dir.x * x1 + Dir.y * y1 <= d1 + g_floorplan_boundary_vertex_dist);
			}
		}
	}
	//or simply fix boundary vertices
	else
	{
		for (int i = 0; i < vertices.size(); i++)
		{
			if (vertices[i]->Border())
			{
				Vec2f pos = vertices[i]->Pos2D();

				vector<GRBVar>& vars = vertex_vars[i];
				model.addConstr(vars[0] == pos.x);
				model.addConstr(vars[1] == pos.y);
			}
		}
	}

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, g_gurobi_time);  //some time limit	
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[Solve] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[Solve] optimize() failed! status:%d", status);
		return false;
	}

	//get results

	float obj_value = model.get(GRB_DoubleAttr_ObjVal);

	mesh->m_selected_edges.clear();
	vector<tuple<Halfedge*, float, bool>> solved_edge_angles;  //edge-angle-assigned_flag
	int num_assigned = 0;
	for (int i = 0; i < edge_vars.size(); i++)
	{
		float angle = -1;
		bool assigned = false;

		vector<GRBVar>& vars = edge_vars[i];
		for (int j = 0; j < vars.size() - 1; j++)
		{
			bool present = round(vars[j].get(GRB_DoubleAttr_X));
			if (present)
			{
				//gotcha!
				if (j < edge_assignment_angles[i].size())
				{
					angle = edge_assignment_angles[i][j];
					assigned = true;
					num_assigned++;
				}
				else
				{
					angle = sample_angles[j - edge_assignment_angles[i].size()];
					assigned = false;

					//select un-assigned edges for visual inspection
					mesh->SelectEdge(interior_edges[i]);
				}
				break;
			}
		}		

		solved_edge_angles.push_back(make_tuple(interior_edges[i], angle, assigned));
	}

	//debug
	/*for (int i = 0; i < solved_edge_angles.size(); i++)
	{
		cout << "e:" << get<0>(solved_edge_angles[i])->serial << " a:" << R2D(get<1>(solved_edge_angles[i])) << " fixed:" << get<2>(solved_edge_angles[i]) << endl;
	}*/

	//save solved edge dirs and vertex poss to m_solutions
	{
		FloorplanSolution fs;

		for (int i = 0; i < vertex_vars.size(); i++)
		{
			float X = vertex_vars[i][0].get(GRB_DoubleAttr_X);
			float Y = vertex_vars[i][1].get(GRB_DoubleAttr_X);
			fs.vertex_poss.push_back(Vec2f(X, Y));
		}

		m_solutions.push_back(fs);
	}

	AddText("[Solve] ok! obj:%f assigned:%d/%d", obj_value, num_assigned, interior_edges.size());

	//now, solve vertex poss by given edge angles
	//SolvePoss(mesh, solved_edge_angles);
	//save to m_solutions
	/*{
		FloorplanSolution fs;
		for (int i = 0; i < vertices.size(); i++)
		{
			fs.vertex_poss.push_back(vertices[i]->Pos2D());
		}
		m_solutions.push_back(fs);
	}*/
	
	mesh->Refresh();

	return true;
}

bool FloorplanSpace::FloorplanIP::SolvePoss(Mesh* mesh, vector<tuple<Halfedge*, float, bool>>& edge_constraints)
{
	ceres::Problem problem;

	//2 double per 2D vertex poss

	vector<Vertex*> vertices = mesh->GetVertices();
	double* vars = new double[vertices.size() * 2];
	map<int, int> vertices_map;  //need a map of vertex serial to indices in vertices
	for (int i = 0; i < vertices.size(); i++)
	{
		vars[i * 2] = vertices[i]->pos.x;
		vars[i * 2 + 1] = vertices[i]->pos.z;
		vertices_map[vertices[i]->serial] = i;
	}

	//create CostFunctions

	const float BIG_WEIGHT = 1000;

	//for every boundary vertex, create two dist-to-line functors (from its two adj boundary edges)

	vector<vector<Halfedge*>> borders;
	mesh->FindBorders(borders);
	if (borders.size() != 1)
	{
		AddTextCritical("[SolvePoss] mesh not single-border?");
		return false;
	}
	vector<Halfedge*>& border = borders[0];

	ceres::CostFunction** functor_dists = new ceres::CostFunction * [border.size() * 2];  //two functor per boundary vertex

	vector<Vertex*> boundary_vertices;  //collect boundary vertices along the border (i.e., every e->o->v)
	for (int i = 0; i < border.size(); i++)
	{
		Halfedge* e_prev = border[(i + border.size() - 1) % border.size()];
		Halfedge* e = border[i];
		boundary_vertices.push_back(e_prev->v);

		//line equations of e and e_prev:
		Vec3f line0 = lineFromPoints(Vec2f(e_prev->o->v->pos[0], e_prev->o->v->pos[2]), Vec2f(e_prev->v->pos[0], e_prev->v->pos[2]));
		Vec3f line1 = lineFromPoints(Vec2f(e->o->v->pos[0], e->o->v->pos[2]), Vec2f(e->v->pos[0], e->v->pos[2]));
		functor_dists[i * 2] = new ceres::AutoDiffCostFunction<FunctorDistToLine, 1/*residual size*/, 2/*2D vertex*/>(
			new FunctorDistToLine(line0.x, line0.y, line0.z, BIG_WEIGHT));
		functor_dists[i * 2 + 1] = new ceres::AutoDiffCostFunction<FunctorDistToLine, 1/*residual size*/, 2/*2D vertex*/>(
			new FunctorDistToLine(line1.x, line1.y, line1.z, BIG_WEIGHT));
	}

	//for every edge in edge_constraints, create a gradient / inverse gradient functor

	ceres::CostFunction** functor_gradients = new ceres::CostFunction * [edge_constraints.size()];

	for (int i = 0; i < edge_constraints.size(); i++)
	{
		float angle = get<1>(edge_constraints[i]);
		float weight = 1;
		if (get<2>(edge_constraints[i]))  //hard flag?
			//weight = BIG_WEIGHT;
			weight = 1;
		else
			weight = 0.1;

		//-45~45 and 135~225: conventionl gradient (y/x)
		if ((angle >= D2R(315) && angle <= D2R(360)) || (angle >= 0 && angle <= D2R(45)) ||
			(angle >= D2R(135) && angle <= D2R(225)))
		{
			functor_gradients[i] = new ceres::AutoDiffCostFunction<FunctorLineGradient, 1/*residual size*/, 2/*2D vertex*/, 2/*2D vertex*/>(
				new FunctorLineGradient(tan(angle), weight));
		}
		//otherwise, inverse gradient (x/y)
		else
		{
			functor_gradients[i] = new ceres::AutoDiffCostFunction<FunctorLineGradientInverse, 1/*residual size*/, 2/*2D vertex*/, 2/*2D vertex*/>(
				new FunctorLineGradientInverse(1 / tan(angle), weight));
		}
	}

	//AddResidualBlocks:

	//two dist-to-line functors per boundary vertex
	for (int i = 0; i < boundary_vertices.size(); i++)
	{
		int v_index = vertices_map[boundary_vertices[i]->serial];
		double* V = &vars[v_index * 2];

		problem.AddResidualBlock(functor_dists[i * 2], NULL, V);
		problem.AddResidualBlock(functor_dists[i * 2 + 1], NULL, V);
	}

	//one line gradient / inverse gradient functor per interior edge
	for (int i = 0; i < edge_constraints.size(); i++)
	{
		Halfedge* e = get<0>(edge_constraints[i]);
		
		int v0_index = vertices_map[e->o->v->serial];
		double* V0 = &vars[v0_index * 2];
		int v1_index = vertices_map[e->v->serial];
		double* V1 = &vars[v1_index * 2];

		problem.AddResidualBlock(functor_gradients[i], NULL, V0, V1);
	}

	// Run the solver!
	ceres::Solver::Options options;
	//options.linear_solver_type = ceres::DENSE_SCHUR;  //by experiments this is faster
	options.max_num_iterations = 1000;
	options.minimizer_progress_to_stdout = true;
	//silent	
	options.logging_type = ceres::SILENT;
	options.minimizer_progress_to_stdout = false;

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	cout << summary.BriefReport() << endl;

	//copy results back

	for (int i = 0; i < vertices.size(); i++)
	{
		vertices[i]->pos = Vec3f(vars[i * 2], 0, vars[i * 2 + 1]);
	}

	return true;
}

bool FloorplanSpace::FloorplanIP::Apply(Mesh *mesh, int sol_index, float ratio)
{
	if (sol_index < 0 || sol_index >= m_solutions.size())
	{
		return false;
	}

	FloorplanSolution &sol = m_solutions[sol_index];

	vector<Vertex*> vertices = mesh->GetVertices();
	if (vertices.size() != sol.vertex_poss.size())
	{
		cout << "[Apply] vertex sizes mismatch?" << endl;
		return false;
	}

	//just apply new vertex poss?
	if (m_ori_poss.size() != vertices.size() || ratio >= 1 || ratio < 0)
	{
		for (int i = 0; i < vertices.size(); i++)
		{
			if (abs(vertices[i]->pos.y) <= 1e-6)
			{
				vertices[i]->pos.x = sol.vertex_poss[i].x;
				vertices[i]->pos.z = sol.vertex_poss[i].y;
			}
			else if (abs(vertices[i]->pos.z) <= 1e-6)
			{
				vertices[i]->pos.x = sol.vertex_poss[i].x;
				vertices[i]->pos.y = sol.vertex_poss[i].y;
			}
			else if (abs(vertices[i]->pos.x) <= 1e-6)
			{
				vertices[i]->pos.y = sol.vertex_poss[i].x;
				vertices[i]->pos.z = sol.vertex_poss[i].y;
			}
		}
	}
	//do interpolation
	else
	{
		for (int i = 0; i < vertices.size(); i++)
		{
			Vec3f ori_pos = m_ori_poss[i];
			Vec2f new_pos = sol.vertex_poss[i];

			if (abs(vertices[i]->pos.y) <= 1e-6)
			{
				vertices[i]->pos.x = ori_pos.x * (1-ratio) + new_pos[0] * ratio;
				vertices[i]->pos.z = ori_pos.z * (1 - ratio) + new_pos[1] * ratio;
			}
			else if (abs(vertices[i]->pos.z) <= 1e-6)
			{
				vertices[i]->pos.x = ori_pos.x * (1 - ratio) + new_pos[0] * ratio;
				vertices[i]->pos.y = ori_pos.y * (1 - ratio) + new_pos[1] * ratio;
			}
			else if (abs(vertices[i]->pos.x) <= 1e-6)
			{
				vertices[i]->pos.y = ori_pos.y * (1 - ratio) + new_pos[0] * ratio;
				vertices[i]->pos.z = ori_pos.z * (1 - ratio) + new_pos[1] * ratio;
			}
		}
	}

	//select some edges for visual debugs:
	if(sol.edge_state_vars.size()>0)
	{
		mesh->m_selected_edges.clear();

		vector<Halfedge*> edges = mesh->GetFulledges();
		vector<Halfedge*> interior_edges;
		for (int i = 0; i < edges.size(); i++)
		{
			if (!edges[i]->Border2())
				interior_edges.push_back(edges[i]);
		}

		//select "collapsed" and no-state edges?
		for (int i = 0; i < sol.edge_state_vars.size(); i++)
		{
			vector<bool> &flags = sol.edge_state_vars[i];
			int state = -1;
			for (int j = 0; j < flags.size(); j++)
			{
				if (flags[j])
				{
					state = j;
					break;
				}
			}
			if (state == flags.size()-1/*collapsed*/ || state == -1/*no-state*/)
			{
				mesh->SelectEdge(interior_edges[i]);
			}
		}
	}

	mesh->Refresh();

	return true;
}

bool FloorplanSpace::FloorplanIP::SolveGrids(Mesh* mesh, float step_size)
{
	//first, separate boundary edges into groups of semi-linear "segments". each group is belong to a common grid
	const float corner_angle_threshold = D2R(30);
	vector<vector<Halfedge*>> borders;
	mesh->FindBorders(borders);
	vector<vector<Halfedge*>> segments;
	{
		for (int i = 0; i < borders.size(); i++)
		{
			vector<Halfedge*>& border = borders[i];

			//collect border edges pointing to a corner:
			map<int, bool> corner_edges;			
			for (int j = 0; j < border.size(); j++)
			{
				Halfedge* e0 = border[j];
				Halfedge* e1 = border[(j + 1) % border.size()];
				float a0 = vectors_angle_signless(e0->Dir2(), Vec2f(1, 0));
				float a1 = vectors_angle_signless(e1->Dir2(), Vec2f(1, 0));
				float diff = AngularDist(a0, a1);
				if (diff >= D2R(corner_angle_threshold))
				{
					//gotcha!
					corner_edges[j] = true;
				}
			}

			if (corner_edges.empty())
			{
				//the whole border is one segment...
				segments.push_back(border);
			}
			else
			{
				//e_begin is the next of a coner edge
				int e_begin = ((*corner_edges.begin()).first + 1) % border.size();

				//let's trace segments starting at e_begin
				vector<Halfedge*> cur_segment;
				for (int offset = 0; offset < border.size(); offset++)
				{
					int i = (e_begin + offset) % border.size();

					cur_segment.push_back(border[i]);
					
					//meet corner?
					if (corner_edges.count(i) > 0)
					{
						//save cur segment and restart
						segments.push_back(cur_segment);
						cur_segment.clear();
					}
				}
			}			
		}
		cout << "segments:" << segments.size() << endl;
	}

	//for each segment, find its avg 2D angle
	//then cluster segments by its 4-way rotational angle 
	vector<vector<int>> segment_clusters;  //as indices in segments
	vector<float> segment_cluster_angles;  //angle of the clusters
	for (int i = 0; i < segments.size(); i++)
	{		
		float avg_angle = 0;
		vector<Halfedge*>& segment = segments[i];
		for (int j = 0; j < segment.size(); j++)
		{
			float angle = vectors_angle_signless(segment[j]->Dir2(), Vec2f(1, 0));
			//turn to the 0~90 angle (i.e., 4-way rot-symmetry)
			while (true)
			{
				if (angle <= D2R(90))
					break;
				angle -= D2R(90);
			}

			avg_angle += angle;
		}
		avg_angle /= (float)segment.size();

		int existing_cluster = -1;  //-1: nope
		if (!segment_clusters.empty())
		{
			//the same as any existing cluster?
			for (int j = 0; j < segment_cluster_angles.size(); j++)
			{
				if (abs(segment_cluster_angles[j] - avg_angle) <= D2R(5))
				{
					existing_cluster = j;
					break;						
				}
			}			
		}

		if (existing_cluster >= 0)
		{
			segment_clusters[existing_cluster].push_back(i);			
		}
		else
		{
			//create a new cluster / angle
			vector<int> cluster;
			cluster.push_back(i);
			segment_clusters.push_back(cluster);
			segment_cluster_angles.push_back(avg_angle);			
		}
	}
	cout << "segment_clusters:" << segment_clusters.size() << endl;

	//all "corner" boundary vertices:
	vector<Vertex*> corners_all;
	{
		for (int i = 0; i < borders.size(); i++)
		{
			vector<Halfedge*>& border = borders[i];
			for (int j = 0; j < border.size(); j++)
			{
				Halfedge* e0 = border[j];
				Halfedge* e1 = border[(j + 1) % border.size()];
				float a0 = vectors_angle_signless(e0->Dir2(), Vec2f(1, 0));
				float a1 = vectors_angle_signless(e1->Dir2(), Vec2f(1, 0));
				float diff = AngularDist(a0, a1);
				if (diff >= D2R(45))
				{
					corners_all.push_back(e0->v);					
				}
			}
		}
	}

	//the domain polygon for PointInsidePolygon() tests - assume it is the corners_all
	vector<Vec2f> domain_polygon;
	for (int i = 0; i < corners_all.size(); i++)
	{
		domain_polygon.push_back(corners_all[i]->Pos2D());
	}

	////let's collect grid points (of possibly multiple grids)

	m_grid_points.clear();

	map<int, vector<int>> boundary_vertex_points;  //key: boundary vertex serial, value: its points (indices in m_grid_points) 

	//for each cluster, solve for a grid (its direction is given and we just solve for the origin (x,y))
	GRBEnv env = GRBEnv();
	for (int c = 0; c < segment_clusters.size(); c++)
	{
		//"corner" boundary vertices of this cluster:
		vector<Vertex*> corners_this;
		{
			map<int, Vertex*> cmap;
			for (int i = 0; i < segment_clusters[c].size(); i++)
			{
				vector<Halfedge*>& segment = segments[segment_clusters[c][i]];
				cmap[segment.front()->o->v->serial] = segment.front()->o->v;
				cmap[segment.back()->v->serial] = segment.back()->v;
			}

			for (map<int, Vertex*>::iterator itr = cmap.begin(); itr != cmap.end(); itr++)
				corners_this.push_back((*itr).second);
		}		

		GRBModel model = GRBModel(env);

		vector<GRBVar> ori_pos;
		ori_pos.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));  //x
		ori_pos.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));  //y

		//"x-y" coordinates for each corner
		vector<vector<GRBVar>> corner_coordiates;
		for (int i = 0; i < corners_this.size(); i++)
		{
			vector<GRBVar> coord;
			coord.push_back(model.addVar(-1000, 1000, 0, GRB_INTEGER));
			coord.push_back(model.addVar(-1000, 1000, 0, GRB_INTEGER));
			corner_coordiates.push_back(coord);
		}

		model.update();

		//objective function:

		GRBQuadExpr qobj;

		//for every corner of this cluster, it should be matched by the grid as much as possible
		float grid_angle = segment_cluster_angles[c];
		Vec2f axis1(cos(grid_angle), sin(grid_angle));
		Vec2f axis2(cos(grid_angle + D2R(90)), sin(grid_angle + D2R(90)));
		for (int i = 0; i < corners_this.size(); i++)
		{
			Vertex* v = corners_this[i];
			Vec2f pos = v->Pos2D();
			vector<GRBVar>& coord = corner_coordiates[i];
			GRBLinExpr eq_X = ori_pos[0] + axis1[0] * step_size * coord[0] + axis2[0] * step_size * coord[1];
			GRBLinExpr eq_Y = ori_pos[1] + axis1[1] * step_size * coord[0] + axis2[1] * step_size * coord[1];
			qobj += (pos.x - eq_X) * (pos.x - eq_X);
			qobj += (pos.y - eq_Y) * (pos.y - eq_Y);
		}
		model.setObjective(qobj);  //to minimize

		model.getEnv().set(GRB_DoubleParam_TimeLimit, 1);  //some time limit
		model.optimize();
		int status = model.get(GRB_IntAttr_Status);
		if (status == 3)
		{
			//infeasible
			AddTextCritical("[SolveGrids] the problem is infeasible.");
			return false;
		}
		else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
		{
			//some other failure
			AddTextCritical("[SolveGrids] optimize() failed! status:%d", status);
			return false;
		}

		//get results

		Vec2f Ori(ori_pos[0].get(GRB_DoubleAttr_X), ori_pos[1].get(GRB_DoubleAttr_X));
		cout << "cluster#" << c << " ori:" << Ori << endl;
		for (int i = 0; i < corners_this.size(); i++)
		{
			Vec2f pos = corners_this[i]->Pos2D();
			vector<GRBVar>& coord = corner_coordiates[i];
			Vec2f Coord(round(coord[0].get(GRB_DoubleAttr_X)), round(coord[1].get(GRB_DoubleAttr_X)));
			Vec2f Pos = Ori + Coord[0] * step_size * axis1 + Coord[1] * step_size * axis2;			
		}

		//let's collect grid points of this grid:
		//first find a discrete BB of all corners, slightly enlarged
		//then include all points that are within the boundary,
		//and also include points that are on the same grid-lines as the segment edges
		Vec2i BB_min(INT_MAX, INT_MAX), BB_max(-INT_MAX, -INT_MAX);
		for (int i = 0; i < corners_all.size(); i++)
		{
			Vec2f pos = corners_all[i]->Pos2D();

			//round to a discrete pos:
			Vec2f vec = pos - Ori;
			//project vec onto axis1 and axis2 and divide by step_size to get coordinate:
			int coord0 = round(vec.dot(axis1) / step_size);
			int coord1 = round(vec.dot(axis2) / step_size);
			if (coord0 < BB_min[0])
				BB_min[0] = coord0;
			if (coord0 > BB_max[0])
				BB_max[0] = coord0;
			if (coord1 < BB_min[1])
				BB_min[1] = coord1;
			if (coord1 > BB_max[1])
				BB_max[1] = coord1;
		}
		
		//for every segment of this cluster, find its corresponding horizontal or vertical grid-line
		vector<pair<Vec2i, Vec2i>> segment_grid_lines;
		for (int i = 0; i < segment_clusters[c].size(); i++)
		{
			vector<Halfedge*>& segment = segments[segment_clusters[c][i]];
			Vec2f p0 = segment.front()->o->v->Pos2D();
			Vec2f p1= segment.back()->v->Pos2D();

			//project to grid coordinates
			Vec2i c0(round((p0 - Ori).dot(axis1) / step_size), round((p0 - Ori).dot(axis2) / step_size));
			Vec2i c1(round((p1 - Ori).dot(axis1) / step_size), round((p1 - Ori).dot(axis2) / step_size));
			segment_grid_lines.push_back(make_pair(c0, c1));
		}

		//now, test every point within the BB to see if it is inside the domain
		vector<int> points_here;
		for (int coord0 = BB_min[0]; coord0 <= BB_max[0]; coord0++)
		{
			for (int coord1 = BB_min[1]; coord1 <= BB_max[1]; coord1++)
			{
				Vec2f pos = Ori + coord0 * step_size * axis1 + coord1 * step_size * axis2;

				//is this point on a grid-line along some boundary edges?
				bool found = false;
				for (int ii = 0; ii < segment_grid_lines.size(); ii++)
				{
					Vec2i c0 = segment_grid_lines[ii].first;
					Vec2i c1 = segment_grid_lines[ii].second;
					if( (coord0 == c0[0] && coord0 == c1[0] &&
						coord1 >= MIN2(c0[1], c1[1]) && coord1 <= MAX2(c0[1], c1[1]))/*horizontal?*/ ||
						(coord1 == c0[1] && coord1 == c1[1] &&
							coord0 >= MIN2(c0[0], c1[0]) && coord0 <= MAX2(c0[0], c1[0])/*vertical?*/)
					  )
					{
						//gocha
						points_here.push_back(m_grid_points.size());
						m_grid_points.push_back(make_tuple(pos, c, true/*boundary*/));
						found = true;
						break;
					}
				}
				if (!found)
				{
					if (point_in_polygon_2d(pos, domain_polygon))
					{
						points_here.push_back(m_grid_points.size());
						m_grid_points.push_back(make_tuple(pos, c, false/*interior*/));
					}
				}
			}
		}

		//for every boundary vertex here, find its associated point as the closest one
		for (int ii = 0; ii < corners_this.size(); ii++)
		{
			int closest = -1;
			float closest_dist = FLT_MAX;
			for (int jj = 0; jj < points_here.size(); jj++)
			{
				float dist = (corners_this[ii]->Pos2D() - get<0>(m_grid_points[points_here[jj]])).length();
				if (dist < closest_dist)
				{
					closest = points_here[jj];
					closest_dist = dist;
				}
			}
			if (closest != -1)
			{
				boundary_vertex_points[corners_this[ii]->serial].push_back(closest);
			}
		}

		//hack: just one cluster?
		//break;
	}

	//merge points associated with the same boundary vertices
	{
		vector<tuple<Vec2f, int, bool>> new_points;  //new points to add
		map<int, bool> points_to_remove;
		for (map<int, vector<int>>::iterator itr = boundary_vertex_points.begin(); itr != boundary_vertex_points.end(); itr++)
		{
			vector<int>& ps = (*itr).second;
			if (ps.size() > 1)
			{
				//merge the associated points to a single one
				tuple<Vec2f, int, bool> new_point;
				//calculate avg pos:
				Vec2f avg(0);
				for (int ii = 0; ii < ps.size(); ii++)
				{
					avg += get<0>(m_grid_points[ps[ii]]);

					points_to_remove[ps[ii]] = true;  //record this point to remove
				}
				avg /= (float)ps.size();
				get<0>(new_point) = avg;
				get<1>(new_point) = -1;  //cluster no. -1 for "merged"
				get<2>(new_point) = true;  //yes it is boundary point

				new_points.push_back(new_point);
			}
		}

		vector<tuple<Vec2f, int, bool>> new_grid_points;  //grid_points after merging
		for (int ii = 0; ii < m_grid_points.size(); ii++)
		{
			if (points_to_remove.count(ii) == 0)  //not a point to remove?
				new_grid_points.push_back(m_grid_points[ii]);  
		}
		for (int ii = 0; ii < new_points.size(); ii++)  //add new points
			new_grid_points.push_back(new_points[ii]);

		cout << "grid points merge: " << m_grid_points.size() << "->" << new_grid_points.size() << endl;
		m_grid_points = new_grid_points;
	}
	
	AddText("[SolveGrids] done. grid_points:%d", m_grid_points.size());
	return true;
}

//HalfEdge for SolveMesh()
struct HE
{
	int p0;  //p0->p1 (as indices in points)
	int p1;
	int o;  //opposite index
	float angle;  //0~2PI angle
	bool is_boundary;
};

//recursively trace faces
//vector<int> cur_face: current tracing face, as indices in corners
//faces_all: all unique faces traced so far. the corner indices are sorted from small to large
void TraceFace(vector<int> &cur_face, vector<HE> &halfedges, vector<tuple<int, int, int, float>> &corners, map<int,vector<int>> &corner_map, 
	Vec2i deg_range, vector<vector<int>> &faces_all)
{
	//have enough degree for a face? - try add to faces_all
	bool formed = false;
	if (cur_face.size() >= deg_range[0])
	{
		//closed? (the last corner's e_next is the first corner's e_prev)
		tuple<int, int, int, float>& first_corner = corners[cur_face.front()]; 
		tuple<int, int, int, float>& last_corner = corners[cur_face.back()];
		if (get<1>(last_corner) == get<0>(first_corner))
		{
			formed = true;

			//new unique face?
			bool new_face = true;
			for (int i = 0; i < faces_all.size(); i++)
			{
				vector<int>& f = faces_all[i];
				if (cur_face.size() != f.size())
					continue;
				//compare at every possible starting pos
				for (int start = 0; start < f.size(); start++)
				{
					bool same = true;
					for (int j = 0; j < f.size(); j++)
					{
						if (cur_face[j] != f[(start + j) % f.size()])
						{
							same = false;
							break;
						}
					}
					if (same)
					{
						new_face = false;
						break;
					}
				}
			}
			if (new_face)
			{
				faces_all.push_back(cur_face);
			}
		}		
	}
	if (formed)
		return;  //have formed a closed face so no need to grow more
	
	//can still enlarge by one corner?
	if (cur_face.size() < deg_range[1])
	{
		//enlarge with every possible next corner:
		
		//grow to every corners at e_next->v
		int e_next = get<1>(corners[cur_face.back()]);
		vector<int>& next_corners = corner_map[halfedges[e_next].p1];
		for (int i = 0; i < next_corners.size(); i++)
		{
			//is this corner's e_prev = e_next?
			int ep = get<0>(corners[next_corners[i]]);
			if (ep == e_next)
			{
				//gotcha!
				vector<int> new_face = cur_face;
				new_face.push_back(next_corners[i]);  //add a new corner to the cur_face
				TraceFace(new_face, halfedges, corners, corner_map, deg_range, faces_all);
			}
		}		
	}
}

bool FloorplanSpace::FloorplanIP::SolveMesh(vector<tuple<Vec2f, int, bool>>& points, float ideal_edge_len)
{
	//1) enumerate edge candidates between two points w/ len within the range	
	//weight each edge w/ its length r.s.t. the ideal length
	//
	////objective function: maximize weighted # of edges
	//
	////constraints:
	//
	//all "boundary" points must be active (exactly one of its config is active)
	//otherwise, at most one of its config is active
	//
	//for every point, enuemrate all its valid edge configurations. a valid config must be:
	//1) valence >= 2 and <= upper_bound
	//2) no bad corner angles w.r.t. 90'
	//
	//"non-overlapping" constraint:
	//two too-close edges cannot appear at the same time
	//
	//"singly-connected constraint":
	//except for boundary points, for a point to be active, at least one of its edges has "asscending depth"
	//(TODO)

	//allowed edge len range:
	const Vec2f len_range(ideal_edge_len * 0.5, ideal_edge_len * 1.3);
	
	//allowed corner angle range (radian):
	const Vec2f corner_angle_range(D2R(45), D2R(135));

	//allowed face degrees
	const Vec2i face_degree_range(3, 5);

	//enumerate half-edges
	
	vector<HE> halfedges;
	//incoming edges of each points in CCW (angle-asscending) order
	vector<vector<int/*half-edge index*/>> point_edges(points.size());
	for (int i = 0; i < points.size()-1; i++)
	{
		Vec2f p0 = get<0>(points[i]);
		for (int j = i + 1; j < points.size(); j++)
		{
			Vec2f p1 = get<0>(points[j]);

			float len = (p1 - p0).length();
			if (len < len_range[0] || len > len_range[1])
				continue;  //len out-of-range

			bool is_boundary = get<2>(points[i]) && get<2>(points[j]);  //e is boundary if both points are boundary

			//create e and eo:
			int e_index = halfedges.size();
			int eo_index = halfedges.size() + 1;
			{
				HE e;
				e.p0 = i;
				e.p1 = j;
				e.o = eo_index;
				e.angle = vectors_angle_signless((p1 - p0).normalize(), Vec2f(1, 0));
				e.is_boundary = is_boundary;
				halfedges.push_back(e);
			}
			{
				HE eo;
				eo.p0 = j;
				eo.p1 = i;
				eo.o = e_index;
				eo.angle = vectors_angle_signless((p0 - p1).normalize(), Vec2f(1, 0));
				eo.is_boundary = is_boundary;
				halfedges.push_back(eo);
			}
						
			//insert to point_edges records of p0 and p1
			for (int Case = 0; Case < 2; Case++)  //p0->p1 or p1->p0
			{
				//angle w.r.t. p0 or p1
				float angle = 0;
				vector<int>* incoming_edges = NULL;
				int ei = -1;
				if (Case == 0)
				{
					angle = halfedges[e_index].angle;
					incoming_edges = &point_edges[j];
					ei = e_index;
				}
				else
				{
					angle = halfedges[eo_index].angle;
					incoming_edges = &point_edges[i];
					ei = eo_index;
				}

				//find the right pos to insert
				int insert_pos = -1;
				for (int k = 0; k < (*incoming_edges).size(); k++)
				{
					int ii = (*incoming_edges)[k];
					if (halfedges[ii].angle >= angle)
					{
						//insert here						
						insert_pos = k;
						break;
					}
				}				
				if (insert_pos >= 0)
					(*incoming_edges).insert((*incoming_edges).begin() + insert_pos, ei);
				else
					(*incoming_edges).push_back(ei);
			}
		}
	}
	
	//enumerate all possible corners (each is two edges at a point)
	vector<tuple<int, int, int, float>> corners;   //{e_prev index, e_next index, point index, CCW angle}
	map<int, vector<int>> corner_map;  //key: point index, value: # of corners at the point
	for (int i = 0; i < points.size(); i++)
	{
		vector<int>& incoming_edges = point_edges[i];

		//check every possible pair:
		for (int x = 0; x < incoming_edges.size(); x++)
		{
			for (int y = 0; y < incoming_edges.size(); y++)
			{
				if (x == y)
					continue;

				int e_prev = incoming_edges[x];
				int e_next = halfedges[incoming_edges[y]].o;

				//CCW corner angle?
				float a0 = halfedges[e_prev].angle + D2R(180);  //opposite dir
				if (a0 >= D2R(360))
					a0 -= D2R(360);
				float a1 = halfedges[e_next].angle;
				float corner_angle = 0;  //CCW from a1 to a0
				if (a0 >= a1)
					corner_angle = a0 - a1;
				else
					corner_angle = D2R(360) - (a1 - a0);

				if (corner_angle >= corner_angle_range[0] && corner_angle <= corner_angle_range[1])
				{
					corner_map[i].push_back(corners.size());
					corners.push_back(make_tuple(e_prev, e_next, i, corner_angle));
				}
			}
		}
	}

	//enumerate all possible faces:
	//starting from a corner, recursively trace several consecutive corners to see if it forms a closed loop
	vector<vector<int>> faces;  //each is a list of corner indices
	for (int i = 0; i < corners.size(); i++)
	{
		vector<int> cur_face;
		cur_face.push_back(i);  //start with a new face w/ a corner
		TraceFace(cur_face, halfedges, corners, corner_map, face_degree_range, faces);
	}
	cout << "[SolveMesh] points:" << points.size() << " halfedges:" << halfedges.size() << 
		" corners:" << corners.size() << " faces:" << faces.size() << endl;

	//build lists of face at every half-edge
	vector<vector<int>> edge_faces(halfedges.size());
	for (int i = 0; i < faces.size(); i++)
	{
		vector<int>& f = faces[i];  //as corner indices
		for(int j = 0; j < f.size(); j++)
		{
			int he = get<1>(corners[f[j]]);  //the "next" half-edge at this corner
			edge_faces[he].push_back(i);
		}		
	}

	//score every face
	vector<float> face_penalities;  //penalities of every face (0 = perfect)
	for (int i = 0; i < faces.size(); i++)
	{
		vector<int> &f = faces[i];

		float penalty = 0;

		//non-quad penalty?
		if (f.size() != 4)
			penalty += 100;

		//sum of angle deviations
		float perfect_angle = MYPI - (2 * MYPI / f.size());
		for (int j = 0; j < f.size(); j++)
		{
			//angle at this corner
			float angle = get<3>(corners[f[j]]);
			float deviation = abs(perfect_angle - angle);
			penalty += deviation;
		}		

		face_penalities.push_back(penalty);
	}

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	//one Boolean for every face
	vector<GRBVar> face_vars;
	for (int i = 0; i < faces.size(); i++)
	{
		face_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
	}

	model.update();

	//objective function:

	//minimize # of faces
	/*GRBLinExpr obj;
	for (int i = 0; i < face_vars.size(); i++)
	{
		obj += face_vars[i];
	}
	model.setObjective(obj);*/

	//minimize sum of face penalities
	GRBLinExpr obj;
	for (int i = 0; i < face_vars.size(); i++)
	{
		obj += face_penalities[i] * face_vars[i];
	}
	model.setObjective(obj);

	////constraints:

	//for every boundary half-edge, exactly one of its adj faces is true
	//else for interior half-edge, at most one of its faces is true
	for (int i = 0; i < edge_faces.size(); i++)
	{
		if (edge_faces[i].empty())
			continue;

		GRBLinExpr sum;
		for (int j = 0; j < edge_faces[i].size(); j++)
		{
			sum += face_vars[edge_faces[i][j]];
		}

		//boundary half-edge? (no faces, both vertices are boundary, and of the same boundary cluster)
		bool boundary = false;
		int eo_index = halfedges[i].o;		
		if (halfedges[i].is_boundary && edge_faces[eo_index].empty())
		{	
			HE& he = halfedges[i];
			if( (get<1>(points[he.p0]) == get<1>(points[he.p1])) ||
				get<1>(points[he.p0]) == -1 || get<1>(points[he.p1]) == -1)
			{
				//cout << "bhe " << edge_faces[i].size() << endl;
				model.addConstr(sum == 1);
				boundary = true;
			}						
		}
		if(!boundary)
		{
			model.addConstr(sum <= 1);
		}
	}

	//for every e/eo pair, # of faces equal (either 0 or 1, simutenously)
	if(true)
	{
		map<int, bool> visited;  //visited halfedge indices
		for (int ei = 0; ei < edge_faces.size(); ei++)
		{
			if (visited.count(ei) > 0)
				continue;

			int eo_index = halfedges[ei].o;
			visited[ei] = true;
			visited[eo_index] = true;

			if (edge_faces[ei].size() > 0 && edge_faces[eo_index].size() > 0)
			{
				GRBLinExpr sum0;
				for (int ii = 0; ii < edge_faces[ei].size(); ii++)
				{
					sum0 += face_vars[edge_faces[ei][ii]];
				}
				GRBLinExpr sum1;
				for (int ii = 0; ii < edge_faces[eo_index].size(); ii++)
				{
					sum1 += face_vars[edge_faces[eo_index][ii]];
				}

				model.addConstr(sum0 == sum1);
			}
		}
	}

	//overlapping faces cannot appear at the same time
	if (false)
	{
		//calculate the polygon of every face. a little shrinked for intersection test
		vector< vector<Vec2f> > face_polys;
		for (int i = 0; i < faces.size(); i++)
		{
			vector<Vec2f> poly;
			Vec2f avg(0);
			for (int j = 0; j < faces[i].size(); j++)
			{
				int pi = get<2>(corners[faces[i][j]]);
				poly.push_back(get<0>(points[pi]));
				avg += get<0>(points[pi]);
			}
			avg /= (float)faces[i].size();

			for (int j = 0; j < poly.size(); j++)
			{
				poly[j] = avg + (poly[j] - avg) * 0.995;
			}

			face_polys.push_back(poly);
		}

		int num_overlaps = 0;  //for debug
		for (int i = 0; i < faces.size() - 1; i++)
		{
			vector<Vec2f>& poly0 = face_polys[i];

			for (int j = i + 1; j < faces.size(); j++)
			{
				vector<Vec2f>& poly1 = face_polys[j];

				if (polygon_polygon_intersection_2d(poly0, poly1))
				{
					GRBVar& f0_var = face_vars[i];
					GRBVar& f1_var = face_vars[j];
					model.addConstr(f0_var + f1_var <= 1);
					num_overlaps++;
				}
			}
		}
		cout << "num_overlaps:" << num_overlaps << endl;
	}

	//solve!

	//model.getEnv().set(GRB_DoubleParam_TimeLimit, 200);  //some time limit
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		AddTextCritical("[SolveConnectivity] the problem is infeasible.");
		return false;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		AddTextCritical("[SolveConnectivity] optimize() failed! status:%d", status);
		return false;
	}

	//get results
	m_solved_faces.clear();
	for (int i = 0; i < faces.size(); i++)
	{
		bool present = round(face_vars[i].get(GRB_DoubleAttr_X));
		if (present)
		{
			vector<Vec2f> face;
			for (int j = 0; j < faces[i].size(); j++)
			{
				int pi = get<2>(corners[faces[i][j]]);
				face.push_back(get<0>(points[pi]));				
			}
			m_solved_faces.push_back(face);
		}
	}

	AddText("[SolveConnectivity] done! faces:%d", m_solved_faces.size());
	return true;
}

void FloorplanSpace::FloorplanIP::Draw()
{
	float point_size = g_draw_vertex_size / g_camera->m_scale;
	float edge_size = g_draw_edge_size / g_camera->m_scale;

	//draw points
	for (int i = 0; i < m_grid_points.size(); i++)
	{
		Vec2f p = get<0>(m_grid_points[i]);
		int c = get<1>(m_grid_points[i]);  //# cluster
		
		Vec3f color(0);
		if (c == 0)
			color = Vec3f(1, 0, 0);
		else if (c == 1)
			color = Vec3f(0, 1, 0);
		else if (c == 2)
			color = Vec3f(0, 0, 1);
		else if (c == 3)
			color = Vec3f(1, 1, 0);
		else if (c == 4)
			color = Vec3f(1, 0, 1);
		else if (c == 5)
			color = Vec3f(0, 1, 1);
		else if (c == 6)
			color = Vec3f(1, 1, 0.5);

		DrawSphere(p.x, 0, p.y, color, point_size);
	}

	//draw polys
	for (int i = 0; i < m_solved_faces.size(); i++)
	{
		vector<Vec2f>& face = m_solved_faces[i];
		for (int j = 0; j < face.size(); j++)
		{
			Vec3f v0(face[j].x, 0, face[j].y);
			Vec3f v1(face[(j + 1) % face.size()].x, 0, face[(j + 1) % face.size()].y);
			DrawBar(v0, v1, g_draw_edge_size * 0.5, g_draw_edge_size * 0.5, Vec3f(0.4), Vec3f(0.4));
		}
	}
}