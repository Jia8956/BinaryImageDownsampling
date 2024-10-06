//#include "Basic.h"
#include <windows.h>
#include <vector>
#include <iostream>
#include <map>
#include <unordered_map>
#include <queue>
#include <omp.h>  //OpenMP
#include "gurobi_c++.h"
#include "lodepng.h"
#include "ILMBase.h"
#include "downsampling.h"

#define MIN2(a,b) (((a) < (b))?(a):(b))
#define MAX2(a,b) (((a) > (b))?(a):(b))
#define MAX3(a,b,c) ( ( (MAX2(a,b)) > (c) ) ? (MAX2(a,b)) : (c) )
#define MIN3(a,b,c) ( ( (MIN2(a,b)) > (c) ) ? (MIN2(a,b)) : (c) )

extern int g_ds_bigpixel_width;
extern int g_ds_bigpixel_height;
extern bool g_ds_save;
extern bool g_ds_print_debug;
extern float g_ds_land_bias_ratio;
extern int g_ds_land_weight;
extern Vec2i g_ds_input_size;  //w,h
extern vector<bool> g_ds_input;  //saved input binary mask
extern Vec2i g_ds_output_size;  //w,h
extern vector<bool> g_ds_output;  //saved output binary mask
extern vector<int> g_ds_output_components;
extern vector<int> g_ds_input_components;
extern bool g_ds_save_components_to_file;
extern vector< vector<tuple<Vec2i, int, int>> > g_ds_output_boundaries;  //pos, type, solved dist
extern bool g_ds_alternative_offsets;
extern int g_ds_h_offset;
extern int g_ds_v_offset;
extern int g_ds_png_treshold;
extern bool g_ds_local_constraint;
extern int g_ds_neighobrhood_offset;

using namespace DSSpace;

//label topology (i.e., lands (=islands) and water (= sea and lakes) of a given mask
//we will mark every pixel of the mask as:
//{type (false=land, true=water), index (# of island; # of sea (0) or lakes (>0)) }
//num_islands = "h0", num_waters = "h1" for betti numbers
//labels (result): key: row-major index, value: pos and label 
//Eulers: Euler characterstics of all components (lands then waters)
void DSSpace::LabelTopology(int width, int height, bool* mask, unordered_map<int, pair<Vec2i,pair<bool, int>>> &labels,
	int & num_lands, int &num_waters, vector<tuple<int,int,int>> &Eulers, vector<bool> &is_boundary_flags)
{
	num_lands = 0;  //current # of land
	num_waters = 0;  //current # of water

	labels.clear();

	while (true)
	{
		bool filled = false;  //filled a new component at this pass?

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int index = y * width + x;

				if (labels.count(index) > 0)
					continue;  //already labelled

				bool is_land = mask[index];

				//let's start a flooding of either a land or a water from here
				unordered_map<int, Vec2i> visited;
				visited[index] = Vec2i(x, y);
				queue<Vec2i> Q;
				Q.push(Vec2i(x, y));
				while (!Q.empty())
				{
					Vec2i p = Q.front();
					Q.pop();

					//try to flood to adjacent land or water
					//land component: consider all the 8 neighbors in a 3x3 mask
					//water component: only 4-neighbors
					for (int dir = 0; dir < 8; dir++)
					{
						if (!is_land)
						{
							if (dir == 0 || dir == 2 || dir == 4 || dir == 6)
								continue;
						}

						Vec2i p2;
						if (dir == 0)
							p2 = p + Vec2i(-1, -1);
						else if (dir == 1)
							p2 = p + Vec2i(0, -1);
						else if (dir == 2)
							p2 = p + Vec2i(1, -1);
						else if (dir == 3)
							p2 = p + Vec2i(1, 0);
						else if (dir == 4)
							p2 = p + Vec2i(1, 1);
						else if (dir == 5)
							p2 = p + Vec2i(0, 1);
						else if (dir == 6)
							p2 = p + Vec2i(-1, 1);
						else if (dir == 7)
							p2 = p + Vec2i(-1, 0);

						int p2_index = p2.y * width + p2.x;

						//is p2 within and not visited yet?
						if (p2.x >= 0 && p2.x < width && p2.y >= 0 && p2.y < height &&
							visited.count(p2_index) == 0)
						{
							bool to_add = false;  //both land or both water?
							if ((is_land && mask[p2_index] == true) ||
								(!is_land && mask[p2_index] == false))
								to_add = true;

							if (to_add)
							{
								//go here!
								visited[p2_index] = p2;
								Q.push(p2);
							}
						}
					}
				}

				//label this mass
				pair<bool, int> type_;
				if (is_land)
				{
					type_.first = false;
					type_.second = num_lands;
					num_lands++;
				}
				else
				{
					type_.first = true;
					type_.second = num_waters;
					num_waters++;
				}

				for (unordered_map<int, Vec2i>::iterator itr = visited.begin(); itr != visited.end(); itr++)
				{
					labels[(*itr).first] = make_pair((*itr).second, type_);
				}

				filled = true;
			}
		}

		if (!filled)
		{
			//all pixels are visited/filled, done!
			break;
		}
	}

	//find Euler characteristics of components

	//for each components:
	vector<unordered_map<int, bool>> faces_all(num_lands + num_waters);  //width * height
	vector<unordered_map<int, bool>> vertices_all(num_lands + num_waters);  //(width+1) * (height + 1)
	vector<unordered_map<int, bool>> hedges_all(num_lands + num_waters);  //width * (height + 1)
	vector<unordered_map<int, bool>> vedges_all(num_lands + num_waters);  //(width+1) * height
	for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
	{
		int c = -1;
		if (!(*itr).second.second.first)
			c = (*itr).second.second.second;
		else
			c = num_lands + (*itr).second.second.second;

		Vec2i XY = (*itr).second.first;  //face's pos

		faces_all[c][XY.y * width + XY.x] = true;

		//mark this faces' adjacent vertices and edges
		vertices_all[c][XY.y * (width + 1) + XY.x] = true;
		vertices_all[c][XY.y * (width + 1) + (XY.x + 1)] = true;
		vertices_all[c][(XY.y + 1) * (width + 1) + (XY.x + 1)] = true;
		vertices_all[c][(XY.y + 1) * (width + 1) + XY.x] = true;

		//mark adjcent hedges and vedges
		hedges_all[c][XY.y * width + XY.x] = true;
		hedges_all[c][(XY.y + 1) * width + XY.x] = true;
		vedges_all[c][XY.y * (width + 1) + XY.x] = true;
		vedges_all[c][XY.y * (width + 1) + (XY.x + 1)] = true;
	}
	//calculate the Euler characteristics:
	for (int c = 0; c < num_lands + num_waters; c++)
	{
		int Euler = vertices_all[c].size() - (hedges_all[c].size() + vedges_all[c].size()) + faces_all[c].size();

		if (g_ds_print_debug)
		{
			if (c < num_lands)
				cout << "land  ";
			else
				cout << "water ";
			cout << "C#" << c << " V:" << vertices_all[c].size() << " E:" << hedges_all[c].size() + vedges_all[c].size() <<
				" F:" << faces_all[c].size() << " Euler:" << Euler << endl;
		}

		Eulers.push_back(make_tuple(vertices_all[c].size(), hedges_all[c].size() + vedges_all[c].size(), 
			faces_all[c].size()));
		//Eulers.push_back(Euler);
	}

	//find is_boundary flags for every component:
	is_boundary_flags.clear();
	is_boundary_flags.resize(num_lands + num_waters, false);
	for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
	{
		//boundary pixel?
		Vec2i xy = (*itr).second.first;
		if (xy.x == 0 || xy.x == width - 1 || xy.y == 0 || xy.y == height - 1)
		{
			int c = -1;
			if (!(*itr).second.second.first)
				c = (*itr).second.second.second;
			else
				c = num_lands + (*itr).second.second.second;

			is_boundary_flags[c] = true;
		}
	}

	if (g_ds_print_debug)
	{
		cout << "[LabelTopology] #labels:" << labels.size() << " num_lands:" << num_lands << " num_waters" << num_waters << endl;
	}
}

void DSSpace::LabelTopology2(int width, int height, bool* mask, unordered_map<int, pair<Vec2i, pair<bool, int>>>& labels,
	int& num_lands, int& num_waters)
{
	num_lands = 0;  //current # of land
	num_waters = 0;  //current # of water

	labels.clear();

	while (true)
	{
		bool filled = false;  //filled a new component at this pass?

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int index = y * width + x;

				if (labels.count(index) > 0)
					continue;  //already labelled

				bool is_land = mask[index];

				//let's start a flooding of either a land or a water from here
				unordered_map<int, Vec2i> visited;
				visited[index] = Vec2i(x, y);
				queue<Vec2i> Q;
				Q.push(Vec2i(x, y));
				while (!Q.empty())
				{
					Vec2i p = Q.front();
					Q.pop();

					//try to flood to adjacent land or water
					//land component: consider all the 8 neighbors in a 3x3 mask
					//water component: only 4-neighbors
					for (int dir = 0; dir < 8; dir++)
					{
						if (!is_land)
						{
							if (dir == 0 || dir == 2 || dir == 4 || dir == 6)
								continue;
						}

						Vec2i p2;
						if (dir == 0)
							p2 = p + Vec2i(-1, -1);
						else if (dir == 1)
							p2 = p + Vec2i(0, -1);
						else if (dir == 2)
							p2 = p + Vec2i(1, -1);
						else if (dir == 3)
							p2 = p + Vec2i(1, 0);
						else if (dir == 4)
							p2 = p + Vec2i(1, 1);
						else if (dir == 5)
							p2 = p + Vec2i(0, 1);
						else if (dir == 6)
							p2 = p + Vec2i(-1, 1);
						else if (dir == 7)
							p2 = p + Vec2i(-1, 0);

						int p2_index = p2.y * width + p2.x;

						//is p2 within and not visited yet?
						if (p2.x >= 0 && p2.x < width && p2.y >= 0 && p2.y < height &&
							visited.count(p2_index) == 0)
						{
							bool to_add = false;  //both land or both water?
							if ((is_land && mask[p2_index] == true) ||
								(!is_land && mask[p2_index] == false))
								to_add = true;

							if (to_add)
							{
								//go here!
								visited[p2_index] = p2;
								Q.push(p2);
							}
						}
					}
				}

				//label this mass
				pair<bool, int> type_;
				if (is_land)
				{
					type_.first = false;
					type_.second = num_lands;
					num_lands++;
				}
				else
				{
					type_.first = true;
					type_.second = num_waters;
					num_waters++;
				}

				for (unordered_map<int, Vec2i>::iterator itr = visited.begin(); itr != visited.end(); itr++)
				{
					labels[(*itr).first] = make_pair((*itr).second, type_);
				}

				filled = true;
			}
		}

		if (!filled)
		{
			//all pixels are visited/filled, done!
			break;
		}
	}

	if (g_ds_print_debug)
	{
		cout << "[LabelTopology] #labels:" << labels.size() << " num_lands:" << num_lands << " num_waters" << num_waters << endl;
	}
}

bool DSSpace::EnumerateBoundaries(int width, int height, unordered_map<int, pair<Vec2i, pair<bool, int>>>& labels,
	int num_lands, int num_waters, vector<pair<Vec2i,int>> &boundaries)
{
	//check every half-edge (horizontal and vertical), look for component boundary edges
	//we assume land is at e->f and water is at eo->f

	//(x,y) is the left or lower vertex position of a horizontal or vertical edge, respectively

	//first collect F0 (he->f) and F1 (he->o-f) pairs
	vector<pair< pair<Vec2i, pair<bool, int>>, pair<Vec2i, pair<bool, int>> >> F0F1s;
	{
		//check horizontal edges (w/ faces below and above):
		for (int vy = 1; vy < height; vy++)
		{
			for (int vx = 0; vx < width; vx++)
			{
				//case0: he is (vx,vy) to (vx+1,vy), face above is land / face below is water
				//case1: he is (vx+1,vy) to (vx,vy), face below is land / face above is water
				for (int Case = 0; Case < 2; Case++)
				{
					int X0 = 0, Y0 = 0;  //he->f
					int X1 = 0, Y1 = 0;  //he->o->f

					if (Case == 0)
					{
						X0 = vx;
						Y0 = vy;
						X1 = vx;
						Y1 = vy - 1;
					}
					else
					{
						X0 = vx;
						Y0 = vy - 1;
						X1 = vx;
						Y1 = vy;
					}

					//locate the two face pixels:

					pair<Vec2i, pair<bool, int>> F0;
					if (labels.count(Y0 * width + X0) == 0)
					{
						cout << "cannot find " << X0 << "," << Y0 << endl;
						return false;
					}
					F0 = labels[Y0 * width + X0];

					pair<Vec2i, pair<bool, int>> F1;
					if (labels.count(Y1 * width + X1) == 0)
					{
						cout << "cannot find " << X1 << "," << Y1 << endl;
						return false;
					}
					F1 = labels[Y1 * width + X1];

					F0F1s.push_back(make_pair(F0, F1));
				}
			}
		}

		//check vertical edges (w/ faces left and right):
		for (int vy = 0; vy < height; vy++)
		{
			for (int vx = 1; vx < width; vx++)
			{
				//case0: he is (vx,vy) to (vx,vy + 1), face left is land / face right is water
				//case1: he is (vx,vy + 1) to (vx,vy), face right is land / face left is water
				for (int Case = 0; Case < 2; Case++)
				{
					int X0 = 0, Y0 = 0;  //he->f
					int X1 = 0, Y1 = 0;  //he->o->f

					if (Case == 0)
					{
						X0 = vx - 1;
						Y0 = vy;
						X1 = vx;
						Y1 = vy;
					}
					else
					{
						X0 = vx;
						Y0 = vy;
						X1 = vx - 1;
						Y1 = vy;
					}

					//locate the two face pixels:

					pair<Vec2i, pair<bool, int>> F0;
					if (labels.count(Y0 * width + X0) == 0)
					{
						cout << "cannot find " << X0 << "," << Y0 << endl;
						return false;
					}
					F0 = labels[Y0 * width + X0];

					pair<Vec2i, pair<bool, int>> F1;
					if (labels.count(Y1 * width + X1) == 0)
					{
						cout << "cannot find " << X1 << "," << Y1 << endl;
						return false;
					}
					F1 = labels[Y1 * width + X1];

					F0F1s.push_back(make_pair(F0, F1));
				}
			}
		}
	}

	//key: 0-based land_index * (num_components) + 0-based water_index. value: <index-index, size>
	unordered_map<int, pair<Vec2i,int>> boundaries_map; 
	for (int i = 0; i < F0F1s.size(); i++)
	{
		pair<Vec2i, pair<bool, int>>& F0 = F0F1s[i].first;
		pair<Vec2i, pair<bool, int>>& F1 = F0F1s[i].second;
		
		//is F0 a land pixel and F1 a water pixel?
		if (F0.second.first == false && F1.second.first == true)
		{
			boundaries_map[F0.second.second * (num_lands+num_waters) + F1.second.second].first = 
				Vec2i(F0.second.second, F1.second.second + num_lands /*actual component index of this water component*/);
			boundaries_map[F0.second.second * (num_lands + num_waters) + F1.second.second].second++;
		}
	}
	//to vector
	for (unordered_map<int, pair<Vec2i,int>>::iterator itr = boundaries_map.begin(); itr != boundaries_map.end(); itr++)
	{
		if(g_ds_print_debug)
			cout << (*itr).second.first << " size:" << (*itr).second.second << endl;
		boundaries.push_back((*itr).second);
	}

	return true;
}


bool DSSpace::DownsamplePng(const char* input_filename, bool calculate_error_metrics)
{
	std::vector<unsigned char> in_buffer; //the raw pixels (RGBA)
	unsigned width = 0, height = 0;
	unsigned error = lodepng::decode(in_buffer, width, height, input_filename);
	if (error)
	{
		cout << input_filename << " lodepng::decode error:" << error << " " << lodepng_error_text(error) << endl;
		return false;
	}

	if ((width % g_ds_bigpixel_width) != 0 || (height % g_ds_bigpixel_height) != 0)
	{
		cout << "error: width / height not dividable!" << endl;
		return false;
	}

	//turn the image buffer to a binary mask buffer
	bool* mask_ori = new bool[width * height];
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int offset = (y * width + x) * 4;
			Vec4<unsigned char> p(in_buffer[offset], in_buffer[offset + 1], in_buffer[offset + 2], in_buffer[offset + 3]);

			//binarilization: <threshold = black, >threshold = white
			if (p[0] >= g_ds_png_treshold && p[1] >= g_ds_png_treshold && p[2] >= g_ds_png_treshold)
			{
				mask_ori[y * width + x] = true;
			}
			else
			{
				mask_ori[y * width + x] = false;
			}
		}
	}

	//test fill holes?
	//FillHoles(width, height, 2, mask_ori);

	//test thinning?
	/*if(false)
	{
		bool* tmp = new bool[width * height];
		Thinning(width, height, mask, tmp);

		g_ds_skeleton.clear();
		g_ds_skeleton.resize(width * height);
		for (int i = 0; i < width * height; i++)
		{
			g_ds_skeleton[i] = tmp[i];
		}

		delete tmp;
	}*/

	//save input buffer to global
	g_ds_input_size = Vec2i(width, height);
	g_ds_input.clear();
	g_ds_input.resize(width * height);
	for (int i = 0; i < width * height; i++)
	{
		g_ds_input[i] = mask_ori[i];
	}

	//multuple input masks w/ alternative offsets?
	vector<bool*> masks;
	vector<Vec2i> offsets;  //corresponding {h-offset and v-offset} for every masks
	
	masks.push_back(mask_ori);  //default:
	offsets.push_back(Vec2i(0, 0));

	if (g_ds_alternative_offsets)
	{
		int h_offset_step = g_ds_bigpixel_width / 4;
		int v_offset_step = g_ds_bigpixel_height / 4;
		/*h_offset_step = 1;
		v_offset_step = 1;*/

		//create more alternative versions of mask buffers w/ horizontal and vertical shifts
		for (int h_offset = 0; h_offset <= g_ds_bigpixel_width / 2; h_offset += h_offset_step)
		{
			for (int v_offset = 0; v_offset <= g_ds_bigpixel_height / 2; v_offset += v_offset_step)
			{
				if (h_offset == 0 && v_offset == 0)
					continue;  //skip default one

				bool* mask_ = new bool[width * height];
				for (int y = 0; y < height; y++)
				{
					for (int x = 0; x < width; x++)
					{
						//corresponding pos in original buffer?
						int x_ = x - h_offset;
						int y_ = y - v_offset;
						if (x_ < 0 || y_ < 0)
						{
							mask_[y * width + x] = false;
						}
						else
						{
							mask_[y * width + x] = mask_ori[y_ * width + x_];
						}
					}
				}

				masks.push_back(mask_);
				offsets.push_back(Vec2i(h_offset, v_offset));
			}
		}
	}

	int new_width = width / g_ds_bigpixel_width;
	int new_height = height / g_ds_bigpixel_height;

	//try all mask buffers, keep the one w/ best error metrics
	float best_score = -1;
	for (int i = 0; i < masks.size(); i++)
	{
		bool* output = new bool[new_width * new_height];
		float IoU = 0, Dice = 0, Precision = 0, Recall = 0;
		if (Downsample(width, height, masks[i], g_ds_bigpixel_width, g_ds_bigpixel_height, output) == 0)
		{
			//calculate error metrics?
			//(must need to do if alternative offsets are taken)
			if(calculate_error_metrics || g_ds_alternative_offsets)
				ErrorMetrics(width, height, masks[i], new_width, new_height, output, 
					IoU, Dice, Precision, Recall);

			//best one so far? save to output
			if (g_ds_output.size() == 0 /*first one?*/ || IoU > best_score /*consider what?*/)
			{
				g_ds_output_size = Vec2i(new_width, new_height);
				g_ds_output.clear();
				g_ds_output.resize(new_width * new_height);
				for (int i = 0; i < new_width * new_height; i++)
				{
					g_ds_output[i] = output[i];
				}

				best_score = IoU;

				//record the best offset:
				g_ds_h_offset = offsets[i].x;
				g_ds_v_offset = offsets[i].y;
			}
		}
		
		delete output;
	}

	for (int i = 0; i < masks.size(); i++)
	{
		delete masks[i];
	}
	
	if (g_ds_output.size() == 0)  //failed?
	{
		return false;
	}
	else
	{
		//save output to a png file
		if (g_ds_save)
		{
			std::vector<unsigned char> out_buffer(new_width * new_height * 4);
			for (int y = 0; y < new_height; y++)
			{
				for (int x = 0; x < new_width; x++)
				{
					int offset = (y * new_width + x) * 4;

					if (g_ds_output[y * new_width + x])
					{
						//white
						out_buffer[offset] = 255;
						out_buffer[offset + 1] = 255;
						out_buffer[offset + 2] = 255;
						out_buffer[offset + 3] = 255;
					}
					else
					{
						//black
						out_buffer[offset] = 0;
						out_buffer[offset + 1] = 0;
						out_buffer[offset + 2] = 0;
						out_buffer[offset + 3] = 255;
					}
				}
			}

			string output_filename = string(input_filename) + "." + to_string(new_width) + "x" + to_string(new_height);
			output_filename += ".lw" + to_string(g_ds_land_weight);
			output_filename += ".png";

			lodepng::encode(output_filename, out_buffer, new_width, new_height);
		}

		//save component flags?
		if (g_ds_save_components_to_file)
		{
			//save input component indices to png file:
			if (g_ds_input_components.size() == width * height)
			{
				std::vector<unsigned char> buffer(width * height * 4);
				for (int y = 0; y < height; y++)
				{
					for (int x = 0; x < width; x++)
					{
						int offset = (y * width + x) * 4;

						//to color:
						int index = g_ds_input_components[y * width + x] * 5;
						Vec3i color(index % 256, (index / 256) % 256, (index / 65536) % 256);

						buffer[offset] = color[0];
						buffer[offset + 1] = color[1];
						buffer[offset + 2] = color[2];
						buffer[offset + 3] = 255;
					}
				}
				string filename = string(input_filename) + "." + to_string(new_width) + "x" + to_string(new_height);
				filename += ".input_components.png";
				lodepng::encode(filename, buffer, width, height);
			}

			//save output component indices to png file:
			if (g_ds_output_components.size() == new_width * new_height)
			{
				std::vector<unsigned char> buffer(new_width* new_height * 4);
				for (int y = 0; y < new_height; y++)
				{
					for (int x = 0; x < new_width; x++)
					{
						int offset = (y * new_width + x) * 4;

						//to color:
						int index = g_ds_output_components[y * new_width + x] * 5;
						Vec3i color(index % 256, (index / 256) % 256, (index / 65536) % 256);

						buffer[offset] = color[0];
						buffer[offset + 1] = color[1];
						buffer[offset + 2] = color[2];
						buffer[offset + 3] = 255;

					}
				}
				string filename = string(input_filename) + "." + to_string(new_width) + "x" + to_string(new_height);
				filename += ".output_components.png";
				lodepng::encode(filename, buffer, new_width, new_height);
			}
		}

		return true;
	}
}

bool DSSpace::Thinning(int width, int height, bool* input/*size = width*height*/, bool* output/*size = width*height */)
{
	//copy input to output. we will work on output buffer only
	for (int i = 0; i < width * height; i++)
	{
		output[i] = input[i];
	}

	//thinning cut-out masks
	//0:must be 0, 1:must be 1, 2:don't care. in left-to-right, bottom-to-top row order
	vector<vector<int>> masks;  
	//default orientation (upward):
	{
		masks.push_back(vector<int>{ 1, 1, 1, 2, 1, 2, 0, 0, 0 });
	}
	{
		masks.push_back(vector<int>{ 2, 1, 2, 1, 1, 0, 2, 0, 0 });
	}
	for (int orientation = 0; orientation < 3; orientation++)
	{
		//source:
		vector<int>& Mask0 = masks[0];
		vector<int>& Mask1 = masks[1];

		vector<int> mask0(9), mask1(9);  //target
		
		//CCW 90':
		if (orientation == 0)
		{
			vector<int>* Mask = NULL;
			vector<int>* mask = NULL;
			for (int Case = 0; Case < 2; Case++)
			{
				if (Case == 0)
				{
					Mask = &Mask0;
					mask = &mask0;
				}
				else
				{
					Mask = &Mask1;
					mask = &mask1;
				}

				(*mask)[6] = (*Mask)[0];
				(*mask)[3] = (*Mask)[1];
				(*mask)[0] = (*Mask)[2];
				(*mask)[7] = (*Mask)[3];
				(*mask)[4] = (*Mask)[4];
				(*mask)[1] = (*Mask)[5];
				(*mask)[8] = (*Mask)[6];
				(*mask)[5] = (*Mask)[7];
				(*mask)[2] = (*Mask)[8];
			}
		}
		//CCW 180':
		else if (orientation == 1)
		{
			vector<int>* Mask = NULL;
			vector<int>* mask = NULL;
			for (int Case = 0; Case < 2; Case++)
			{
				if (Case == 0)
				{
					Mask = &Mask0;
					mask = &mask0;
				}
				else
				{
					Mask = &Mask1;
					mask = &mask1;
				}

				(*mask)[8] = (*Mask)[0];
				(*mask)[7] = (*Mask)[1];
				(*mask)[6] = (*Mask)[2];
				(*mask)[5] = (*Mask)[3];
				(*mask)[4] = (*Mask)[4];
				(*mask)[3] = (*Mask)[5];
				(*mask)[2] = (*Mask)[6];
				(*mask)[1] = (*Mask)[7];
				(*mask)[0] = (*Mask)[8];
			}
		}
		//CCW 270':
		else if (orientation == 2)
		{
			vector<int>* Mask = NULL;
			vector<int>* mask = NULL;
			for (int Case = 0; Case < 2; Case++)
			{
				if (Case == 0)
				{
					Mask = &Mask0;
					mask = &mask0;
				}
				else
				{
					Mask = &Mask1;
					mask = &mask1;
				}

				(*mask)[2] = (*Mask)[0];
				(*mask)[5] = (*Mask)[1];
				(*mask)[8] = (*Mask)[2];
				(*mask)[1] = (*Mask)[3];
				(*mask)[4] = (*Mask)[4];
				(*mask)[7] = (*Mask)[5];
				(*mask)[0] = (*Mask)[6];
				(*mask)[3] = (*Mask)[7];
				(*mask)[6] = (*Mask)[8];
			}
		}

		masks.push_back(mask0);
		masks.push_back(mask1);
	}

	int iter_max = 1e5;
	for (int iter = 0; iter < iter_max; iter++)
	{
		//find pixels w/ 3x3 neighborhood exactly matching a mask, setting it to 0
		//do in-place refinements

		int num_modified = 0;

		//we thin by each mask
		for (int m = 0; m < 8; m++)
		{
			vector<int>& mask = masks[m];

			for (int y = 1; y < height - 1; y++)
			{
				for (int x = 1; x < width - 1; x++)
				{
					//compare:
					bool match = true;
					for (int yy = 0; yy <= 2; yy++)
					{
						for (int xx = 0; xx <= 2; xx++)
						{
							int s = (int)output[(y + yy - 1) * width + (x + xx - 1)];
							int t = mask[yy * 3 + xx];
							if (t != 2 && (s != t))
							{
								match = false;
								break;
							}
						}
					}

					if (match)
					{
						//match! set to 0:
						output[y * width + x] = 0;
						num_modified++;;
					}
				}
			}

		}

		if (num_modified == 0)
		{
			//done!
			cout << "[Thinning] done at iter#" << iter << endl;
			break;
		}
		else if (iter == iter_max - 1)
		{
			cout << "[Thinning] ?? reached iter_max " << iter_max << endl;
		}

		cout << "[Thinning] iter#" << iter << " num_modified:" << num_modified << endl;
	}

	return true;
}

bool DSSpace::FillHole(int width, int height, bool* mask)
{
	//let's analyze the topology of input buffer first
	//let's label the topology of current mask first
	unordered_map<int, pair<Vec2i, pair<bool, int>>> labels;
	int num_lands = 0;
	int num_waters = 0;
	vector<tuple<int,int,int>> Eulers;
	vector<bool> is_boundary_flags;
	LabelTopology(width, height, mask, labels, num_lands, num_waters, Eulers, is_boundary_flags);

	//find the smallest water component
	int smallest_water_size = -1;
	int smallest_index = -1;
	for (int i = num_waters; i < Eulers.size()/*same as # of components*/; i++)
	{	
		int num_pixels = get<2>(Eulers[i]);
		if (smallest_water_size <0 || num_pixels < smallest_water_size)
		{
			//new smallest!
			smallest_water_size = num_pixels;
			smallest_index = i;
		}
	}
	cout << "[FillHoles] smallest_index:" << smallest_index << " size:" << smallest_water_size << endl;
	
	//now, fill pixels of the smallest component components
	for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
	{
		Vec2i xy = (*itr).second.first;

		//the actual component index:
		int c = -1;
		if (!(*itr).second.second.first)
			c = (*itr).second.second.second;
		else
			c = num_lands + (*itr).second.second.second;

		if (c == smallest_index)
		{
			//turn this water pixel to land!
			cout << "fill " << xy << endl;
			mask[xy.y * width + xy.x] = true;
		}
	}

	return true;
}

//calculate pixel weights for land and water pixels, respectively
//total_pixels: # pixels in a big-pixel. threshold_land_count: the amount of land pixels to "just" win
//e.g., in a 16-bit pixel, threshold = 8, then 8 land pixels win 8 water pixels, but 7 land pixels shall lose 9 water pixels
bool DSSpace::CalculatePixelWeights(int total_pixels, int threshold_land_count, int &land_weight, int &water_weight)
{
	//name threshold_lands as X and threshold_waters as Y
	float X = threshold_land_count;
	float Y = total_pixels - threshold_land_count;
	//name land_weight as a and water_weight as b

	//we have:
	//X*a > Y*b
	//(X-1)*a < (Y+1)*b
	
	// Y/X * b < a < (Y+1)/(X-1) * b

	//we want to find the smallest solution, so we try b=1, b=2,... until we find a feasible solution fo a
	int a = 0;
	int b = 1;
	const int iter_max = 10000;
	for(int iter=0; iter<iter_max; iter++)
	{
		float LHS = Y / X * b;
		float RHS = (Y + 1) / (X - 1) * b;

		int lower = floor(LHS);
		int upper = ceil(RHS);

		if ((upper - 1) >= (lower +1))
		{
			//this b works!
			a = lower + 1;
			cout << "[CalculatePixelWeights]X:" << X << " Y:" << Y << " found!LHS:" << LHS << " RHS : " << RHS << " a = " << a << " b = " << b << endl;
			land_weight = a;
			water_weight = b;
			return true;
		}
		else
		{
			//nope? try bigger b
			b = b + 1;
		}
	}

	cout << "[CalculatePixelWeights] failed?" << endl;
	return false;
}

//generate a neighborhood mask according to Euclidean distance threshold
void NeighborhoodMask(float dist, vector<Vec2i>& mask)
{
	//cout << "NeighborhoodMask: ";
	for (int x = -(int)dist; x <= (int)dist; x++)
	{
		for (int y = -(int)dist; y <= (int)dist; y++)
		{
			Vec2f p(x, y);
			if (p.length() <= dist)
			{
				mask.push_back(p);
			}
		}
	}
	//cout << endl;
}

int DSSpace::DownsampleByEuler(int width, int height, bool* mask, int bigpixel_width, int bigpixel_height, bool* output)
{
	if (width % bigpixel_width != 0 || height % bigpixel_width != 0)
	{
		cout << "[Downsampling] error: width / height not dividable" << endl;
		return 1;
	}

	DWORD time_begin = timeGetTime();
	
	GRBEnv env;
	
	int cur_width = width;
	int cur_height = height;
	bool* cur_mask = new bool[cur_width * cur_height];
	memcpy(cur_mask, mask, sizeof(bool) * cur_width * cur_height);

	//let's label the topology of current mask first
	unordered_map<int, pair<Vec2i, pair<bool, int>>> labels;
	int num_lands = 0;
	int num_waters = 0;
	vector<tuple<int,int,int>> Eulers;
	vector<bool> is_boundary_flags;

	DWORD time = time_begin;

	LabelTopology(cur_width, cur_height, cur_mask, labels, num_lands, num_waters, Eulers, is_boundary_flags);

	if (g_ds_print_debug)
	{
		cout << "[Downsampling] LabelTopology time:" << timeGetTime() - time << endl;
		time = timeGetTime();
	}

	//size of the new mask:
	int new_width = cur_width / bigpixel_width;
	int new_height = cur_height / bigpixel_height;

	GRBModel model(env);

	//for every connected compoent (islands and waters), collect its overlapping big-pixels
	//and create the Boolean vars for each of its big-pixels
	struct BigPixel  //a big-pixel of a component
	{
		Vec2i pos;  //in output grid
		int score;  //"score" of this big-pixel w.s.t. the component
		GRBVar var;

		BigPixel()
		{
			pos = Vec2i(0);
			score = 0;
		}
	};
	vector<unordered_map<int, BigPixel>> components;  //big-pixels of every component (lands first then waters)
	components.resize(num_lands + num_waters);
		
	//also record the occupying potential big-pixels of each big-pixel locations 
	vector<unordered_map<int,bool>> bigpixels_map;  //row-major order. each is a map of component#
	bigpixels_map.resize(new_width * new_height);

	//"pixel-neighborhood" way:
	if (true)
	{
		const int land_weight = g_ds_land_weight;
		const int water_weight = 1;

		//get dist-based neighborhood mask
		vector<Vec2i> neighborhood;
		NeighborhoodMask(bigpixel_width / 4, neighborhood);

		for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
		{
			Vec2i xy = (*itr).second.first;

			//the actual component index:
			int c = -1;
			if (!(*itr).second.second.first)
				c = (*itr).second.second.second;
			else
				c = num_lands + (*itr).second.second.second;

			unordered_map<int, BigPixel>& component = components[c];

			int weight = 0;
			if (c < num_lands)  //land pixel?
				weight = land_weight;
			else
				weight = water_weight;

			//propagate scores within a neighborood, add to bigpixel's scores
			//assume the neighborhood is half size of a bigpixel?
			for (int ii = 0; ii < neighborhood.size(); ii++)
			{	
				Vec2i xy_ = xy + neighborhood[ii];

				//within?
				if (xy_.x >= 0 && xy_.x < width && xy_.y >= 0 && xy_.y < height)
				{
					//to residing big-pixel coordinate:
					Vec2i XY(xy_.x / bigpixel_width, xy_.y / bigpixel_height);

					//first time for this big-pixel?
					if (component.count(XY.y * new_width + XY.x) == 0)
					{
						component[XY.y * new_width + XY.x].pos = XY;
						component[XY.y * new_width + XY.x].score = weight;
						component[XY.y * new_width + XY.x].var = model.addVar(0, 1, 0, GRB_BINARY);

						bigpixels_map[XY.y * new_width + XY.x][c] = true;
					}
					else
					{
						//accumulate a score
						component[XY.y * new_width + XY.x].score += weight;
					}
				}
			}
		}
	}
	//"CalculatePixelWeights" way:
	else
	{
		//let's build components and bigpixels_map:
		int weight_inside_land = 0;
		int weight_inside_water = 0;
		int weight_outside = 0;
		int land_bias = ceil((float)(bigpixel_width * bigpixel_height) / 2) * g_ds_land_bias_ratio;
		CalculatePixelWeights(bigpixel_width * bigpixel_height,
			ceil((float)(bigpixel_width * bigpixel_height) / 2) - land_bias,
			weight_inside_land, weight_inside_water);
		for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
		{
			Vec2i xy = (*itr).second.first;

			//the actual component index:
			int c = -1;
			if (!(*itr).second.second.first)
				c = (*itr).second.second.second;
			else
				c = num_lands + (*itr).second.second.second;

			unordered_map<int, BigPixel>& component = components[c];

			//to residing big-pixel coordinate:
			Vec2i XY(xy.x / bigpixel_width, xy.y / bigpixel_height);

			int weight = 0;
			if (c < num_lands)  //land pixel?
				weight = weight_inside_land;
			else
				weight = weight_inside_water;

			//first time for this big-pixel?
			if (component.count(XY.y * new_width + XY.x) == 0)
			{
				component[XY.y * new_width + XY.x].pos = XY;
				component[XY.y * new_width + XY.x].score = weight;
				component[XY.y * new_width + XY.x].var = model.addVar(0, 1, 0, GRB_BINARY);

				bigpixels_map[XY.y * new_width + XY.x][c] = true;
			}
			else
			{
				//accumulate a score
				component[XY.y * new_width + XY.x].score += weight;
			}

		}

		//"grow" each component's big-pixels to include adjacent big-pixels?
		if (true)
		{
			for (int c = 0; c < components.size(); c++)
			{
				unordered_map<int, BigPixel>& component = components[c];

				//collect new big-pixel locations
				unordered_map<int, Vec2i> new_bigpixels;
				for (unordered_map<int, BigPixel>::iterator itr = component.begin(); itr != component.end(); itr++)
				{
					Vec2i XY_ = (*itr).second.pos;

					//try neighbors in a 3x3 promixity
					for (int Y_ = -1; Y_ <= 1; Y_++)
					{
						for (int X_ = -1; X_ <= 1; X_++)
						{
							Vec2i XY = XY_ + Vec2i(X_, Y_);

							if (XY.x < 0 || XY.x >= new_width || XY.y < 0 || XY.y >= new_height)
								continue;

							//new big-pixel?
							if (component.count(XY.y * new_width + XY.x) == 0)
							{
								new_bigpixels[XY.y * new_width + XY.x] = XY;
							}
						}
					}
				}

				//create big-pixel records for the new ones
				for (unordered_map<int, Vec2i>::iterator itr = new_bigpixels.begin(); itr != new_bigpixels.end(); itr++)
				{
					Vec2i XY = (*itr).second;

					//create a new big-pixel for the component, w/ outside_weight
					component[XY.y * new_width + XY.x].pos = XY;
					component[XY.y * new_width + XY.x].score = weight_outside;
					component[XY.y * new_width + XY.x].var = model.addVar(0, 1, 0, GRB_BINARY);

					bigpixels_map[XY.y * new_width + XY.x][c] = true;
				}
			}
		}
	}

	//to model the Euler characteristics of each component, we need to model 
	//the Boolean vars of the vertices and edges

	vector< vector</*for every vertex*/ tuple< Vec2i/*itself position*/, 
		vector<Vec2i>/*adjacent vars of big-pixels of a vertex*/, GRBVar>> > component_vertices;
	for (int i = 0; i < components.size(); i++)
	{
		vector<tuple<Vec2i, vector<Vec2i>, GRBVar>> vertices;

		//let's check every possible vertex positions (vx,vy) to find its adjacent big-pixels of this component
		//note: some vertices will NOT have adjacent active big-pixels so they are not present
		for (int vy = 0; vy <= new_height; vy++)
		{
			for (int vx = 0; vx <= new_width; vx++)
			{
				vector<Vec2i> faces;  //adjacent face (X,Y) of this vertex

				//this vertex #(vx,vy) have at most 4 neighbor faces (vx-1 or vx, vy-1 or vy)
				for (int Y = vy - 1; Y <= vy; Y++)
				{
					for (int X = vx - 1; X <= vx; X++)
					{
						if (Y < 0 || Y >= new_height || X < 0 || X >= new_width)
							continue;  //out of bound

						//does the component have this face (X,Y)?
						if (components[i].count(Y * new_width + X) > 0)
						{
							faces.push_back(Vec2i(X, Y));
						}
					}
				}

				if (faces.size() > 0)
				{
					//this vertex do has adjacent valid faces of the component. let's create it
					vertices.push_back(make_tuple(Vec2i(vx,vy), faces, model.addVar(0, 1, 0, GRB_BINARY)));
				}
			}
		}
			
		component_vertices.push_back(vertices);
	}

	vector< vector<tuple<Vec2i, vector<Vec2i>/*adjacent big-pixels of an edge*/, GRBVar>> > component_hedges;
	//vertical edges:
	vector< vector<tuple<Vec2i, vector<Vec2i>/*adjacent big-pixels of an edge*/, GRBVar>> > component_vedges;
	for (int i = 0; i < components.size(); i++)
	{
		//create horizontal edges:
		vector<tuple<Vec2i, vector<Vec2i>, GRBVar>> hedges;
		//(x0,y0) is the left or lower vertex position of a horizontal or vertical edge, respectively
		for (int vy = 0; vy <= new_height; vy++)
		{
			for (int vx = 0; vx < new_width; vx++)
			{
				vector<Vec2i> faces;  //adjacent face (X,Y) of this edge

				//face below?
				if (vy > 0)
				{
					int X = vx;
					int Y = vy - 1;

					//does the component have this face (X,Y)?
					if (components[i].count(Y * new_width + X) > 0)
					{
						faces.push_back(Vec2i(X, Y));
					}
				}
				//face above?
				if (vy < new_height)
				{
					int X = vx;
					int Y = vy;

					//does the component have this face (X,Y)?
					if (components[i].count(Y * new_width + X) > 0)
					{
						faces.push_back(Vec2i(X, Y));
					}
				}

				if (faces.size() > 0)
				{
					hedges.push_back(make_tuple(Vec2i(vx,vy), faces, model.addVar(0, 1, 0, GRB_BINARY)));
				}
			}
		}
		component_hedges.push_back(hedges);

		//vertical edges:
		vector<tuple<Vec2i, vector<Vec2i>, GRBVar>> vedges;
		//(vx,vy) is the lower vertex
		for (int vy = 0; vy < new_height; vy++)
		{
			for (int vx = 0; vx <= new_width; vx++)
			{
				vector<Vec2i> faces;  //adjacent face (X,Y) of this edge

				//face left?
				if (vx > 0)
				{
					int X = vx - 1;
					int Y = vy;

					//does the component have this face (X,Y)?
					if (components[i].count(Y * new_width + X) > 0)
					{
						faces.push_back(Vec2i(X, Y));
					}
				}
				//face right?
				if (vx < new_width)
				{
					int X = vx;
					int Y = vy;

					//does the component have this face (X,Y)?
					if (components[i].count(Y * new_width + X) > 0)
					{
						faces.push_back(Vec2i(X, Y));
					}
				}

				if (faces.size() > 0)
				{
					vedges.push_back(make_tuple(Vec2i(vx,vy), faces, model.addVar(0, 1, 0, GRB_BINARY)));
				}
			}
		}
		component_vedges.push_back(vedges);
	}

	model.update();

	//objective function: maximize scores of active big pixels
	GRBLinExpr obj;
	for (int i = 0; i < components.size(); i++)
	{
		unordered_map<int, BigPixel>& component = components[i];
		for (unordered_map<int, BigPixel>::iterator itr = component.begin(); itr != component.end(); itr++)
		{
			obj += (*itr).second.score * (*itr).second.var;
		}
	}
	model.setObjective(obj, GRB_MAXIMIZE);

	//for every big-pixel, it is occupied by exactly one big-pixel var (of different components)
	for (int Y = 0; Y < new_width; Y++)
	{
		for (int X = 0; X < new_width; X++)
		{
			GRBLinExpr sum;

			unordered_map<int,bool>& components_here = bigpixels_map[Y * new_width + X];
			for (unordered_map<int, bool>::iterator itr = components_here.begin(); itr != components_here.end(); itr++)
			{
				sum += components[(*itr).first][Y * new_width + X].var;
			}

			model.addConstr(sum == 1);
		}
	}

	//for every connected component, at least one of its big-pixels needs to be active
	for (int i = 0; i < components.size(); i++)
	{			
		GRBLinExpr sum;

		unordered_map<int, BigPixel>& component = components[i];
		for (unordered_map<int, BigPixel>::iterator itr = component.begin(); itr != component.end(); itr++)
		{
			sum += (*itr).second.var;
		}
		model.addConstr(sum >= 1);
	}

	//for water components, because they are 4-connectivity, diagonal faces are not allowed?
	if(true)
	{
		for (int c = 0; c < components.size(); c++)
		{
			if (c < num_lands)
				continue;

			unordered_map<int, BigPixel>& C = components[c];

			//check every interior vertex position:
			for (int vy = 1; vy < new_height; vy++)
			{
				for (int vx = 1; vx < new_width; vx++)
				{
					//this vertex (vx,vy)'s 4 adjacent face positions:
					Vec2i f0(vx - 1, vy - 1);
					Vec2i f1(vx, vy - 1);
					Vec2i f2(vx, vy);
					Vec2i f3(vx - 1, vy);

					//the presences of the 4 faces in this component:
					bool F0 = C.count(f0.y * new_width + f0.x) > 0;
					bool F1 = C.count(f1.y * new_width + f1.x) > 0;
					bool F2 = C.count(f2.y * new_width + f2.x) > 0;
					bool F3 = C.count(f3.y * new_width + f3.x) > 0;
					int count = F0 + F1 + F2 + F3;

					if (count <= 1)
						continue;  //no worry
					else if (count == 2)
					{
						//two possible diagonal ways:
						if (F0 && F2)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							model.addConstr(FF0 + FF2 <= 1);
						}
						else if (F1 && F3)
						{
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF1 + FF3 <= 1);
						}
					}
					else if (count == 3)
					{
						//4 possible ways:
						if (!F0)
						{
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF1 + (1 - FF2) + FF3 <= 2);
						}
						else if (!F1)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF0 + (1 - FF3) + FF2 <= 2);
						}
						else if (!F2)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF1 + (1 - FF0) + FF3 <= 2);
						}
						else if (!F3)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							model.addConstr(FF0 + (1 - FF1) + FF2 <= 2);
						}
					}
					else if (count == 4) //all present:
					{
						GRBVar FF0 = C[f0.y * new_width + f0.x].var;
						GRBVar FF1 = C[f1.y * new_width + f1.x].var;
						GRBVar FF2 = C[f2.y * new_width + f2.x].var;
						GRBVar FF3 = C[f3.y * new_width + f3.x].var;
						model.addConstr(FF0 + (1 - FF1) + FF2 + (1 - FF3) <= 3);
						model.addConstr((1 - FF0) + FF1 + (1 - FF2) + FF3 <= 3);
					}
				}
			}
		}
	}

	//for each land or water pixel, its cannot be adajcent to big-pixels of any other land or water pixels
	//let's check every interior half-edges, sum{ he->f + sum{all of he->f's incompatible faces) } <= 1
	//v-connect=true case: (for land only?) also check each of the two diagonal ways
	{
		//collect candidates at "adjacent" faces - either sharing an edge, or sharing a vertex (diagonally)
		vector<tuple<Vec2i, Vec2i,bool>> candidates_pairs;  //<face_pos,face_pos,diagonal-flag>

		//horizontal interior half-edges:
		for (int vy = 1; vy < new_height; vy++)
		{
			for (int vx = 0; vx < new_width; vx++)
			{
				//two ways:
				for (int Case = 0; Case < 2; Case++)
				{
					Vec2i F, FO;
					if (Case == 0)
					{
						F = Vec2i(vx, vy);
						FO = Vec2i(vx, vy - 1);
					}
					else
					{
						F = Vec2i(vx, vy - 1);
						FO = Vec2i(vx, vy);
					}

					candidates_pairs.push_back(make_tuple(F, FO, false));
				}
			}
		}

		//vertical interior half-edges
		for (int vy = 0; vy < new_height; vy++)
		{
			for (int vx = 1; vx < new_width; vx++)
			{
				//two ways:
				for (int Case = 0; Case < 2; Case++)
				{
					Vec2i F, FO;
					if (Case == 0)
					{
						F = Vec2i(vx, vy);
						FO = Vec2i(vx - 1, vy);
					}
					else
					{
						F = Vec2i(vx - 1, vy);
						FO = Vec2i(vx, vy);
					}

					candidates_pairs.push_back(make_tuple(F, FO, false/*not diagonal*/));
				}
			}
		}
		
		//for every interior vertex, check its two diagonal ways:
		for (int vy = 1; vy < new_height; vy++)
		{
			for (int vx = 1; vx < new_width; vx++)
			{
				//2 cases of two diagonal ways:
				for (int Case = 0; Case < 4; Case++)
				{
					Vec2i F, FO;
					if (Case == 0)
					{
						F = Vec2i(vx - 1, vy - 1);
						FO = Vec2i(vx, vy);
					}
					else if (Case == 1)
					{
						F = Vec2i(vx, vy);
						FO = Vec2i(vx - 1, vy - 1);
					}
					else if (Case == 2)
					{
						F = Vec2i(vx, vy - 1);
						FO = Vec2i(vx - 1, vy);
					}
					else if (Case == 3)
					{
						F = Vec2i(vx - 1, vy);
						FO = Vec2i(vx, vy - 1);
					}

					candidates_pairs.push_back(make_tuple(F, FO, true/*diagonal*/));
				}
			}
		}

		for (int i = 0; i < candidates_pairs.size(); i++)
		{
			Vec2i F = get<0>(candidates_pairs[i]);
			Vec2i FO = get<1>(candidates_pairs[i]);
			bool is_diagonal = get<2>(candidates_pairs[i]);

			unordered_map<int, bool>& F_cs = bigpixels_map[F.y * new_width + F.x];
			unordered_map<int, bool>& FO_cs = bigpixels_map[FO.y * new_width + FO.x];

			//for every big-pixel candidate at F:
			for (unordered_map<int, bool>::iterator itr = F_cs.begin(); itr != F_cs.end(); itr++)
			{
				int c0 = (*itr).first;
				bool c0_is_land = (c0 < num_lands);

				//collect incompatible big-pixel candidates at FO:
				GRBLinExpr sum;
				for (unordered_map<int, bool>::iterator itr = FO_cs.begin(); itr != FO_cs.end(); itr++)
				{
					int c1 = (*itr).first;
					bool c1_is_land = (c1 < num_lands);

					if (!is_diagonal)
					{
						//land cannot connect to any another land and
						//water cannot connect to any another water
						if ((c0_is_land && c1_is_land && c0 != c1) ||
							(!c0_is_land && !c1_is_land && c0 != c1))
						{
							sum += components[c1][FO.y * new_width + FO.x].var;
						}
					}
					else
					{
						//diagonal: only for land-land case
						if ((c0_is_land && c1_is_land && c0 != c1))
						{
							sum += components[c1][FO.y * new_width + FO.x].var;
						}
					}
				}

				GRBVar F_var = components[c0][F.y * new_width + F.x].var;
				model.addConstr(F_var + sum <= 1);
			}
		}
	}

	//vertex and edges constraints:
	for (int i = 0; i < components.size(); i++)
	{
		unordered_map<int, BigPixel>& component = components[i];

		for (int Case = 0; Case < 3; Case++)  //vertices, hedges, vedges
		{
			vector<tuple<Vec2i, vector<Vec2i>, GRBVar>>* elements;
			if (Case == 0)
				elements = &component_vertices[i];
			else if(Case == 1)
				elements = &component_hedges[i];
			else if(Case == 2)
				elements = &component_vedges[i];

			//check every pos of elements:
			for (int pos = 0; pos < elements->size(); pos++)
			{
				vector<Vec2i>& adj_big_pixels = std::get<1>( (*elements)[pos] );
				int N = adj_big_pixels.size();
				GRBLinExpr sum;
				for (int j = 0; j < adj_big_pixels.size(); j++)
				{
					Vec2i XY = adj_big_pixels[j];
					sum += component[XY.y * new_width + XY.x].var;
				}

				GRBVar var = std::get<2>((*elements)[pos]);
				
				
				model.addConstr(-N + 1 <= sum - N * var);
				model.addConstr(sum - N * var <= 0);
			}
		}
	}

	//Euler characteristics constraints (per-components):
	for (int i = 0; i < components.size(); i++)
	{
		//test: skip this for "oceans" (waters touching boundaries)?
		if (i >= num_lands && is_boundary_flags[i])
		{
			//cout << "skip Euler constraint for component#" << i << endl;
			continue;
		}

		GRBLinExpr sum_F;
		for (unordered_map<int, BigPixel>::iterator itr= components[i].begin(); itr != components[i].end(); itr++)
		{
			sum_F += (*itr).second.var;
		}

		GRBLinExpr sum_V;
		for (int j = 0; j < component_vertices[i].size(); j++)
		{
			sum_V += std::get<2>(component_vertices[i][j]);
		}

		GRBLinExpr sum_E;
		for (int j = 0; j < component_hedges[i].size(); j++)
		{
			sum_E += std::get<2>(component_hedges[i][j]);
		}
		for (int j = 0; j < component_vedges[i].size(); j++)
		{
			sum_E += std::get<2>(component_vedges[i][j]);
		}

		//shall equal to the component's expected Euler characteristic
		int Euler = get<0>(Eulers[i]) - get<1>(Eulers[i]) + get<2>(Eulers[i]);
		model.addConstr(sum_V - sum_E + sum_F == Euler);
	}

	if (g_ds_print_debug)
	{
		cout << "[Downsampling] program setup time:" << timeGetTime() - time << endl;
		time = timeGetTime();
	}

	//solve!
	model.getEnv().set(GRB_IntParam_OutputFlag, false);  //silent
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible
		printf("[Downsampling] the problem is infeasible.");
		return 2;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		printf("[Downsampling] optimize failed! status:%d", status);
		return 1;
	}

	float runtime = model.get(GRB_DoubleAttr_Runtime);

	printf("[Downsampling] optimize done! time:%f", runtime);

	////get results

	memset(output, false, sizeof(bool)* new_width * new_height);
	for (int i = 0; i < components.size(); i++)
	{
		for (unordered_map<int, BigPixel>::iterator itr = components[i].begin(); itr != components[i].end(); itr++)
		{
			Vec2i XY = (*itr).second.pos;

			bool flag = (bool)((*itr).second.var.get(GRB_DoubleAttr_X));
			if (flag)
			{
				if (i < num_lands)
					output[XY.y * new_width + XY.x] = true;
				else
					output[XY.y * new_width + XY.x] = false;
			}
		}
	}
		
	if (g_ds_print_debug)
	{
		/*cout << "output:" << endl;
		for (int y = 0; y < new_height; y++)
		{
			for (int x = 0; x < new_width; x++)
			{
				cout << output[y * new_width + x] << " ";
			}
			cout << endl;
		}*/

		/*cout << "solved_components:" << endl;
		for (int y = 0; y < new_height; y++)
		{
			for (int x = 0; x < new_width; x++)
			{
				cout << solved_components[y * new_width + x] << " ";
			}
			cout << endl;
		}*/

		for (int i = 0; i < components.size(); i++)
		{
			int sum_F = 0;

			for (unordered_map<int, BigPixel>::iterator itr = components[i].begin(); itr != components[i].end(); itr++)
			{
				bool flag = (*itr).second.var.get(GRB_DoubleAttr_X);

				if (flag)
					sum_F++;
			}

			int sum_V = 0;			
			vector<tuple< Vec2i, vector<Vec2i>, GRBVar>>& vertices = component_vertices[i];
			for (int j = 0; j < vertices.size(); j++)
			{
				bool flag = std::get<2>(vertices[j]).get(GRB_DoubleAttr_X);

				if (flag)
					sum_V++;
			}

			int sum_E = 0;

			vector<tuple< Vec2i, vector<Vec2i>, GRBVar>>& hedges = component_hedges[i];
			for (int j = 0; j < hedges.size(); j++)
			{
				bool flag = std::get<2>(hedges[j]).get(GRB_DoubleAttr_X);

				if (flag)
					sum_E++;
			}
			vector<tuple< Vec2i, vector<Vec2i>, GRBVar>>& vedges = component_vedges[i];
			for (int j = 0; j < vedges.size(); j++)
			{
				bool flag = std::get<2>(vedges[j]).get(GRB_DoubleAttr_X);

				if (flag)
					sum_E++;
			}
				
			int Euler = sum_V - sum_E + sum_F;
			cout << "C#" << i << "Euler:" << Euler << " V:" << sum_V << " E:" << sum_E << " sum_F:" << sum_F << endl;
		}
	}

	cout << "[Downsample] done. total time:" << timeGetTime() - time_begin << endl;
	return 0;
}


int DSSpace::Downsample(int width, int height, bool* mask, int bigpixel_width, int bigpixel_height, bool* output)
{
	const float BIG_NUM_MULTIPLIER = 1.5;

	if (width % bigpixel_width != 0 || height % bigpixel_width != 0)
	{
		cout << "[Downsampling] error: width / height not dividable" << endl;
		return 1;
	}

	DWORD time_begin = timeGetTime();

	GRBEnv env;

	int cur_width = width;
	int cur_height = height;
	bool* cur_mask = new bool[cur_width * cur_height];
	memcpy(cur_mask, mask, sizeof(bool) * cur_width * cur_height);

	//let's label the topology of current mask first
	unordered_map<int, pair<Vec2i, pair<bool, int>>> labels;
	int num_lands = 0;
	int num_waters = 0;

	DWORD time = time_begin;

	LabelTopology2(cur_width, cur_height, cur_mask, labels, num_lands, num_waters);

	//find land-water boundaries
	vector<pair<Vec2i,int>> boundaries;  //<land_index-water_index, size>
	EnumerateBoundaries(cur_width, cur_height, labels, num_lands, num_waters, boundaries);

	if (g_ds_print_debug)
	{
		cout << "[Downsampling] LabelTopology/EnumerateBoundaries time:" << timeGetTime() - time << endl;
		time = timeGetTime();
	}

	//size of the new mask:
	int new_width = cur_width / bigpixel_width;
	int new_height = cur_height / bigpixel_height;

	GRBModel model(env);

	//for every connected compoent (islands and waters), collect its overlapping big-pixels
	//and create the Boolean vars for each of its big-pixels
	struct BigPixel  //a big-pixel of a component
	{
		Vec2i pos;  //in output grid
		int score;  //"score" of this big-pixel w.s.t. the component
		GRBVar var;

		BigPixel()
		{
			pos = Vec2i(0);
			score = 0;
		}
	};
	vector<unordered_map<int, BigPixel>> components;  //big-pixels of every component (lands first then waters)
	components.resize(num_lands + num_waters);

	//also record the occupying potential big-pixels of each big-pixel locations 
	vector<unordered_map<int, bool>> bigpixels_map;  //row-major order. each is a map of component#
	bigpixels_map.resize(new_width * new_height);

	//"pixel-neighborhood" way:
	{
		const int land_weight = g_ds_land_weight;
		const int water_weight = 1;

		//get dist-based neighborhood mask
		vector<Vec2i> neighborhood;
		NeighborhoodMask(bigpixel_width / 4 + g_ds_neighobrhood_offset, neighborhood);

		for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
		{
			Vec2i xy = (*itr).second.first;

			//the actual component index:
			int c = -1;
			if (!(*itr).second.second.first)
				c = (*itr).second.second.second;
			else
				c = num_lands + (*itr).second.second.second;

			unordered_map<int, BigPixel>& component = components[c];

			int weight = 0;
			if (c < num_lands)  //land pixel?
				weight = land_weight;
			else
				weight = water_weight;

			//propagate scores within a neighborood, add to bigpixel's scores
			//assume the neighborhood is half size of a bigpixel?
			for (int ii = 0; ii < neighborhood.size(); ii++)
			{
				Vec2i xy_ = xy + neighborhood[ii];

				//within?
				if (xy_.x >= 0 && xy_.x < width && xy_.y >= 0 && xy_.y < height)
				{
					//to residing big-pixel coordinate:
					Vec2i XY(xy_.x / bigpixel_width, xy_.y / bigpixel_height);

					//first time for this big-pixel?
					if (component.count(XY.y * new_width + XY.x) == 0)
					{
						component[XY.y * new_width + XY.x].pos = XY;
						component[XY.y * new_width + XY.x].score = weight;
						component[XY.y * new_width + XY.x].var = model.addVar(0, 1, 0, GRB_BINARY);

						bigpixels_map[XY.y * new_width + XY.x][c] = true;
					}
					else
					{
						//accumulate a score
						component[XY.y * new_width + XY.x].score += weight;
					}
				}
			}
		}
	}

	////for each land-water boundary, we need to constraint its single-loop topology

	//(constant) adjacent land and water faces of each of the 12 vertex configurations
	vector< pair< vector<Vec2i>/*land poss, relative*/, vector<Vec2i>/*water poss*/ > /*12 total*/> VCFaces;
	for (int type = 0; type < 12; type++)
	{
		vector<Vec2i> land_offsets;  //values are -1 or 0
		vector<Vec2i> water_offsets;

		//incoming half-edge is +x:
		if (type == 0)
		{
			//ougoing half-edge is +y
			land_offsets.push_back(Vec2i(-1, 0));
			water_offsets.push_back(Vec2i(-1, -1));
			water_offsets.push_back(Vec2i(0, 0));
		}
		else if (type == 1)
		{
			//ougoing half-edge is +x
			land_offsets.push_back(Vec2i(-1, 0));
			land_offsets.push_back(Vec2i(0, 0));
			water_offsets.push_back(Vec2i(-1, -1));
			water_offsets.push_back(Vec2i(0, -1));
		}
		else if (type == 2)
		{
			//ougoing half-edge is -y
			land_offsets.push_back(Vec2i(-1, 0));
			land_offsets.push_back(Vec2i(0, -1));
			water_offsets.push_back(Vec2i(-1, -1));
		}

		//incoming he is +y
		else if (type == 3)
		{
			//ougoing half-edge is -x
			land_offsets.push_back(Vec2i(-1, -1));
			water_offsets.push_back(Vec2i(0, -1));
			water_offsets.push_back(Vec2i(-1, 0));
		}
		else if (type == 4)
		{
			//ougoing half-edge is +y
			land_offsets.push_back(Vec2i(-1, -1));
			land_offsets.push_back(Vec2i(-1, 0));
			water_offsets.push_back(Vec2i(0, -1));
			water_offsets.push_back(Vec2i(0, 0));
		}
		else if (type == 5)
		{
			//ougoing half-edge is +x
			land_offsets.push_back(Vec2i(-1, -1));
			land_offsets.push_back(Vec2i(0, 0));
			water_offsets.push_back(Vec2i(0, -1));
		}

		//incoming he is -x
		else if (type == 6)
		{
			//ougoing half-edge is -y
			land_offsets.push_back(Vec2i(0, -1));
			water_offsets.push_back(Vec2i(0, 0));
			water_offsets.push_back(Vec2i(-1, -1));
		}
		else if (type == 7)
		{
			//ougoing half-edge is -x
			land_offsets.push_back(Vec2i(0, -1));
			land_offsets.push_back(Vec2i(-1, -1));
			water_offsets.push_back(Vec2i(0, 0));
			water_offsets.push_back(Vec2i(-1, 0));
		}
		else if (type == 8)
		{
			//ougoing half-edge is +y
			land_offsets.push_back(Vec2i(0, -1));
			land_offsets.push_back(Vec2i(-1, 0));
			water_offsets.push_back(Vec2i(0, 0));
		}

		//incoming he is -y
		else if (type == 9)
		{
			//ougoing half-edge is +x
			land_offsets.push_back(Vec2i(0, 0));
			water_offsets.push_back(Vec2i(-1, 0));
			water_offsets.push_back(Vec2i(0, -1));
		}
		else if (type == 10)
		{
			//ougoing half-edge is -y
			land_offsets.push_back(Vec2i(0, 0));
			land_offsets.push_back(Vec2i(0, -1));
			water_offsets.push_back(Vec2i(-1, 0));
			water_offsets.push_back(Vec2i(-1, -1));
		}
		else if (type == 11)
		{
			//ougoing half-edge is -x
			land_offsets.push_back(Vec2i(0, 0));
			land_offsets.push_back(Vec2i(-1, -1));
			water_offsets.push_back(Vec2i(-1, 0));
		}

		VCFaces.push_back(make_pair(land_offsets, water_offsets));
	}

	//first collect candidate boundary vertices of every boundary
	vector< vector<pair<int/*0~11 type*/, Vec2i/*pos of vertex*/>> /*for every boundary*/> boundary_VCs(boundaries.size());
	for (int b = 0; b < boundaries.size(); b++)
	{
		unordered_map<int, BigPixel> &land = components[boundaries[b].first.x];
		unordered_map<int, BigPixel> &water = components[boundaries[b].first.y];

		//search all new-map interior vertices
		for (int vy = 1; vy < new_height; vy++)
		{
			for (int vx = 1; vx < new_width; vx++)
			{
				//there are 4 (incoming half-edge's dir) * 3 (outgoing half-edge's dir) = 12 possible vertex configurations
				//check each VC's eligibility

				//we assume he->f is land and he->o->f is water

				//for each VC, there are several land faces and water face candidates need to exist
				for (int type = 0; type < 12; type++)
				{
					vector<Vec2i> &land_offsets = VCFaces[type].first;  //values are -1 or 1
					vector<Vec2i> &water_offsets = VCFaces[type].second;

					//now, check if all the land and water faces bigpixels exist
					bool ok = true;
					for (int i = 0; i < land_offsets.size(); i++)
					{
						Vec2i p = Vec2i(vx, vy) + land_offsets[i];
						if (land.count(p.y * new_width + p.x) == 0)
						{
							ok = false;
							break;
						}
					}
					if (ok)
					{
						for (int i = 0; i < water_offsets.size(); i++)
						{
							Vec2i p = Vec2i(vx, vy) + water_offsets[i];
							if (water.count(p.y * new_width + p.x) == 0)
							{
								ok = false;
								break;
							}
						}
					}

					if (ok)
					{
						//gotcha! add a VC candidate to this boundary's VC list

						/*if (b == 2)
							cout << "new VC type:" << type << " pos:" << vx << "," << vy << endl;*/

						boundary_VCs[b].push_back(make_pair(type, Vec2i(vx, vy)));
					}

				}
			}
		}
	}

	//for every boundary:
	//1. create a Boolean var for every VC candidate
	//2. create a Integer "distance" var for every VC candidate
	//3. create a Boolean "last" var for every VC candidate
	vector< vector<GRBVar> > boundary_VC_vars;
	vector< vector<GRBVar> > boundary_dist_vars;
	vector< vector<GRBVar> > boundary_last_vars;
	for (int b = 0; b < boundary_VCs.size(); b++)
	{
		//the length of this boundary in original resolution
		int boundary_size_ori = boundaries[b].second;
		//so a good "BIG_NUM" for this boundary is the approximate length in new resolution multiply a scalar:
		//and cannot smaller than a value
		int BIG_NUM = MAX2((boundary_size_ori / bigpixel_width) * BIG_NUM_MULTIPLIER, 10);
		
		if(g_ds_print_debug)
			cout << "boundary_size_ori:" << boundary_size_ori << " BIG_NUM:" << BIG_NUM << endl;

		vector<GRBVar> VC_vars;
		vector<GRBVar> dist_vars;
		vector<GRBVar> last_vars;

		for (int i = 0; i < boundary_VCs[b].size(); i++)
		{
			VC_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
			dist_vars.push_back(model.addVar(0, BIG_NUM, 0, GRB_INTEGER));
			last_vars.push_back(model.addVar(0, 1, 0, GRB_BINARY));
		}

		boundary_VC_vars.push_back(VC_vars);
		boundary_dist_vars.push_back(dist_vars);
		boundary_last_vars.push_back(last_vars);
	}

	model.update();

	//objective function: maximize scores of active big pixels
	GRBLinExpr obj;
	for (int i = 0; i < components.size(); i++)
	{
		unordered_map<int, BigPixel>& component = components[i];
		for (unordered_map<int, BigPixel>::iterator itr = component.begin(); itr != component.end(); itr++)
		{
			obj += (*itr).second.score * (*itr).second.var;
		}
	}
	model.setObjective(obj, GRB_MAXIMIZE);

	//for every big-pixel, it is occupied by exactly one big-pixel var (of different components)
	for (int Y = 0; Y < new_height; Y++)
	{
		for (int X = 0; X < new_width; X++)
		{
			GRBLinExpr sum;

			unordered_map<int, bool>& components_here = bigpixels_map[Y * new_width + X];
			for (unordered_map<int, bool>::iterator itr = components_here.begin(); itr != components_here.end(); itr++)
			{
				sum += components[(*itr).first][Y * new_width + X].var;
			}

			model.addConstr(sum == 1);
		}
	}

	//for every connected component, at least one of its big-pixels needs to be active
	for (int i = 0; i < components.size(); i++)
	{
		GRBLinExpr sum;

		unordered_map<int, BigPixel>& component = components[i];
		for (unordered_map<int, BigPixel>::iterator itr = component.begin(); itr != component.end(); itr++)
		{
			sum += (*itr).second.var;
		}
		model.addConstr(sum >= 1);
	}

	//for water components, because they are 4-connectivity, diagonal faces are not allowed?
	if (true)
	{
		for (int c = 0; c < components.size(); c++)
		{
			if (c < num_lands)
				continue;

			unordered_map<int, BigPixel>& C = components[c];

			//check every interior vertex position:
			for (int vy = 1; vy < new_height; vy++)
			{
				for (int vx = 1; vx < new_width; vx++)
				{
					//this vertex (vx,vy)'s 4 adjacent face positions:
					Vec2i f0(vx - 1, vy - 1);
					Vec2i f1(vx, vy - 1);
					Vec2i f2(vx, vy);
					Vec2i f3(vx - 1, vy);

					//the presences of the 4 faces in this component:
					bool F0 = C.count(f0.y * new_width + f0.x) > 0;
					bool F1 = C.count(f1.y * new_width + f1.x) > 0;
					bool F2 = C.count(f2.y * new_width + f2.x) > 0;
					bool F3 = C.count(f3.y * new_width + f3.x) > 0;
					int count = F0 + F1 + F2 + F3;

					if (count <= 1)
						continue;  //no worry
					else if (count == 2)
					{
						//two possible diagonal ways:
						if (F0 && F2)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							model.addConstr(FF0 + FF2 <= 1);
						}
						else if (F1 && F3)
						{
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF1 + FF3 <= 1);
						}
					}
					else if (count == 3)
					{
						//4 possible ways:
						if (!F0)
						{
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF1 + (1 - FF2) + FF3 <= 2);
						}
						else if (!F1)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF0 + (1 - FF3) + FF2 <= 2);
						}
						else if (!F2)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF3 = C[f3.y * new_width + f3.x].var;
							model.addConstr(FF1 + (1 - FF0) + FF3 <= 2);
						}
						else if (!F3)
						{
							GRBVar FF0 = C[f0.y * new_width + f0.x].var;
							GRBVar FF1 = C[f1.y * new_width + f1.x].var;
							GRBVar FF2 = C[f2.y * new_width + f2.x].var;
							model.addConstr(FF0 + (1 - FF1) + FF2 <= 2);
						}
					}
					else if (count == 4) //all present:
					{
						GRBVar FF0 = C[f0.y * new_width + f0.x].var;
						GRBVar FF1 = C[f1.y * new_width + f1.x].var;
						GRBVar FF2 = C[f2.y * new_width + f2.x].var;
						GRBVar FF3 = C[f3.y * new_width + f3.x].var;
						model.addConstr(FF0 + (1 - FF1) + FF2 + (1 - FF3) <= 3);
						model.addConstr((1 - FF0) + FF1 + (1 - FF2) + FF3 <= 3);
					}
				}
			}
		}
	}

	//for each land or water pixel, its cannot be adajcent to big-pixels of any other land or water pixels
	//let's check every interior half-edges, sum{ he->f + sum{all of he->f's incompatible faces) } <= 1
	//v-connect=true case: (for land only?) also check each of the two diagonal ways
	if(g_ds_local_constraint)
	{
		//collect candidates at "adjacent" faces - either sharing an edge, or sharing a vertex (diagonally)
		vector<tuple<Vec2i, Vec2i, bool>> candidates_pairs;  //<face_pos,face_pos,diagonal-flag>

		//horizontal interior half-edges:
		for (int vy = 1; vy < new_height; vy++)
		{
			for (int vx = 0; vx < new_width; vx++)
			{
				//two ways:
				for (int Case = 0; Case < 2; Case++)
				{
					Vec2i F, FO;
					if (Case == 0)
					{
						F = Vec2i(vx, vy);
						FO = Vec2i(vx, vy - 1);
					}
					else
					{
						F = Vec2i(vx, vy - 1);
						FO = Vec2i(vx, vy);
					}

					candidates_pairs.push_back(make_tuple(F, FO, false));
				}
			}
		}

		//vertical interior half-edges
		for (int vy = 0; vy < new_height; vy++)
		{
			for (int vx = 1; vx < new_width; vx++)
			{
				//two ways:
				for (int Case = 0; Case < 2; Case++)
				{
					Vec2i F, FO;
					if (Case == 0)
					{
						F = Vec2i(vx, vy);
						FO = Vec2i(vx - 1, vy);
					}
					else
					{
						F = Vec2i(vx - 1, vy);
						FO = Vec2i(vx, vy);
					}

					candidates_pairs.push_back(make_tuple(F, FO, false/*not diagonal*/));
				}
			}
		}

		//for every interior vertex, check its two diagonal ways:
		for (int vy = 1; vy < new_height; vy++)
		{
			for (int vx = 1; vx < new_width; vx++)
			{
				//2 cases of two diagonal ways:
				for (int Case = 0; Case < 4; Case++)
				{
					Vec2i F, FO;
					if (Case == 0)
					{
						F = Vec2i(vx - 1, vy - 1);
						FO = Vec2i(vx, vy);
					}
					else if (Case == 1)
					{
						F = Vec2i(vx, vy);
						FO = Vec2i(vx - 1, vy - 1);
					}
					else if (Case == 2)
					{
						F = Vec2i(vx, vy - 1);
						FO = Vec2i(vx - 1, vy);
					}
					else if (Case == 3)
					{
						F = Vec2i(vx - 1, vy);
						FO = Vec2i(vx, vy - 1);
					}

					candidates_pairs.push_back(make_tuple(F, FO, true/*diagonal*/));
				}
			}
		}

		for (int i = 0; i < candidates_pairs.size(); i++)
		{
			Vec2i F = get<0>(candidates_pairs[i]);
			Vec2i FO = get<1>(candidates_pairs[i]);
			bool is_diagonal = get<2>(candidates_pairs[i]);

			unordered_map<int, bool>& F_cs = bigpixels_map[F.y * new_width + F.x];
			unordered_map<int, bool>& FO_cs = bigpixels_map[FO.y * new_width + FO.x];

			//for every big-pixel candidate at F:
			for (unordered_map<int, bool>::iterator itr = F_cs.begin(); itr != F_cs.end(); itr++)
			{
				int c0 = (*itr).first;
				bool c0_is_land = (c0 < num_lands);

				//collect incompatible big-pixel candidates at FO:
				GRBLinExpr sum;
				for (unordered_map<int, bool>::iterator itr = FO_cs.begin(); itr != FO_cs.end(); itr++)
				{
					int c1 = (*itr).first;
					bool c1_is_land = (c1 < num_lands);

					if (!is_diagonal)
					{
						//land cannot connect to any another land and
						//water cannot connect to any another water
						if ((c0_is_land && c1_is_land && c0 != c1) ||
							(!c0_is_land && !c1_is_land && c0 != c1))
						{
							sum += components[c1][FO.y * new_width + FO.x].var;
						}
					}
					else
					{
						//diagonal: only for land-land case
						if ((c0_is_land && c1_is_land && c0 != c1))
						{
							sum += components[c1][FO.y * new_width + FO.x].var;
						}
					}
				}

				GRBVar F_var = components[c0][F.y * new_width + F.x].var;
				model.addConstr(F_var + sum <= 1);
			}
		}
	}

	////model the boundary constraints:

	for (int b = 0; b < boundaries.size(); b++)
	{
		unordered_map<int, BigPixel>& land = components[boundaries[b].first.x];
		unordered_map<int, BigPixel>& water = components[boundaries[b].first.y];

		//the length of this boundary in original resolution
		int boundary_size_ori = boundaries[b].second;
		//so a good "BIG_NUM" for this boundary is the approximate length in new resolution multiply a scalar:
		//and cannot smaller than a value
		int BIG_NUM = MAX2((boundary_size_ori / bigpixel_width) * BIG_NUM_MULTIPLIER, 10);

		vector<pair<int/*0~11 type*/, Vec2i/*pos of vertex*/>> &VCs = boundary_VCs[b];
		vector<GRBVar>& VC_vars = boundary_VC_vars[b];
		vector<GRBVar>& dist_vars = boundary_dist_vars[b];
		vector<GRBVar>& last_vars = boundary_last_vars[b];

		//build a quick map of VCs from Vec2i to index
		map<int, vector<int>> VCs_map;  //key: y*new_width+x, value: list of indices in VCs at the position
		for (int i = 0; i < VCs.size(); i++)
		{
			VCs_map[VCs[i].second.y * new_width + VCs[i].second.x].push_back(i);
		}

		for (int i = 0; i < VCs.size(); i++)
		{
			int type = VCs[i].first;  //0~11
			Vec2i pos = VCs[i].second;
			
			//locate the pointing-to pos according to type:
			Vec2i pos_to;
			if (type == 0 || type == 4 || type == 8)  //out-going he is +y
			{
				pos_to = pos + Vec2i(0, 1);
			}
			else if (type == 1 || type == 5 || type == 9)  //+x
			{
				pos_to = pos + Vec2i(1, 0);
			}
			else if (type == 2 || type == 6 || type == 10)  //-y
			{
				pos_to = pos + Vec2i(0, -1);
			}
			else if (type == 3 || type == 7 || type == 11)  //-x
			{
				pos_to = pos + Vec2i(-1, 0);
			}

			GRBVar VC_var = VC_vars[i];
			GRBVar dist_var = dist_vars[i];

			//if this VC's pointing-to pos has no existing VC, this VC cannot be active
			if (VCs_map.count(pos_to.y * new_width + pos_to.x) == 0)
			{
				model.addConstr(VC_var == 0);
			}

			//if a VC is inactive, its dist var must be 0
			model.addConstr(dist_var <= BIG_NUM * VC_var);
			
			//connectivity-based constraints:
			{
				//for a VC to be active, the needed faces all shall match exactly
				vector<Vec2i>& land_offsets = VCFaces[VCs[i].first].first;  //value: relative 0 or 1
				vector<Vec2i>& water_offsets = VCFaces[VCs[i].first].second;

				GRBLinExpr sum;  //sum of the needed face vars
				for (int j = 0; j < land_offsets.size(); j++)
				{
					Vec2i p = pos + land_offsets[j];
					sum += land[p.y * new_width + p.x].var;
				}
				for (int j = 0; j < water_offsets.size(); j++)
				{
					Vec2i p = pos + water_offsets[j];
					sum += water[p.y * new_width + p.x].var;
				}
				int N = land_offsets.size() + water_offsets.size();
				model.addConstr(sum - N * VC_var >= 0);
				model.addConstr(sum - N * VC_var <= N-1);

				//for a VC from vertex#i to vertex#j:
				//VC_(i,j) - (sum<d_j>-d_i) ¡V LAST_i * BIG_NUM <= 0
				//d_i is the dist var of vertex#i
				//sum<d_j> is the sum of dist vars of compatible vertex#j
				//LAST_i is the "last" flag for vertex#i

				//find the VC at pos_to w/ compatible types (i.e., this VC's out dir equals next VC's in dir)
				vector<int> to_indices;
				vector<int>& candidates = VCs_map[pos_to.y * new_width + pos_to.x];
				
				int out_dir = 0;  //+y, +x, -y, -x
				if (type == 0 || type == 4 || type == 8)
					out_dir = 0;
				else if (type == 1 || type == 5 || type == 9)
					out_dir = 1;
				else if (type == 2 || type == 6 || type == 10)
					out_dir = 2;
				else if (type == 3 || type == 7 || type == 11)
					out_dir = 3;

				for (int k = 0; k < candidates.size(); k++)
				{
					int to_type = VCs[candidates[k]].first;

					int in_dir = 0;  //+y, +x, -y, -x
					if (to_type == 0 || to_type == 1 || to_type == 2)  //+x
						in_dir = 1;
					else if (to_type == 3 || to_type == 4 || to_type == 5)  //+y
						in_dir = 0;
					else if (to_type == 6 || to_type == 7 || to_type == 8)  //-x
						in_dir = 3;
					else if (to_type == 9 || to_type == 10 || to_type == 11)  //-y
						in_dir = 2;

					if (out_dir == in_dir)
					{
						//compatible!
						to_indices.push_back(candidates[k]);
					}
				}

				
				/*if (b == 2)
				{
					cout << "VC#" << i << " type:" << type << " pos:" << pos << " pos_to:" << pos_to << " tis: ";
					for (int k = 0; k < to_indices.size(); k++)
					{
						cout << to_indices[k] << ",";
					}
					cout << endl;
				}*/

				GRBLinExpr sum_dist_to_vars;
				for (int k = 0; k < to_indices.size(); k++)
				{
					sum_dist_to_vars += dist_vars[to_indices[k]];
				}
				
				GRBVar last_var = last_vars[i];

				model.addConstr(VC_var - (sum_dist_to_vars - dist_var) - last_var * BIG_NUM <= 0);
			}
		}

		//exactly one last flag for this boundary is true
		{
			GRBLinExpr sum;
			for (int i = 0; i < last_vars.size(); i++)
			{
				sum += last_vars[i];
			}
			model.addConstr(sum == 1);
		}
	}

	if (g_ds_print_debug)
	{
		cout << "[Downsampling] program setup time:" << timeGetTime() - time << endl;
		time = timeGetTime();
	}

	//solve!
	model.getEnv().set(GRB_DoubleParam_TimeLimit, 60);
	if(!g_ds_print_debug)
		model.getEnv().set(GRB_IntParam_OutputFlag, false);  //silent
	model.optimize();
	int status = model.get(GRB_IntAttr_Status);
	if (status == 3)
	{
		//infeasible or timeout
		printf("[Downsampling] the problem is infeasible");
		return 2;
	}
	else if (status != 9/*time-out*/ && status != 2 && status != 11 && status != 13)
	{
		//some other failure
		printf("[Downsampling] optimize failed! status:%d", status);
		return 1;
	}

	float runtime = model.get(GRB_DoubleAttr_Runtime);

	printf("[Downsampling] optimize done! time:%f", runtime);

	//get results

	if (g_ds_save_components_to_file)
	{
		g_ds_input_components.clear();
		g_ds_input_components.resize(width* height, 999999);

		g_ds_output_components.clear();
		g_ds_output_components.resize(new_width * new_height, 999999);
	}

	memset(output, false, sizeof(bool) * new_width * new_height);
	for (int i = 0; i < components.size(); i++)
	{
		for (unordered_map<int, BigPixel>::iterator itr = components[i].begin(); itr != components[i].end(); itr++)
		{
			Vec2i XY = (*itr).second.pos;

			bool flag = (bool)((*itr).second.var.get(GRB_DoubleAttr_X));
			if (flag)
			{
				if (i < num_lands)
					output[XY.y * new_width + XY.x] = true;
				else
					output[XY.y * new_width + XY.x] = false;

				//save output component indices?
				if (g_ds_save_components_to_file)
					g_ds_output_components[XY.y * new_width + XY.x] = i;
			}
		}
	}

	//save input component indices?
	if (g_ds_save_components_to_file)
	{
		for (unordered_map<int, pair<Vec2i, pair<bool, int>>>::iterator itr = labels.begin(); itr != labels.end(); itr++)
		{
			Vec2i xy = (*itr).second.first;

			int c = -1;
			if (!(*itr).second.second.first)
				c = (*itr).second.second.second;
			else
				c = num_lands + (*itr).second.second.second;

			g_ds_input_components[xy.y * width + xy.x] = c;
		}
	}

	if (g_ds_print_debug)
	{
		//save solved boundaries
		g_ds_output_boundaries.clear(); 
		for (int b = 0; b < boundaries.size(); b++)
		{
			vector<pair<int/*0~11 type*/, Vec2i/*pos of vertex*/>>& VCs = boundary_VCs[b];
			vector<GRBVar>& VC_vars = boundary_VC_vars[b];
			vector<GRBVar>& dist_vars = boundary_dist_vars[b];
			vector<GRBVar>& last_vars = boundary_last_vars[b];

			vector<tuple<Vec2i, int/*type*/, int/*dist*/>> output_boundaries;

			for (int i = 0; i < VC_vars.size(); i++)
			{
				bool flag = (bool)(VC_vars[i].get(GRB_DoubleAttr_X));
				if (flag)
				{
					int dist = (int)(dist_vars[i].get(GRB_DoubleAttr_X));

					bool last = (bool)(last_vars[i].get(GRB_DoubleAttr_X));

					if(b==2)
						cout << "active VC! type:" << VCs[i].first << " " << VCs[i].second << " dist:" << dist << " last:" << last << endl;

					output_boundaries.push_back(make_tuple(VCs[i].second, VCs[i].first, dist));
				}
			}

			cout << "boundary#" << b << ":" << output_boundaries.size() << endl;

			g_ds_output_boundaries.push_back(output_boundaries);
		}
	}

	cout << "[Downsample] done. total time:" << timeGetTime() - time_begin << endl;
	return 0;
}


bool DSSpace::ErrorMetricsPng(const char* input_filename, const char* output_filename, 
	float& IoU, float& Dice, float& Precision, float& Recall)
{
	std::vector<unsigned char> in_buffer; //the raw pixels (RGBA)
	unsigned width = 0, height = 0;
	unsigned error = lodepng::decode(in_buffer, width, height, input_filename);
	if (error)
	{
		cout << input_filename << " lodepng::decode error:" << error << " " << lodepng_error_text(error) << endl;
		return false;
	}

	std::vector<unsigned char> out_buffer; //the raw pixels (RGBA)
	unsigned new_width = 0, new_height = 0;
	error = lodepng::decode(out_buffer, new_width, new_height, output_filename);
	if (error)
	{
		cout << output_filename << " lodepng::decode error:" << error << " " << lodepng_error_text(error) << endl;
		return false;
	}

	//turn the image buffers to binary masks
	bool* input = new bool[width * height];
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int offset = (y * width + x) * 4;
			Vec4<unsigned char> p(in_buffer[offset], in_buffer[offset + 1], in_buffer[offset + 2], in_buffer[offset + 3]);

			//binarilization: <128 = black, >128 = white
			if (p[0] >= 128 && p[1] >= 128 && p[2] >= 128)
			{
				input[y * width + x] = true;
			}
			else
			{
				input[y * width + x] = false;
			}
		}
	}

	bool* output = new bool[new_width * new_height];
	for (int y = 0; y < new_height; y++)
	{
		for (int x = 0; x < new_width; x++)
		{
			int offset = (y * new_width + x) * 4;
			Vec4<unsigned char> p(out_buffer[offset], out_buffer[offset + 1], out_buffer[offset + 2], out_buffer[offset + 3]);

			//binarilization: <128 = black, >128 = white
			if (p[0] >= 128 && p[1] >= 128 && p[2] >= 128)
			{
				output[y * new_width + x] = true;
			}
			else
			{
				output[y * new_width + x] = false;
			}
		}
	}

	ErrorMetrics(width, height, input, new_width, new_height, output, IoU, Dice, Precision, Recall);

	delete input;
	delete output;
	return true;
}

bool DSSpace::ErrorMetrics(int width, int height, bool* input/*size = width*height */,
	int new_width, int new_height, bool* output/*size = new_width*new_height */, 
	float& IoU, float& Dice, float& Precision, float& Recall)
{
	if (!input || !output || (width % new_width) != 0 || (height % new_height) != 0 || width == 0 || height == 0)
	{
		cout << "[CalculateErrors] wrong inputs" << endl;
		return false;
	}

	//count # of small pixels for the following cases: 
	//In=1 && Out=1  (true positive)
	//In=0 && Out=1  (false positive)
	//In=1 && Out=0  (false negative)
	//In=0 && Out=0  (true negative)
	int num_TP = 0, num_FP = 0, num_FN = 0, num_TN = 0;

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			bool In = input[y * width + x];
			bool Out = output[(int)(y * (float)new_height / (float)height ) * new_width +
				(int)(x * (float)new_width / (float)width)];

			if (In && Out)
				num_TP++;
			else if (!In && Out)
				num_FP++;
			else if (In && !Out)
				num_FN++;
			else if (!In && !Out)
				num_TN++;
		}
	}

	IoU = (float)num_TP / ((float)num_TP + (float)num_FP + (float)num_FN);

	//Dice = F1 score
	Dice = 2 * (float)num_TP / (2 * (float)num_TP + (float)num_FP + (float)num_FN);

	Precision = (float)num_TP / ((float)num_TP + (float)num_FP);

	Recall = (float)num_TP / (float)(num_TP + num_FN);

	printf("[ErrorMetrics] IoU:%f Dice:%f Precision:%f Recall:%f", IoU, Dice, Precision, Recall);
	return true;
}

//utility function for counting # of connected components of val in a N8-mask or a N4-mask
//buffer: 3x3 neighborhood in left-to-right lower-to-upper row-major buffer, starts at the lower-left corner
int ConnectedComponentsN8(vector<bool> &buffer)
{
	bool val = buffer[4];

	//to CCW order starts from lower-left corner:
	bool f0 = !(buffer[0] != val);
	bool f1 = !(buffer[1] != val);
	bool f2 = !(buffer[2] != val);
	bool f3 = !(buffer[5] != val);
	bool f4 = !(buffer[8] != val);
	bool f5 = !(buffer[7] != val);
	bool f6 = !(buffer[6] != val);
	bool f7 = !(buffer[3] != val);

	int sum = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7;

	if (sum == 0)
		return 5;  //zero components
	else if (sum == 1)
		return 1;
	else if (sum == 2)
	{
		//one component?
		if ((f0 && f1) || (f1 && f2) || (f2 && f3) || (f3 && f4) || (f4 && f5) || (f5 && f6) ||
			(f6 && f7) || (f7 && f0))
		{
			return 1;
		}
		//otherwise, two components
		else
			return 2;
	}
	else if (sum == 3)
	{
		//one 3-component?
		if ((f0 && f1 && f2) || (f1 && f2 && f3) || (f2 && f3 && f4) || (f3 && f4 && f5) || (f4 && f5 && f6) ||
			(f5 && f6 && f7) ||	(f6 && f7 && f0) || (f7 && f0 && f1))
		{
			return 1;
		}
		//two components?
		//(one 2-component, another single component)
		else if ((f0 && f1) || (f1 && f2) || (f2 && f3) || (f3 && f4) || (f4 && f5) || (f5 && f6) ||
			(f6 && f7) || (f7 && f0))
		{
			return 2;
		}
		//otherwise 3 1-components
		else
		{
			return 3;
		}
	}
	else if (sum == 4)
	{
		//one 4-component?
		if ((f0 && f1 && f2 && f3) || (f1 && f2 && f3 && f4) || (f2 && f3 && f4 && f5) || (f3 && f4 && f5 && f6) ||
			(f4 && f5 && f6 && f7) || (f5 && f6 && f7 && f0 ) || (f6 && f7 && f0 && f1) || (f7 && f0 && f1 && f2))
		{
			return 1;
		}
		//one 3-component and one 1-component?
		else if ((f0 && f1 && f2) || (f1 && f2 && f3) || (f2 && f3 && f4) || (f3 && f4 && f5) || (f4 && f5 && f6) ||
			(f5 && f6 && f7) || (f6 && f7 && f0) || (f7 && f0 && f1))
		{
			return 2;
		}
		//one 2-component and one 2-component?
		else if ((f0 && f1 && f3 && f4) || (f0 && f1 && f4 && f5) || (f0 && f1 && f5 && f6) ||
			(f1 && f2 && f4 && f5) || (f1 && f2 && f5 && f6) || (f1 && f2 && f6 && f7) || 
			(f2 && f3 && f5 && f6) || (f2 && f3 && f6 && f7) || (f2 && f3 && f7 && f0) || 
			(f3 && f4 && f6 && f7) || (f3 && f4 && f7 && f0) || 
			(f4 && f5 && f7 && f0) )
		{
			return 2;
		}
		//4 1-components? (interleaved, two possible cases)
		else if( (f0 && !f1 && f2 && !f3 && f4 && !f5 && f6 && !f7) ||
			(!f0 && f1 && !f2 && f3 && !f4 && f5 && !f6 && f7) )
		{
			return 4;
		}
		//otherwise, must be: one 2-component, one 1-component, one 1-component
		else
		{
			return 3;
		}
	}
	else if (sum == 5)
	{
		//check false flags instead
		bool F0 = !f0;
		bool F1 = !f1;
		bool F2 = !f2;
		bool F3 = !f3;
		bool F4 = !f4;
		bool F5 = !f5;
		bool F6 = !f6;
		bool F7 = !f7;

		//one 3-component?
		if ((F0 && F1 && F2) || (F1 && F2 && F3) || (F2 && F3 && F4) || (F3 && F4 && F5) || (F4 && F5 && F6) ||
			(F5 && F6 && F7) || (F6 && F7 && F0) || (F7 && F0 && F1))
		{
			return 1;
		}
		//two components?
		//(one 2-component, another single component)
		else if ((F0 && F1) || (F1 && F2) || (F2 && F3) || (F3 && F4) || (F4 && F5) || (F5 && F6) ||
			(F6 && F7) || (F7 && F0))
		{
			return 2;
		}
		//otherwise 3 1-components
		else
		{
			return 3;
		}
	}
	else if (sum == 6)
	{
		//check false flags instead
		bool F0 = !f0;
		bool F1 = !f1;
		bool F2 = !f2;
		bool F3 = !f3;
		bool F4 = !f4;
		bool F5 = !f5;
		bool F6 = !f6;
		bool F7 = !f7;

		//one component?
		if ((F0 && F1) || (F1 && F2) || (F2 && F3) || (F3 && F4) || (F4 && F5) || (F5 && F6) ||
			(F6 && F7) || (F7 && F0))
		{
			return 1;
		}
		//otherwise, two components
		else
			return 2;
	}
	else if (sum == 7)
		return 1;
	else if (sum == 8)
		return 1;

	cout << "[ConnectedComponentsN8] what?" << endl;
	return -1;
}
int ConnectedComponentsN4(vector<bool> &buffer)
{
	bool val = buffer[4];

	//to CCW order:
	bool f0 = !(buffer[1] != val);
	bool f1 = !(buffer[5] != val);
	bool f2 = !(buffer[7] != val);
	bool f3 = !(buffer[3] != val);
	
	int sum = f0 + f1 + f2 + f3;

	if (sum == 0)
		return 5;
	else if (sum == 1)
		return 1;
	else if (sum == 2)
	{
		if ((f0 && f1) || (f1 && f2) || (f2 && f3) || (f3 && f0))
			return 1;
		else
			return 2;
	}
	else //sum == 3 or == 4
		return 1;

	return 0;
}

bool DSSpace::DownsampleACN(int width, int height, bool* mask/*size = width*height */, bool* output)
{
	if (width % 2 != 0 || height % 2 != 0)
	{
		cout << "[DownsampleACN] error: width / height not dividable" << endl;
		return false;
	}

	//calculate per-pixel adaptive crossing numbers
	vector<int> ACNs(width * height, 0);

	//note: just assume boundary pixels' ACNs are 0. we don't tackle them
	for (int Y = 1; Y < height - 1; Y++)
	{
		for (int X = 1; X < width - 1; X++)
		{
			//prepare the 3x3 neighborhood:
			vector<bool> neighborhood;  //row-major left-to-right bottom-to-top
			for (int y_ = -1; y_ <= 1; y_++)
			{
				for (int x_ = -1; x_ <= 1; x_++)
				{
					neighborhood.push_back(mask[(Y + y_)*width + (X + x_)]);
				}
			}

			bool val = neighborhood[4];

			//p is this point (x,y)
			// 
			//calculate n^I(p)
			int nIp = 0;
			for (int y_ = 0; y_ <= 2; y_++)
			{
				for (int x_ = 0; x_ <= 2; x_++)
				{
					if (x_ == 1 && y_ == 1)
						continue;  //skip itself
					
					if (neighborhood[y_ * 3 + x_] == val)
						nIp++;
				}
			}

			//depends on nI(p):
			if (nIp < 4)  //count # of N4-connected components of p's N8-neighbors w/ opposite value
			{
				ACNs[Y * width + X] = ConnectedComponentsN8(neighborhood);
			}
			else  //count # of N8-connected components of p's N4-neighbors
			{
				ACNs[Y * width + X] = ConnectedComponentsN4(neighborhood);
			}
		}
	}

	//now, write to output buffer according to ACNs
	//(for every 2x2 bigpixel, pick the pixel with the highest ACN. pick first one in tie)
	int new_width = width / 2;
	int new_height = height / 2;
	for (int Y = 0; Y < new_height; Y++)
	{
		for (int X = 0; X < new_width; X++)
		{
			//(X,Y) is bigpixel coord

			//the ACNs of the 4 pixels:
			int ACN0 = ACNs[(Y * 2) * width + (X * 2)];
			int ACN1 = ACNs[(Y * 2) * width + (X * 2 + 1)];
			int ACN2 = ACNs[(Y * 2 + 1) * width + (X * 2 + 1)];
			int ACN3 = ACNs[(Y * 2 + 1) * width + (X * 2)];

			bool val = false;

			//ACN0 win?
			if (ACN0 >= ACN1 && ACN0 >= ACN2 && ACN0 >= ACN3)
			{
				val = mask[(Y * 2) * width + (X * 2)];
			}
			else if (ACN1 > ACN0 && ACN1 >= ACN2 && ACN1 >= ACN3)
			{
				val = mask[(Y * 2) * width + (X * 2 + 1)];
			}
			else if (ACN2 > ACN0 && ACN2 > ACN1 && ACN2 >= ACN3)
			{
				val = mask[(Y * 2 + 1) * width + (X * 2 + 1)];
			}
			else
			{
				val = mask[(Y * 2 + 1) * width + (X * 2)];
			}

			output[Y * new_width + X] = val;
		}
	}

	return true;
}

bool DSSpace::DownsampleACNPng(const char* input_filename)
{
	std::vector<unsigned char> in_buffer; //the raw pixels (RGBA)
	unsigned width = 0, height = 0;
	unsigned error = lodepng::decode(in_buffer, width, height, input_filename);
	if (error)
	{
		cout << input_filename << " lodepng::decode error:" << error << " " << lodepng_error_text(error) << endl;
		return false;
	}

	if ((width % 2) != 0 || (height % 2) != 0)
	{
		cout << "error: width / height not dividable!" << endl;
		return false;
	}

	//turn the image buffer to a binary mask buffer
	bool* mask = new bool[width * height];
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int offset = (y * width + x) * 4;
			Vec4<unsigned char> p(in_buffer[offset], in_buffer[offset + 1], in_buffer[offset + 2], in_buffer[offset + 3]);

			//binarilization: <threshold = black, >threshold = white
			if (p[0] >= g_ds_png_treshold && p[1] >= g_ds_png_treshold && p[2] >= g_ds_png_treshold)
			{
				mask[y * width + x] = true;
			}
			else
			{
				mask[y * width + x] = false;
			}
		}
	}

	int new_width = width / 2;
	int new_height = height / 2;
	bool* output = new bool[new_width * new_height];

	bool ret = DownsampleACN(width, height, mask, output);
	if(ret)
	{
		//save to a png file!
		std::vector<unsigned char> out_buffer(new_width * new_height * 4);
		for (int y = 0; y < new_height; y++)
		{
			for (int x = 0; x < new_width; x++)
			{
				int offset = (y * new_width + x) * 4;

				if (output[y * new_width + x])
				{
					//white
					out_buffer[offset] = 255;
					out_buffer[offset + 1] = 255;
					out_buffer[offset + 2] = 255;
					out_buffer[offset + 3] = 255;
				}
				else
				{
					//black
					out_buffer[offset] = 0;
					out_buffer[offset + 1] = 0;
					out_buffer[offset + 2] = 0;
					out_buffer[offset + 3] = 255;
				}
			}
		}

		string output_filename = string(input_filename) + "." + to_string(new_width) + "x" + to_string(new_height);
		output_filename += ".ACN.png";

		lodepng::encode(output_filename, out_buffer, new_width, new_height);
	}
	
	return ret;
}

int DSSpace::DownsamplePassat2022(int width, int height, bool* mask, int bigpixel_size, 
	bool* output, bool print_debug)
{

	//dimension of the result grid
	const int Width = width / bigpixel_size;
	const int Height = height / bigpixel_size;

	const int iter_max = width * height / 4;

	//for every big-pixel, calculate its "Sigma" (=1 if >50% are black, =-1 if <=50% are black)
	bool *Sigmas = new bool[Width * Height];  //true=1, false=-1
	int threshold = bigpixel_size * bigpixel_size / 2;
	for (int Y = 0; Y < Height; Y++)
	{
		for (int X = 0; X < Width; X++)
		{
			int count = 0;  //# of blacks small pixels?
			for (int y = Y * bigpixel_size; y < Y * bigpixel_size + bigpixel_size; y++)
			{
				for (int x = X * bigpixel_size; x < X * bigpixel_size + bigpixel_size; x++)
				{
					if (mask[y * width + x])
						count++;
				}
			}
			if (count > threshold)
				Sigmas[Y * Width + X] = true;
			else
				Sigmas[Y * Width + X] = false;
		}
	}

	//prepare the result "H" grid of the original dimension
	bool* H = new bool[width * height];
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			//initial with current grid
			H[y * width + x] = mask[y * width + x];
		}
	}

	//collect single-pixel "flip" candidates 
	//first: (x,y), second: true=black-to-white false=white-to-black, third: score (i.e., energy)
	typedef tuple<Vec2i, bool, float> Candidate;
	list<Candidate> candidates;

	//to calculate "zeta" values for every bigpixel, we need current # of black/white pixels 
	//within every current bigpixel of H
	int* BlackCounts = new int[Width * Height];
	for (int Y = 0; Y < Height; Y++)
	{
		for (int X = 0; X < Width; X++)
		{
			int count = 0;  //# of blacks small pixels?
			for (int y = Y * bigpixel_size; y < Y * bigpixel_size + bigpixel_size; y++)
			{
				for (int x = X * bigpixel_size; x < X * bigpixel_size + bigpixel_size; x++)
				{
					if (H[y * width + x])
						count++;
				}
			}

			BlackCounts[Y * Width + X] = count;
		}
	}

	class Worker
	{
	public:

		void TryAddCandidate(bool* H, int x, int y, int width, int height, 
			int Width, int Height, int bigpixel_size, int* BlackCounts, bool *Sigmas, 
			list<Candidate>& candidates)
		{
			//this pixel is a candidate if its 8-neighbors has exactly
			//1 black component and 1 white component

			//cout << x << "," << y << ": ";
			vector<bool> neighborhood;  //9 elements!
			for (int yy = y - 1; yy <= y + 1; yy++)
			{
				for (int xx = x - 1; xx <= x + 1; xx++)
				{
					//cout << H[yy * width + xx] << " ";
					neighborhood.push_back(H[yy * width + xx]);
				}
			}
			//cout << endl;

			neighborhood[4] = false;
			int black_components = ConnectedComponentsN8(neighborhood);
			neighborhood[4] = true;
			int white_components = ConnectedComponentsN8(neighborhood);
			if (black_components == 1 && white_components == 1)  //yes?
			{
				//calculate its score (Eq85):
				float score = 0;
				{
					int X = (int)((float)x / (float)bigpixel_size);
					int Y = (int)((float)y / (float)bigpixel_size);

					int BlackCount = BlackCounts[Y * Width + X];
					int WhiteCount = bigpixel_size * bigpixel_size - BlackCount;

					//to turn black to white:
					if (H[y * width + x])
					{
						float sign = 0;
						//if the big-pixel's Sigma is 1, sign is -1
						//otherwise sign is +1
						float Zeta = 0;
						if (Sigmas[Y * Width + X])
						{
							sign = -1;
							Zeta = (BlackCount - 1) / (float)(bigpixel_size * bigpixel_size);
						}
						else
						{
							sign = 1;
							Zeta = (WhiteCount + 1) / (float)(bigpixel_size * bigpixel_size);
						}

						score = sign * Zeta;
					}
					//to turn white to black:
					else
					{
						float sign = 0;

						float Zeta = 0;
						if (Sigmas[Y * Width + X])
						{
							sign = 1;
							Zeta = (BlackCount + 1) / (float)(bigpixel_size * bigpixel_size);
						}
						else
						{
							sign = -1;
							Zeta = (WhiteCount - 1) / (float)(bigpixel_size * bigpixel_size);
						}

						score = sign * Zeta;
					}
				}

				//cout << "C " << x << "," << y << " |" << H[y * width + x] << " score:" << score << endl;

				//insert into candidates list in energy highest to lowest order
				bool inserted = false;
				for (list<Candidate>::iterator itr = candidates.begin(); itr != candidates.end(); itr++)
				{
					if (score >= get<2>((*itr)))
					{
						//insert here!
						candidates.insert(itr, make_tuple(Vec2i(x, y), H[y * width + x], score));
						inserted = true;
						break;
					}
				}
				if (!inserted)
					candidates.push_front(make_tuple(Vec2i(x, y), H[y * width + x], score));
			}
			else
			{
				//cout << "nope " << x << "," << y << " bc:" << black_components << " wc:" << white_components << endl;
			}
		}

		//check if H is a correct result with zero digi-errors?
		bool CheckDigiCorrectness(bool* H, int width, int height, int Width, int Height, int bigpixel_size)
		{
			for (int Y = 0; Y < Height; Y++)
			{
				for (int X = 0; X < Width; X++)
				{
					//count the black and white small pixels in this bigpixel
					int num_blacks = 0;
					int num_whites = 0;

					for (int y = Y * bigpixel_size; y < Y * bigpixel_size + bigpixel_size; y++)
					{
						for (int x = X * bigpixel_size; x < X * bigpixel_size + bigpixel_size; x++)
						{
							if (H[y * width + x])
								num_blacks++;
							else
								num_whites++;
						}
					}

					if (num_blacks != 0 && num_whites != 0)
						return false;  //oops!
				}
			}
			return true;
		}
	};
	Worker worker;

	//enumerate candidates for the first time
	//note: just skip mesh-boundary pixels. let's create candidates
	for (int y = 1; y < height - 1; y++)
	{
		for (int x = 1; x < width - 1; x++)
		{
			worker.TryAddCandidate(H, x, y, width, height, Width, Height, bigpixel_size,
				BlackCounts, Sigmas, candidates);
		}
	}

	if(print_debug)
		cout << "init. candidates: " << candidates.size() << " iter_max:" << iter_max << endl;

	//record Sigmas that have been flipped. to early stop
	map<int, int> FlippedSigmas_map;  //key: Y*Width+X, value: number of flips
	int max_flips_threshold = 4;  //maximum allowed # of flips at a bigpixel?
	int num_flips = 0;

	for (int iter = 0; iter < iter_max; iter++)
	{
		if (candidates.size() == 0)
		{
			if(print_debug)
				cout << "candidates empty. failed" << endl;
			return 2;
		}

		Candidate best = candidates.front();
		candidates.pop_front();

		//if (print_debug)
		//{
		//	cout << "iter#" << iter << " best score:" << get<2>(best) <<
		//		" " << get<0>(best) << " " << get<1>(best) << endl;
		//}
		
		//if the best candidate's score is negative, may stop?
		if (get<2>(best)<0)
		{
			//let's check for digi-errors?
			if (!worker.CheckDigiCorrectness(H, width, height, Width, Height, bigpixel_size))
			{
				//oops. let's flip the Sigma of the best candidate's bigpixel and try again
				Vec2i p = get<0>(best); 
				int PX = (int)((float)p.x / (float)bigpixel_size);
				int PY = (int)((float)p.y / (float)bigpixel_size);

				if(print_debug)
					cout << "digi-error! let's flip Sigma@" << PX << "," << PY << endl;	

				//check if # of flips at the bigpixel is too many. then abort
				int num_existing_flips = FlippedSigmas_map[PY * Width + PX];
				if (num_existing_flips >= max_flips_threshold)
				{
					//stop!
					if (print_debug)
						cout << "#flips is already too many. abort." << endl;
					return 2;
				}
				FlippedSigmas_map[PY * Width + PX] = num_existing_flips + 1;

				Sigmas[PY * Width + PX] = !Sigmas[PY * Width + PX];
				num_flips++;

				//update costs of candidates in the bigpixel
				for (list<Candidate>::iterator itr = candidates.begin(); itr != candidates.end(); itr++)
				{
					Vec2i pp = get<0>((*itr));
					int XX = (int)((float)pp.x / (float)bigpixel_size);
					int YY = (int)((float)pp.y / (float)bigpixel_size);
					if (XX == PX && YY == PY)
					{
						//calculate its score again
						float score = 0;
						{
							int BlackCount = BlackCounts[YY * Width + XX];
							int WhiteCount = bigpixel_size * bigpixel_size - BlackCount;

							//to turn black to white:
							if (H[pp.y * width + pp.x])
							{
								float sign = 0;
								//if the big-pixel's Sigma is 1, sign is -1
								//otherwise sign is +1
								float Zeta = 0;
								if (Sigmas[YY * Width + XX])
								{
									sign = -1;
									Zeta = (BlackCount - 1) / (float)(bigpixel_size * bigpixel_size);
								}
								else
								{
									sign = 1;
									Zeta = (WhiteCount + 1) / (float)(bigpixel_size * bigpixel_size);
								}

								score = sign * Zeta;
							}
							//to turn white to black:
							else
							{
								float sign = 0;

								float Zeta = 0;
								if (Sigmas[YY * Width + XX])
								{
									sign = 1;
									Zeta = (BlackCount + 1) / (float)(bigpixel_size * bigpixel_size);
								}
								else
								{
									sign = -1;
									Zeta = (WhiteCount - 1) / (float)(bigpixel_size * bigpixel_size);
								}

								score = sign * Zeta;
							}
						}

						get<2>((*itr)) = score;
					}
				}
				
				continue;
			}

			if (print_debug)
				cout << "digi-correct! done" << endl;
			break;
		}
		
		//carry out the best candidate
		Vec2i p = get<0>(best);
		bool action = get<1>(best);
		if (action)
		{
			H[p.y * width + p.x] = false;
		}
		else
		{
			H[p.y * width + p.x] = true;
		}

		//post steps:
		{
			//big-pixel (X,Y) of the candidate:
			int PX = (int)((float)p.x / (float)bigpixel_size);
			int PY = (int)((float)p.y / (float)bigpixel_size);

			//1. update BlackCounts value of the big-pixel
			if (action)
				BlackCounts[PY * Width + PX] = BlackCounts[PY * Width + PX] - 1;
			else
				BlackCounts[PY * Width + PX] = BlackCounts[PY * Width + PX] + 1;

			//2. update the candidates queue: remove existing candidates within the same bigpixel
			list<Candidate> candidates_new;
			for (list<Candidate>::iterator itr = candidates.begin(); itr != candidates.end(); itr++)
			{
				Vec2i pp = get<0>((*itr));
				int XX = (int)((float)pp.x / (float)bigpixel_size);
				int YY = (int)((float)pp.y / (float)bigpixel_size);

				if (XX != PX || YY != PY)
				{
					candidates_new.push_back((*itr));  //keep the same order
				}
			}
			candidates = candidates_new;
			
			//and create new candidates for every smallpixel in the bigpixel
			for (int yy = PY * bigpixel_size; yy < PY * bigpixel_size + bigpixel_size; yy++)
			{
				for (int xx = PX * bigpixel_size; xx < PX * bigpixel_size + bigpixel_size; xx++)
				{
					worker.TryAddCandidate(H, xx, yy, width, height, Width, Height, bigpixel_size,
						BlackCounts, Sigmas, candidates);
				}
			}
		}
	}

	//save the final H to output
	bool success = true;
	for (int Y = 0; Y < Height; Y++)
	{
		for (int X = 0; X < Width; X++)
		{
			//count the black and white small pixels in this bigpixel
			int num_blacks = 0;
			int num_whites = 0;

			for (int y = Y * bigpixel_size; y < Y * bigpixel_size + bigpixel_size; y++)
			{
				for (int x = X * bigpixel_size; x < X * bigpixel_size + bigpixel_size; x++)
				{
					if (H[y * width + x])
						num_blacks++;
					else
						num_whites++;
				}
			}

			if (num_blacks != 0 && num_whites != 0)
			{
				//digitization error not zero here!
				if (print_debug)
					cout << "oops. digi-err @" << X << "," << Y << " nb:" << num_blacks << " nw:" << num_whites << endl;
				success = false;
				break;
			}

			bool val = num_blacks >= num_whites;
			output[Y * Width + X] = val;
		}
		if (!success)
			break;
	}

	if (!success)
		return 2;  //failure
	else
	{
		//clean success w/o flips?
		if (num_flips == 0)
			return 0;
		else
			return 1;
	}
}

bool DSSpace::DownsamplePassat2022Png(const char* input_filename, int bigpixel_size, bool print_debug)
{
	std::vector<unsigned char> in_buffer; //the raw pixels (RGBA)
	unsigned width = 0, height = 0;
	unsigned error = lodepng::decode(in_buffer, width, height, input_filename);
	if (error)
	{
		cout << input_filename << " lodepng::decode error:" << error << " " << lodepng_error_text(error) << endl;
		return false;
	}

	if ((width % 2) != 0 || (height % 2) != 0)
	{
		cout << "error: width / height not dividable!" << endl;
		return false;
	}

	//turn the image buffer to a binary mask buffer
	bool* mask = new bool[width * height];
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			int offset = (y * width + x) * 4;
			Vec4<unsigned char> p(in_buffer[offset], in_buffer[offset + 1], in_buffer[offset + 2], in_buffer[offset + 3]);

			//binarilization: <threshold = black, >threshold = white
			if (p[0] >= g_ds_png_treshold && p[1] >= g_ds_png_treshold && p[2] >= g_ds_png_treshold)
			{
				mask[y * width + x] = true;
			}
			else
			{
				mask[y * width + x] = false;
			}
		}
	}

	int new_width = width / bigpixel_size;
	int new_height = height / bigpixel_size;
	bool* output = new bool[new_width * new_height];

	int ret = DownsamplePassat2022(width, height, mask, bigpixel_size, output, print_debug);
	cout << "ret: " << ret << endl;
	if (ret != 2)
	{
		//save to a png file!
		std::vector<unsigned char> out_buffer(new_width * new_height * 4);
		for (int y = 0; y < new_height; y++)
		{
			for (int x = 0; x < new_width; x++)
			{
				int offset = (y * new_width + x) * 4;

				if (output[y * new_width + x])
				{
					//white
					out_buffer[offset] = 255;
					out_buffer[offset + 1] = 255;
					out_buffer[offset + 2] = 255;
					out_buffer[offset + 3] = 255;
				}
				else
				{
					//black
					out_buffer[offset] = 0;
					out_buffer[offset + 1] = 0;
					out_buffer[offset + 2] = 0;
					out_buffer[offset + 3] = 255;
				}
			}
		}

		string output_filename = string(input_filename) + "." + to_string(new_width) + "x" + to_string(new_height);
		output_filename += ".Passat.png";

		lodepng::encode(output_filename, out_buffer, new_width, new_height);
	}
	else
	{
		cout << "failed" << endl;
	}

	return ret;
}