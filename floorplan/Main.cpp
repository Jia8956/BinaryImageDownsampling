#include <stdio.h>
#include <windows.h>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
#include <omp.h>  //OpenMP
#include "ILMBase.h"
#include "downsampling.h"
#include "lodepng.h"

using namespace std;

char g_filename[260] = { NULL };
char g_foldername[260] = { NULL };

//io state from ImGui
bool g_WantCaptureKeyboard = false;
bool g_WantCaptureMouse = false;
bool g_WantCaptureTextInput = false;
//global stuff
bool g_minimize_menu = false;  //minimize imGui menu drawing?

//downsampling stuffs:
int g_ds_bigpixel_width = 4;  //width and height of a big-pixel in terms of small pixels
int g_ds_bigpixel_height = 4;
bool g_ds_save = true;  //save result image?
bool g_ds_print_debug = false;
float g_ds_land_bias_ratio = 0;  //(for UI only) to calculate g_ds_land_bias as % of max-possible value
int g_ds_land_weight = 2;  //w.r.t water_weight=1, land_weight=?
int g_ds_png_treshold = 26;  //threshold for determining a grey pixel is black (0~256)
bool g_ds_alternative_offsets = false;  //final optimal alternative offsets?
int g_ds_h_offset = 0;  //best horizontal and vertical buffer offsets. for visualization
int g_ds_v_offset = 0;
Vec2i g_ds_input_size;  //w,h
vector<bool> g_ds_input;  //saved input binary mask (0 or 1)
Vec2i g_ds_output_size;  //w,h
vector<bool> g_ds_output;  //saved output binary mask
vector<int> g_ds_input_components;  //component indices of input
vector<int> g_ds_output_components;  //component indices of output
//boundaries for debug. each is a VC, type, and solved dist
vector< vector<tuple<Vec2i, int, int>> > g_ds_output_boundaries;  
bool g_ds_save_components_to_file = false;  //output component index map
bool g_ds_local_constraint = true;
int g_ds_neighobrhood_offset = 0;

int main(int argc, char* argv[])
{
	//first, get module folder location	
	if (GetCurrentDirectoryA(sizeof(g_foldername), g_foldername) > 0)
	{
		//add "/" to the end
		strcpy(g_foldername, (string(g_foldername) + "/").c_str());
	}

	//command line mode:
	if (argc > 1)
	{
		//do downsample png
		string input_filename(argv[1]);
		bool to_calculate_errors = false;

		if (argc == 4)
		{
			g_ds_bigpixel_width = std::stoi(argv[2]);
			g_ds_bigpixel_height = std::stoi(argv[3]);
		}
		else if (argc == 5)
		{
			g_ds_bigpixel_width = std::stoi(argv[2]);
			g_ds_bigpixel_height = std::stoi(argv[3]);
			g_ds_land_weight = std::stof(argv[4]);
		}
		else if (argc == 6)
		{
			g_ds_bigpixel_width = std::stoi(argv[2]);
			g_ds_bigpixel_height = std::stoi(argv[3]);
			g_ds_land_weight = std::stof(argv[4]);
			to_calculate_errors = std::stoi(argv[5]);
		}
		else if (argc == 7)
		{
			g_ds_bigpixel_width = std::stoi(argv[2]);
			g_ds_bigpixel_height = std::stoi(argv[3]);
			g_ds_land_weight = std::stof(argv[4]);
			to_calculate_errors = std::stoi(argv[5]);
			g_ds_save_components_to_file = std::stoi(argv[6]);
		}
		else if (argc == 8)
		{
			g_ds_bigpixel_width = std::stoi(argv[2]);
			g_ds_bigpixel_height = std::stoi(argv[3]);
			g_ds_land_weight = std::stof(argv[4]);
			to_calculate_errors = std::stoi(argv[5]);
			g_ds_save_components_to_file = std::stoi(argv[6]);
			g_ds_neighobrhood_offset = std::stoi(argv[7]);
		}

		cout << "filename: " << input_filename << " bigpixel:" << g_ds_bigpixel_width <<
			"X" << g_ds_bigpixel_height << " land_weight:" << g_ds_land_weight << 
			" calculate_errors:" << to_calculate_errors << " save_components:" << g_ds_save_components_to_file << 
			" neighobrhood_offset:" << g_ds_neighobrhood_offset << endl;
		DSSpace::DownsamplePng(input_filename.c_str(), to_calculate_errors);
		return 0;
	}

    return 0;
}