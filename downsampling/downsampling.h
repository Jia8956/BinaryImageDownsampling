#pragma once

using namespace std;

namespace DSSpace
{

	//label topology (i.e., lands (=islands) and water (= sea and lakes) of a given mask
	//we will mark every pixel of the mask as:
	//{type (0=land, 1=water), index (# of island; # of sea (0) or lakes (>0)) }
	//num_lands = "h0", num_waters = "h1" for betti numbers
	//labels (result): key: row-major index, value: pos and label
	//Eulers: each component's #V, #E, #F
	//is_boundary_flags: "touching-boundary" flags of components 
	void LabelTopology(int width, int height, bool* mask, unordered_map<int, pair<Vec2i, pair<bool, int>>>& labels,
		int& num_lands, int& num_waters, vector<tuple<int,int,int>>& Eulers, vector<bool>& is_boundary_flags);

	//simpler version w/o Eulers calculation etc
	void LabelTopology2(int width, int height, bool* mask, unordered_map<int, pair<Vec2i, pair<bool, int>>>& labels,
		int& num_lands, int& num_waters);

	//enumerate land-water boundaries given a label map
	//boundaries: <list of land component-index and water component index, and boundary size> 
	bool EnumerateBoundaries(int width, int height, unordered_map<int, pair<Vec2i, pair<bool, int>>>& labels,
		int num_lands, int num_waters, vector<pair<Vec2i,int>> &boundaries);

	//load a black-and-white png file, down downsampling, save the result to another png file
	//output filename = input_file.WxH.png (e.g., "input.png.256x256.png")
	//png_threshold: upper threshold (out of 256) for a pixel to be considered as black   
	bool DownsamplePng(const char* input_filename, bool calculate_error_metrics);

	//(old way) topology preserving downsampling of a binary mask of a grid
	//we turn every bigpixel_width X bigpixel_height small pixels into a big pixel
	//this means that width and height should be dividable by them
	//return: 0=success, 2=failure by infeasible, 1=other failures
	int DownsampleByEuler(int width, int height, bool* mask/*size = width*height */,
		int bigpixel_width, int bigpixel_height, bool* output);

	//using new boundary-based topology constraints
	int Downsample(int width, int height, bool* mask/*size = width*height */,
		int bigpixel_width, int bigpixel_height, bool* output);

	//calculate error metrics between an input buffer and an output (smaller) buffer
	//width should be dividable by new_width, etc
	bool ErrorMetrics(int width, int height, bool* input/*size = width*height */,
		int new_width, int new_height, bool* output/*size = new_width*new_height */,
		float &IoU, float &Dice, float &Precision, float &Recall);
	//file-input version:
	bool ErrorMetricsPng(const char* input_filename, const char* output_filename, 
		float& IoU, float& Dice, float& Precision, float& Recall);

	//calculate pixel weights for land and water pixels, respectively
	//total_pixels: # pixels in a big-pixel. threshold_land_count: the amount of land pixels to "just" win
	//e.g., in a 16-bit pixel, threshold = 8, then 8 land pixels win 8 water pixels, but 7 land pixels shall lose 9 water pixels
	bool CalculatePixelWeights(int total_pixels, int threshold_land_count, int& land_weight, int& water_weight);

	//binay image thinning 
	//https://homepages.inf.ed.ac.uk/rbf/HIPR2/thin.htm
	bool Thinning(int width, int height, bool* input/*size = width*height*/, bool *output/*size = width*height */);

	//fill a smallest hole
	bool FillHole(int width, int height, bool *mask/*size = width*height*/);

	//ADAPTIVE CROSSING NUMBERS AND THEIR APPLICATION TO BINARY DOWNSAMPLING
	//always do 2x2 -> 1x1 downsampling
	bool DownsampleACN(int width, int height, bool* mask/*size = width*height */, 
		bool* output/*size = (width/2)*(height/2) */);
	//do many levels of DownsampleACN to an input png file. until 1x1 or failure
	bool DownsampleACNPng(const char* input_filename);

	//Passat2022 (Homotopic affine transformations in the 2D Cartesian Grid)
	//return: 0=success. 1=success after "sigma-flips". 2=failed
	int DownsamplePassat2022(int width, int height, bool* mask/*size = width*height */,
		int bigpixel_size, bool* output/*size = (width/bsize)*(height/bsize) */, bool print_debug);
	bool DownsamplePassat2022Png(const char* input_filename, int bigpixel_size, bool print_debug);
}