# Topology-Preserving Binary Image Downsampling

## How to use topology-preserving downsampling tool

1. Open the "downsampling" project, and select the "x64" platform and "Release" configurtion.
2. Build the project, and the "downsampling,exe" excutable will be generated in the "x64/Release" folder.
3. Run the fllowing command to downsample an image.

   ```
   downsampling.exe <image_filename> [<bigpixel_width> <bigpixel_height> <[land_weight> [<calculate_errors> [<save_components> [<neighborhood_offset>]]]]]
   ```
   ### Parameter Discription
   - <input_filename>: Path to the input binary image file (required)
   - <bigpixel_width> and <bigpixel_height>: Width and height of a bigpixel (optional, integer)
     - downsampling factor
     - Default value: 4x4
   - <land_weight>: Land (Foreground) weight (optional, integer)
     - Default value: 2
   - <calculate_errors>: Whether to calculate errors (opptional, 0 or 1)
     - Default value: 0 (do not calculate)
   - <save_components>: Whether to output component index map (optional, 0 or 1)
     - Default value: 0 (do not output)
   - <neighborhood_offset>: Coverage offset
      - Default: 0


<div style="display: flex; justify-content: center;">
  <table>
    <tr>
      <td><img src="https://github.com/chenchiachia/BinaryImageDownsampling/blob/main/downsampling/images/BinaryMask/19_112.png" alt="Image 1" width="200"></td>
      <td><img src="https://github.com/chenchiachia/BinaryImageDownsampling/blob/main/downsampling/images/BinaryMask/19_112.png.64x64.png" alt="Image 2" width="200"></td>
    </tr>
    <tr>
      <td>Original binary image</td>
      <td>8x dowwnsampling results</td>
    </tr>
  </table>
</div>
