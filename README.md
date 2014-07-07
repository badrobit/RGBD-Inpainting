RGBD-Inpainting
===============
This module handles the reconstruction of an RGB-D point clould. This software is meant as the last and final step in the NaN-Clustering/inpainting project. It will provide the following functionality: 

## Functionality

1. Creating a depth gradient map from the original depth map. 
2. Inpainting the depth map based on all provided NaN-Cluster Masks. 
3. Reconstruction of a pcl::PointCloud<pcl::XYZRGB>. 

## Required Input

1. A Depth Map. 
2. An RGB image that corresponds to the Depth Map which is provided. 
 * **Note:** The Depth Map & RGB image **MUST** be the same size! 
3. A vector of all areas in the Depth Map which should be infilled. 
 * **Note:** The main focus of this is to support my masters thesis research and as such it will only perform the inpainting on the depth map. If you also need to perform inpainting operations on the RGB image you need to add this functionality yourself. 

## Output

1. An organized PCL point cloud. 

### Debugging Output

If the program is compiled in a debugging mode you will also get the following output files: 

1. Depth Map + Inpainted portion at each stage of the inpainting process. 
2. Source Patch selections highlighted on the original Depth Map for visual inspection.

## Compiling

Everything you need to compile and install the source code is provided in install.sh. You only need to modify the permissions and then run the file: 

```
chmod a+x install.sh
sh install.sh
```

## Running the program

## Running tests.
