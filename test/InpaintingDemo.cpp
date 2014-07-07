#include <opencv2/highgui/highgui.hpp>

#include <iostream>

using namespace cv; 
using namespace std; 

int main( int argc, char** argv )
{
    string filename = argc >= 2 ? argv[1] : "data/test.png"; 
    Mat original_image = imread( filename, -1 ); 
    if( original_image.empty() )
    {
    	cerr << "Could not load image " << filename << ".\n" << endl; 	
    	cerr << "\nProper Usage:\n\n \t\t\t InpaintingDemo <image_name>" << endl; 
    	return 1; 
    }
}
