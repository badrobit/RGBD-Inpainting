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

    /**
    namedWindow( "image", 1 );

    img = img0.clone();
    inpaintMask = Mat::zeros(img.size(), CV_8U);

    imshow("image", img);
    setMouseCallback( "image", onMouse, 0 );

    for(;;)
    {
        char c = (char)waitKey();

        if( c == 27 )
            break;

        if( c == 'r' )
        {
            inpaintMask = Scalar::all(0);
            img0.copyTo(img);
            imshow("image", img);
        }

        if( c == 'i' || c == ' ' )
        {
            Mat inpainted;
            inpaint(img, inpaintMask, inpainted, 3, INPAINT_TELEA);
            imshow("inpainted image", inpainted);
        }
    }

	if( event == EVENT_LBUTTONUP || !(flags & EVENT_FLAG_LBUTTON) )
        prevPt = Point(-1,-1);
    else if( event == EVENT_LBUTTONDOWN )
        prevPt = Point(x,y);
    else if( event == EVENT_MOUSEMOVE && (flags & EVENT_FLAG_LBUTTON) )
    {
        Point pt(x,y);
        if( prevPt.x < 0 )
            prevPt = pt;
        line( inpaintMask, prevPt, pt, Scalar::all(255), 5, 8, 0 );
        line( img, prevPt, pt, Scalar::all(255), 5, 8, 0 );
        prevPt = pt;
        imshow("image", img);
    }
	*/
}
