#include "ndt_grid.h"
#include "laserSimulation.h"

void laserData2Container( const slam::sensor::LaserScan &scan, slam::ScanContainer &container )
{
        size_t size = 1440;

        float angle = -3.14159f;
        container.clear();

        for( int i = 0; i < size; i ++ ){
                float dist = scan.ranges[ i ];

                if( dist >= 0.0099999998f && dist <= 14.0000000000f ){
                        //dist *= scaleToMap;
                        container.addData( Eigen::Vector2f( cos(angle) * dist, sin(angle) * dist ) );
                }

                angle += 0.0043633231f;
        }

        std::cout<<"Scan Container Size: "<<container.getSize()<<std::endl;
}

void showTwoScanFrame( const slam::ScanContainer &container1, const slam::ScanContainer &container2, const Eigen::Vector3f &delta, const float scale = 20 )
{
	cv::Mat image = cv::Mat::zeros( 1200, 1200, CV_8UC3 );

	cv::Point2d center( 600, 600 );
        cv::circle(image, center, 1, cv::Scalar(0, 255, 0), 1);
        cv::line( image, cv::Point( 600, 0 ), cv::Point( 600, 1200 ), cv::Scalar( 67, 128, 94 ), 1 );
        cv::line( image, cv::Point( 0, 600 ), cv::Point( 1200, 600 ), cv::Scalar( 67, 128, 94 ), 1 );

	for( int i = 0; i < container1.getSize(); i ++ ){
		cv::Point2d point( container1.getIndexData(i)[0] * scale + 600, container1.getIndexData(i)[1] * scale + 600 );
                cv::circle(image, point, 1, cv::Scalar(0, 0, 255), -1);
		//std::cout<<"x: "<<container1.getIndexData(i)[0]<<", y: "<<container1.getIndexData(i)[1]<<std::endl;
	}

	for( int i = 0; i < container2.getSize(); i ++ ){
		Eigen::Matrix<float, 2, 2> rotate;
		rotate << ::cos( delta[2] ), -::sin( delta[2] ),
                          ::sin( delta[2] ),  ::cos( delta[2] );		
		Eigen::Matrix<float, 2, 1> trans( delta[0], delta[1] );
		
		Eigen::Vector2f point_old( container2.getIndexData(i)[0], container2.getIndexData(i)[1] );
		Eigen::Vector2f point_new = rotate * point_old + trans;

                cv::Point2d point( point_new[0] * scale + 600, point_new[1] * scale + 600 );
                cv::circle(image, point, 1, cv::Scalar(0, 255, 0), -1);
        }

	cv::imshow( "scan", image );
        cv::waitKey(0);
}

int main()
{
	std::cout<<"--------------- NDT TEST ---------------"<<std::endl;	

	ndt::NDT ndt;

	slam::simulation::Simulation simulation1, simulation2;
        simulation1.openSimulationFile( "frame1.txt" );
	simulation2.openSimulationFile( "frame2.txt" );

	slam::sensor::LaserScan scan1, scan2;

        simulation1.readAFrameData( scan1 );
	simulation2.readAFrameData( scan2 );

	slam::ScanContainer scanContainer1, scanContainer2;
	laserData2Container( scan1, scanContainer1 );
	laserData2Container( scan2, scanContainer2 );

	Eigen::Vector3f p( 0.2, -0.3, 0.09 );
	showTwoScanFrame( scanContainer1, scanContainer2, p );

	
	ndt.caculateNDTByFirstScan( scanContainer1 );

	Eigen::Matrix<float, 3, 3> H;
        Eigen::Matrix<float, 3, 1> b;

	ndt.getHessianDerived( scanContainer2, p, H, b );
	//ndt.estimateTransformationOnce( scanContainer2, p );

	std::cout<<"result p : "<<std::endl<<p<<std::endl;	

	return 0;
}
