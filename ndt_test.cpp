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

	Eigen::Vector3f p( 0, 0, 0 );

	
	ndt.ndtProcess( scanContainer1, scanContainer2, p );

	std::cout<<"result p : "<<std::endl<<p<<std::endl;	

	return 0;
}
