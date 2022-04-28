#include "ndt_grid.h"
#include "laserSimulation.h"

#include <chrono>

#include <fstream>

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

	std::ofstream pose_record( "./pose.txt", std::ios::app );
	if( !pose_record.is_open() ){
		std::cerr<<"can not open file !"<<std::endl;
		return 0;
	}

	ndt::NDT ndt;

	slam::simulation::Simulation simulation;
        simulation.openSimulationFile( "../../test_data/laser_data.txt" );

	std::vector<Eigen::Vector3f> keyPoses;

	Eigen::Vector3f robot_pose( 0.0f, 0.0f, 0.0f );

	slam::ScanContainer pre_scan_container;
	bool is_init = false;
	while( !simulation.endOfFile() ){
		slam::sensor::LaserScan scan;
	
                simulation.readAFrameData( scan ); // read the laser data
		
		std::cout<<"-------------- "<<simulation.getFrameCount()<<" -------------"<<std::endl;
		slam::ScanContainer cur_scan_container;
		laserData2Container( scan, cur_scan_container );

		if( !is_init ){
			pre_scan_container = cur_scan_container;

			is_init = true;
		}
		else {
				
			Eigen::Vector3f pose_delta( 0.0f, 0.0f, 0.0f );
			std::cout<<"pre scan container size : "<<pre_scan_container.getSize()<<std::endl;
			ndt.ndtProcess( pre_scan_container, cur_scan_container, pose_delta, 5 );
			std::cout<<"pose delta = "<<std::endl<<pose_delta<<std::endl;
			if( pose_delta[0] < -1 || pose_delta[0] > 1 || pose_delta[1] < -1 || pose_delta[1] > -1 ){
				
				robot_pose += pose_delta;
			}
			std::cout<<"robot pose now = "<<std::endl<<robot_pose<<std::endl;
			pose_record<<robot_pose[0]<<" "<<robot_pose[1]<<" "<<robot_pose[2]<<std::endl;
			pre_scan_container = cur_scan_container;
		}
		
	}


	simulation.closeSimulationFile();
	pose_record.close();
	return 0;
}
