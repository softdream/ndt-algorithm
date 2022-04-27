#ifndef __NDT_GRID_H
#define __NDT_GRID_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <type_traits>

#include "scanContainer.h"


#define LASER_DIST 14.0
#define ROWS 28
#define COLUMNS 28
#define HALF_ROWS 14
#define HALF_COLUMNS 14


namespace ndt{

template<typename T>
struct is_double_or_float
{
        static const bool value = false;
};

template<>
struct is_double_or_float<float>
{
        static const bool value = true;
};

template<>
struct is_double_or_float<double>
{
        static const bool value = true;
};

template<typename T, int Rows, int Cols,
        template<typename U, int R, int C, int Option, int MaxR, int MaxC>
        class EigenType>
struct is_Eigen_type
{
        static const bool value = false;
};

template<typename T, int Rows, int Cols>
struct is_Eigen_type<T, Rows, Cols, Eigen::Matrix>
{
        using type = Eigen::Matrix<T, Rows, Cols>;
        static const bool value = true;
};

template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static constexpr T CONST_PI()
{
        return static_cast<T>( 3.14159265358979323846 );
}

template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static constexpr T CONST_TWO_PI()
{
	return CONST_PI<T>() * 2.0;
}


template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static constexpr T laser_max_dist()
{
	return static_cast<T>( LASER_DIST );
}

template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static constexpr T get_row_range()
{
	return laser_max_dist<T>() / HALF_ROWS;
}

template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static constexpr T get_column_range()
{
	return laser_max_dist<T>() / HALF_COLUMNS;
}

template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static void angleNormalize( T &angle )
{
	if( angle >= CONST_PI<T>() ) {
		angle -= CONST_TWO_PI<T>();
	}

	if( angle <= -CONST_PI<T>() ){
		angle += CONST_TWO_PI<T>();
	}
}


template<typename DataType,
                typename = typename std::enable_if<is_Eigen_type<typename DataType::value_type, DataType::RowsAtCompileTime, DataType::ColsAtCompileTime, Eigen::Matrix>::value>::type>
struct GridCell
{
	using type = DataType;
	using value_type = typename DataType::value_type;
	using covarince_type = typename Eigen::Matrix<value_type, 2, 2>;
	
	GridCell()
	{

	}
	
	~GridCell()
	{

	}

	GridCell( const DataType &mean, const DataType &covarince, const int number ) : mean_( mean ), covarince_( covarince ), number_( number )
	{
		
	}
	
	DataType mean_ = DataType::Zero();;
	covarince_type covarince_ = covarince_type::Zero();
	int number_ = 0;
	
	std::vector<DataType> points;

	typename DataType::value_type probablity_ = 0; 
};

template<typename T>
using GridCellType = GridCell<Eigen::Matrix<T, 2, 1>>;

template<typename T, template<typename U> class GridType = GridCellType>
class NdtGrid
{
public:
	using PointType = typename GridType<T>::type;
	using DataType = typename GridType<T>::value_type;
	using PoseType = typename Eigen::Matrix<DataType, 3, 1>;
	using CovarinceType = typename GridType<T>::covarince_type;

	NdtGrid() {  }
	~NdtGrid() {  } 

	bool ndt_process( const slam::ScanContainer &first_scan,
			  const slam::ScanContainer &second_scan,
			  PoseType &p,
			  const int max_iterations = 10 )
	{
		if( first_scan.isEmpty() || second_scan.isEmpty() ){
			return false;
		}	

		caculateNDTByFirstScan( first_scan );
		
		int iteration = 0;
		for( ; iteration < max_iterations; iteration ++ ){
			estimateTransformationOnce( second_scan, p );
		}
		
		angleNormalize( p[2] );
		std::cout<<"estimated p = "<<std::endl<<p<<std::endl;		

		return true;
	}	

//private:
	void caculateNDTByFirstScan( const slam::ScanContainer &scan )
	{
		for( size_t i = 0; i < scan.getSize(); i ++ ){
			//std::cout<<"--------------"<<std::endl;
			PointType point = scan.getIndexData( i );
			int index = pointMapToGrid( point );
			//std::cout<<"point "<<i<<": "<<std::endl<<point<<std::endl;
			//std::cout<<"index : "<<index<<std::endl;	
			grid[index].number_ ++;
			grid[index].mean_ += point;
			grid[index].points.push_back( point );
		}	

		for( size_t i = 0; i < grid.size(); i ++ ){
			if( grid[i].number_ >= 3 ){
				PointType average = grid[i].mean_ / static_cast<DataType>( grid[i].number_ );
				grid[i].mean_ = average;
				for( auto item : grid[i].points ){
					CovarinceType sigma = ( item - grid[i].mean_ ) * ( item - grid[i].mean_ ).transpose();
					grid[i].covarince_ += sigma;
				}
				
				grid[i].covarince_ /= static_cast<DataType>( grid[i].number_  - 1 );
				//std::cout<<"covarince : "<<std::endl<<it.covarince_<<std::endl;
			}
			else {
				//grid[i].mean_ = PointType( 0, 0 );	
				CovarinceType cov;
				cov << 65536, 65536, 65536, 65536;
				grid[i].covarince_ = cov;
			}
			//std::cout<<"covarince : "<<it.number_<<", i : "<<i<<std::endl<<it.covarince_<<std::endl;
		}
		
	}

	void getHessianDerived( const slam::ScanContainer &scan, 
	    	           const PoseType &p,
			   Eigen::Matrix<DataType, 3, 3> &H,
			   Eigen::Matrix<DataType, 3, 1> &b )
	{
		H = Eigen::Matrix<DataType, 3, 3>::Zero();
		b = Eigen::Matrix<DataType, 3, 1>::Zero();

		for( size_t i = 0; i < scan.getSize(); i ++ ){
			// for all points in second scan frame, transform
                        PointType point = scan.getIndexData( i );
			PointType point_in_first_frame = pointCoordinateTransform( point, p );                
			
		        int index = pointMapToGrid( point_in_first_frame );
			//std::cout<<"index : "<<index<<std::endl;				

			PointType e = point_in_first_frame - grid[index].mean_;
			CovarinceType sigma_inverse = ( grid[index].covarince_ ).inverse();
			//std::cout<<"e : "<<std::endl<<e<<std::endl;
			//std::cout<<"sigma inverse : "<<std::endl<<sigma_inverse<<std::endl;		

			DataType tmp1 = -point[0] * ::sin( p[2] ) - point[1] * ::cos( p[2] );
			DataType tmp2 =  point[0] * ::cos( p[2] ) - point[1] * ::sin( p[2] );
			Eigen::Matrix<DataType, 2, 3> Jacobian;
			Jacobian << 1, 0,  tmp1,
				    0, 1,  tmp2;

			std::cout<<"Jacobian : "<<std::endl<<Jacobian<<std::endl;
			//std::cout<<"b = "<<std::endl<<( e.transpose() * sigma_inverse * Jacobian ).transpose();

			b += ( e.transpose() * sigma_inverse * Jacobian ).transpose();
			H += Jacobian.transpose() * sigma_inverse * Jacobian;			
			//std::cout<<"b : "<<std::endl<<b<<std::endl;
		}
		
		std::cout<<"b : "<<std::endl<<b<<std::endl;
                std::cout<<"H : "<<std::endl<<H<<std::endl;
	}	

	void estimateTransformationOnce( const slam::ScanContainer &scan,
                           		 PoseType &p )	
	{
		getHessianDerived( scan, p, H, b );

		PoseType delta_p = -H.inverse() * b;
		
		p += delta_p; 	
	}

private:	
	const int pointMapToGrid( const PointType &point ) const
	{
		int x = ( point[0] / get_column_range<DataType>() );
		int y = ( point[1] / get_row_range<DataType>() );
	
		if( point[1] >= 0 ){
			if( point[0] >= 0 ){
				return ( HALF_ROWS - y - 1 ) * COLUMNS + ( x + HALF_COLUMNS + 1 );
			}
			else {
				return ( HALF_ROWS - y - 1 ) * COLUMNS + ( x + HALF_COLUMNS );
			}
		}
		else {
			if( point[0] >= 0 ){
				return ( HALF_ROWS - y ) * COLUMNS + ( HALF_COLUMNS + x ) + 1;
			}
			else {
				return ( HALF_ROWS - y ) * COLUMNS + ( HALF_COLUMNS + x );
			}
		}
	}
	
	const DataType caculateProbability( const PointType &x, const PointType &mean, const PointType &covarince ) const
	{
		return ::exp( -( ( x - mean ).transpose() * covarince.inverse() * ( x - mean ) ) * 0.5 );
	}	

	const PoseType poseCoordinateTransform( const PoseType &pose_old, const PoseType &delta ) const
	{
		Eigen::Matrix<DataType, 3, 3> trans;
		trans << ::cos( delta[2] ), -::sin( delta[2] ), delta[0],
			 ::sin( delta[2] ),  ::cos( delta[2] ), delta[1],
				0	  ,	    0	      ,    1	;
	
		return trans * pose_old;
	}

	const PointType pointCoordinateTransform( const PointType &point_old, const PoseType &delta ) const
	{
		Eigen::Matrix<DataType, 2, 2> rotate;
		rotate << ::cos( delta[2] ), -::sin( delta[2] ), 
			  ::sin( delta[2] ),  ::cos( delta[2] );
	
		Eigen::Matrix<DataType, 2, 1> trans( delta[0], delta[1] );
		return rotate * point_old + trans;
	}

//private:

public:
	std::vector<GridType<T>> grid = std::vector<GridType<T>>( COLUMNS * ROWS + 1 );
	
	Eigen::Matrix<DataType, 3, 3> H;
        Eigen::Matrix<DataType, 3, 1> b;
};

using NDT = NdtGrid<float>;

}

#endif


