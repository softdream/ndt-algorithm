#ifndef __NDT_GRID_H
#define __NDT_GRID_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <type_traits>

#include "scanContainer.h"


#define LASER_DIST 15.0
#define ROWS 30
#define COLUMNS 30
#define HALF_ROWS ROWS / 2
#define HALF_COLUMNS COLUMNS / 2


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
	return laser_max_dist<T>() / ROWS;
}

template<typename T, typename = typename std::enable_if<is_double_or_float<T>::value>::type>
static constexpr T get_column_range()
{
	return laser_max_dist<T>() / COLUMNS;
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

	GridCell( const DataType &mean, const DataType &covarince, const int number ) : mean_( mean ), covarince_( covarince ), number_( number )
	{
		
	}
	
	DataType mean_ = 0;
	covarince_type covarince_ = 0;
	int number_ = 0;
	
	std::vector<DataType> points;

	typename DataType::value_type probablity_;
};

template<typename T>
using GridCellType = GridCell<Eigen::Matrix<T, 2, 1>>;

template<typename T, template<template U> class GridType = GridCellType>
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
		if( first_scan.empty() || second_scan.empty() ){
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

private:
	void caculateNDTByFirstScan( const slam::ScanContainer &scan )
	{
		for( size_t i = 0; i < scan.getSize(); i ++ ){
			PointType point = scan.getIndexData( i );
			int index = pointMapToGrid( point );
			
			grid[index].number_ ++;
			grid[index].mean_ ++;
			grid[index].points.push_back( point );
		}	

		for( auto it : grid ){
			if( it.number >= 3 ){
				DataType average = it.mean_ / static_cast<DataType>( it.number );
				it.mean_ = average;
				for( auto item : it.points ){
					CovarinceType sigma = ( item - it.mean_ ) * ( item - it.mean_ ).transpose();
					it.covarince_ += sigma;
				}
				
				it.covarince_ /= static_cast<DataType>( it.number );
			}
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
				
			PointType e = point_in_first_frame - grid[index].mean_;
			CovarinceType sigma_inverse = ( grid[index].covarince_ ).inverse();

			Eigen::Matrix<DataType, 2, 3> Jacobian;
			Jacobian << 1, 0, -point_in_first_frame[0] * ::sin( point_in_first_frame[2] ) - point_in_first_frame[1] * ::cos( point_in_first_frame[2] ),
				    0, 1,  point_in_first_frame[0] * ::cos( point_in_first_frame[2] ) - point_in_first_frame[1] * ::sin( point_in_first_frame[2] );

			b += ( e.transpose() * sigma_inverse * Jacobian ).transpose();
			H += Jacobian.transpose() * sigma_inverse * Jacobian;			
		}
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
		DataType x = point[0] / get_column_range<DataType>();
		DataType y = point[1] / get_row_range<DataType>();

		if( point[1] >= 0 ){
			if( point[0] >= 0 ){
				return ( HALF_ROWS - y ) * 10 + ( x + HALF_COLUMNS - 1 );
			}
			else {
				return ( HALF_ROWS - y ) * 10 + ( 0 - x - 1 );
			}
		}
		else {
			if( point[0] >= 0 ){
				return ( HALF_ROWS - y - 1 ) * 10 + ( x + HALF_COLUMNS - 1 );
			}
			else {
				return ( HALF_ROWS - y - 1 ) * 10 + ( 0 - x - 1 );
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
	
		return rotate * point_old + delta.head<2>();
	}

private:

	std::vector<GridType<T>> grid = std::vector<GridType<T>>( COLUMNS * ROWS, GridType());
	
	Eigen::Matrix<DataType, 3, 3> H;
        Eigen::Matrix<DataType, 3, 1> b;
};

}

#endif
