#ifndef __NDT_GRID_H
#define __NDT_GRID_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "scanContainer.h"


#define LASER_DIST 15.0
#define ROWS 30
#define COLUMNS 30
#define HALF_ROWS ROWS / 2
#define HALF_COLUMNS COLUMNS / 2


namespace ndt{

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


template<typename T>
constexpr T laser_max_dist()
{
	return static_cast<T>(LASER_DIST);
}

template<typename T>
constexpr T get_row_range()
{
	return laser_max_dist<T>() / ROWS;
}

template<typename T>
constexpr T get_column_range()
{
	return laser_max_dist<T>() / COLUMNS;
}

template<typename DataType,
                typename = typename std::enable_if<is_Eigen_type<typename DataType::value_type, DataType::RowsAtCompileTime, DataType::ColsAtCompileTime, Eigen::Matrix>::value>::type>
struct GridCell
{
	using type = DataType;
	using value_type = typename DataType::value_type;

	GridCell( const DataType &mean, const DataType &covarince, const int number ) : mean_( mean ), covarince_( covarince ), number_( number )
	{
		
	}
	
	DataType mean_ = 0;
	DataType covarince_ = 0;
	int number_ = 0;
	
	std::vector<DataType> points;

	typename DataType::value_type probablity_;
};


template<typename GridType>
class NdtGrid
{
	using PointType = typename GridType::type;
	using DataType = typename GridType::value_type;
public:
	NdtGrid(  );
	~NdtGrid();

	bool caculateNDTByOneScan( const slam::ScanContainer &scan )
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
					PointType sigma = ( item - it.mean_ ) * ( item - it.mean_ ).transpose();
					it.covarince_ += sigma;
				}
				
				it.covarince_ /= static_cast<DataType>( it.number );
			}
		}
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

private:

	std::vector<GridType> grid = std::vector<GridType>( COLUMNS * ROWS, GridType());
	
};

}

#endif
