CMAKE_MINIMUM_REQUIRED( VERSION 3.21.1 )

ADD_DEFINITIONS( -std=c++14 )

SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" )

#find required eigen
FIND_PACKAGE( Eigen3 )

#directory of eigen headers
INCLUDE_DIRECTORIES( ${EIGEN3_INCLUDE_DIR} )

#find required opencv
FIND_PACKAGE( OpenCV REQUIRED )

#directory of opencv headers
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/src src_file )
#AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR} src_file )

INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/include" )

#directory of opencv link libs
LINK_DIRECTORIES(${OpenCV_LIBRARY_DIRS})

# 找到所有含有mian函数的cpp文件，保存在mains里
file(
    GLOB_RECURSE mains RELATIVE "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/ndt_test.cpp"
    "${PROJECT_SOURCE_DIR}/ndt_laser_odometry.cpp"
)

# 循环遍历每个含main的文件，分别编译生成exe文件
foreach(mainfile IN LISTS mains)
    # Get file name without directory
    get_filename_component(mainname ${mainfile} NAME_WE)
    add_executable(${mainname} ${mainfile} ${src_file})
    # 添加库文件目录
    target_link_libraries(${mainname} ${OpenCV_LIBS})# 链接库文件
endforeach()

