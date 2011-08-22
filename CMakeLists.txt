INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#Change PROJECT_NAME to the name of your project
PROJECT(RGBImageColorTransformations)

cmake_minimum_required(VERSION 2.8)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
SET(Required_Packages
ITK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
ITKCommon
ITKIO
ITKBasicFilters
ITKAlgorithms
)

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project
SET(CurrentExe "ImageCompare")
ADD_EXECUTABLE(${CurrentExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
SET(CurrentExe "RGBImageColorTransformations")
ADD_EXECUTABLE(${CurrentExe} RGBImageColorTransformations.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png HSI .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSI_H.png ${CMAKE_SOURCE_DIR}/images/peppers_HSI_H.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSI_S.png ${CMAKE_SOURCE_DIR}/images/peppers_HSI_S.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSI_I.png ${CMAKE_SOURCE_DIR}/images/peppers_HSI_I.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png XYZ .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_XYZ_X.png ${CMAKE_SOURCE_DIR}/images/peppers_XYZ_X.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_XYZ_Y.png ${CMAKE_SOURCE_DIR}/images/peppers_XYZ_Y.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_XYZ_Z.png ${CMAKE_SOURCE_DIR}/images/peppers_XYZ_Z.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png Yuv .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Yuv_Y.png ${CMAKE_SOURCE_DIR}/images/peppers_Yuv_Y.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Yuv_u.png ${CMAKE_SOURCE_DIR}/images/peppers_Yuv_u.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Yuv_v.png ${CMAKE_SOURCE_DIR}/images/peppers_Yuv_v.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png YUV .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_YUV_Y.png ${CMAKE_SOURCE_DIR}/images/peppers_YUV_Y.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_YUV_U.png ${CMAKE_SOURCE_DIR}/images/peppers_YUV_U.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_YUV_V.png ${CMAKE_SOURCE_DIR}/images/peppers_YUV_V.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png HSV .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSV_Y.png ${CMAKE_SOURCE_DIR}/images/peppers_HSV_H.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSV_S.png ${CMAKE_SOURCE_DIR}/images/peppers_HSV_S.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSV_V.png ${CMAKE_SOURCE_DIR}/images/peppers_HSV_V.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png Lab .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Lab_L.png ${CMAKE_SOURCE_DIR}/images/peppers_Lab_L.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Lab_a.png ${CMAKE_SOURCE_DIR}/images/peppers_Lab_a.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Lab_b.png ${CMAKE_SOURCE_DIR}/images/peppers_Lab_b.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png Luv .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Luv_L.png ${CMAKE_SOURCE_DIR}/images/peppers_Luv_L.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Luv_u.png ${CMAKE_SOURCE_DIR}/images/peppers_Luv_u.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_Luv_v.png ${CMAKE_SOURCE_DIR}/images/peppers_Luv_v.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png HSL .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSL_H.png ${CMAKE_SOURCE_DIR}/images/peppers_HSL_H.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSL_S.png ${CMAKE_SOURCE_DIR}/images/peppers_HSL_S.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_HSL_L.png ${CMAKE_SOURCE_DIR}/images/peppers_HSL_L.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png CMY .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMY_C.png ${CMAKE_SOURCE_DIR}/images/peppers_CMY_C.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMY_M.png ${CMAKE_SOURCE_DIR}/images/peppers_CMY_M.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMY_Y.png ${CMAKE_SOURCE_DIR}/images/peppers_CMY_Y.png)
ADD_TEST(RGBImageColorTransformations RGBImageColorTransformations ${CMAKE_SOURCE_DIR}/images/peppers .png CMYK .png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMYK_C.png ${CMAKE_SOURCE_DIR}/images/peppers_CMYK_C.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMYK_M.png ${CMAKE_SOURCE_DIR}/images/peppers_CMYK_M.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMYK_Y.png ${CMAKE_SOURCE_DIR}/images/peppers_CMYK_Y.png)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/images/peppers_base_CMYK_K.png ${CMAKE_SOURCE_DIR}/images/peppers_CMYK_K.png)

#Once you're done modifying this template, you should rename it to "CMakeLists.txt"
