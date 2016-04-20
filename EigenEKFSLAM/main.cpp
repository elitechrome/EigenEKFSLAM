//
//  main.cpp
//  EigenEKFSLAM
//
//  Created by jaemin on 1/26/16.
//  Copyright © 2016 jaemin. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <cstring>
#include <opencv2/opencv.hpp>

#include "EigenEKFSLAM.h"
#include "CLC/CLC.h"

using namespace cv;
using namespace std;
using namespace Eigen;

// helper function:
// finds a cosine of angle between vectors
// from pt0->pt1 and from pt0->pt2


std::vector<cv::Point2f> keypoints1;
int i = 0;

//callback function
void mouseEvent1(int evt, int x, int y, int flags, void* param){
    cv::Mat *src1 = (cv::Mat*)param;
    cv::Point pot;
    //cv::imshow("src1",*src1);
    
    if (evt == CV_EVENT_LBUTTONDOWN && i<4){
        //keypoints1[i].pt.x = x;
        //keypoints1[i].pt.y = y;
        pot = cv::Point(x, y);
        cv::circle(*src1, pot, 5, cv::Scalar(0, 255, 0), 4, 5);
        
        keypoints1.push_back(cv::Point(x, y));
        printf("사각형의 %d번째 꼭지점의 좌표(%d, %d)\n", i + 1, x, y);
        cv::imshow("Image", *src1);
        i++;
    }
}

int main(int argc, char** argv)
{
    //input parameters of CLC : (fx, fy, cx, cy)
    CLC clc(300,300,512,384);
    std::vector<vector<Point2f> > squares;
    cv::Mat image;
    argv[1] ="/Users/jaemin/Dropbox/EigenEKFSLAM/EigenEKFSLAM/IMG_2658.JPG";
    argv[2] ="/Users/jaemin/Dropbox/EigenEKFSLAM/EigenEKFSLAM/IMG_2660.JPG";
    argv[3] ="/Users/jaemin/Dropbox/EigenEKFSLAM/EigenEKFSLAM/IMG_2664.JPG";
    argv[4] ="/Users/jaemin/Dropbox/EigenEKFSLAM/EigenEKFSLAM/IMG_2658.JPG";
    argv[5] =0;

    EigenEKFSLAM ekf_clc(Eigen::Vector3d(0, 0, 0), Eigen::Quaternion<double>(1, 0, 0, 0));
    
    for( int index = 1; argv[index] != 0; index++ )
    {
        cv::Mat original_image = cv::imread(argv[index], 1);
        if( original_image.empty() ){
            std::cout << "Couldn't load " << argv[index] << std::endl;
            return -1;
        }
        cv::resize(original_image, original_image, Size(1024, 768));
        //1. EKF Prediction
        double dt = 0.1;
        ControlVector tmpControl;
        tmpControl << 0,0,0,1,0,0,0;
        ekf_clc.aprioriEstimate(dt, tmpControl);
        
        //2. Rectangle Selection
        std::vector<MeasurementVector> vecLandmarks;
        char inputFlag;
        do{
            keypoints1.clear();
            i = 0;
            original_image.copyTo(image);
            cv::namedWindow("Image", CV_WINDOW_AUTOSIZE);
            cv::imshow("Image", image);
                cv::setMouseCallback("Image", mouseEvent1, &image);

            inputFlag = cv::waitKey();
            if (keypoints1.empty()){
                std::cout << "error, no points are selected.\n";
                i=0;
                continue;
            }
            
            if (inputFlag == 'd'){
                continue;
            } else{
        //3. CLC Pose Calculation for each rectangle
                int ID_rect=0;
                cv::imshow("Image", original_image);
                std::cout << "Input ID of this rectangle : ";
                std::cin >> ID_rect;
                clc.SetOffCenteredQuad(keypoints1);
                clc.FindProxyQuadrilateral();
                Vector3d trans; Quaternion<double> q;
                if(!clc.CalcCLC(trans, q))
                {
                    std::cerr << "CLC is NaN" << std::endl;
                    continue;
                }
                MeasurementVector tmpLandmark(8);
                tmpLandmark << 0,0,0,q.w(),q.x(),q.y(),q.z(),ID_rect;
                tmpLandmark.block(0,0,3,1) = trans;
                std::cout << tmpLandmark << std::endl;
                vecLandmarks.push_back(tmpLandmark);
                clc.Visualization(image);
                cv::imshow("Image", image);
                cv::waitKey(1);
                inputFlag = cv::waitKey();
            }
        }while(inputFlag != 'f');
        inputFlag = 0;
        //4. EKF Correction
        ekf_clc.measurementEstimate(vecLandmarks);
        
        std::cout << "EKF state :" << ekf_clc.getState().transpose() << std::endl;
        
    }
    return 0;
}