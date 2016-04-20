#include "CLC.h"
bool _compare_min_x(Point2f const &p1, Point2f const &p2) { return p1.x < p2.x; }
bool _compare_min_y(Point2f const &p1, Point2f const &p2) { return p1.y < p2.y; }
CLC::CLC(double _fx, double _fy, double _cx, double _cy)
{
    fx = _fx;
    fy = _fy;
    cx = _cx;
    cy = _cy;
    um =Point2d(_cx, _cy);
}
bool Quadrilateral::SortPoints()
{
    //    Point2d min_x = *std::min_element(points.begin(), points.end(), &_compare_min_x);
    //    Point2d min_y = *std::min_element(points.begin(), points.end(), &_compare_min_y);
    //    Point2d max_x = *std::max_element(points.begin(), points.end(), &_compare_min_x);
    //    Point2d max_y = *std::max_element(points.begin(), points.end(), &_compare_min_y);
    if(points.size()!=4){
        std::cout<<"The number of Points is must be 4"<<std::endl;
        return false;
    }
    Point2d center;
    Point2d bottom_l, bottom_r, top_r, top_l;
    bool isFoundBl=false, isFoundBr=false, isFoundTr=false, isFoundTl=false;
    
    center.x = (points[0].x+points[1].x+points[2].x+points[3].x)/4;
    center.y = (points[0].y+points[1].y+points[2].y+points[3].y)/4;
    
    for(int i = 0; i < points.size();i++){
        if(((points[i].x-center.x)<0)&&((points[i].y-center.y)>0)&&(!isFoundBl)){
            bottom_l = points[i];
            isFoundBl = true;
            continue;
        }
        else if(((points[i].x-center.x)>0)&&((points[i].y-center.y)>0)&&(!isFoundBr)){
            bottom_r = points[i];
            isFoundBr = true;
            continue;
        }
        else if(((points[i].x-center.x)>0)&&((points[i].y-center.y)<0)&&(!isFoundTr)){
            top_r = points[i];
            isFoundTr = true;
            continue;
        }
        else if(((points[i].x-center.x)<0)&&((points[i].y-center.y)<0)&&(!isFoundTl)){
            top_l = points[i];
            isFoundTl = true;
            continue;
        }
        else{
            std::cout<<"Point sorting error : it's not quadrilateral."<<std::endl;
            return false;
        }
    }
    points.clear();
    points.push_back(bottom_l);
    points.push_back(bottom_r);
    points.push_back(top_r);
    points.push_back(top_l);
    return true;
}
Point2f CLC::GetIntersectPoint(const Point2f &p1, const Point2f &p2, const Point2f &p3, const Point2f &p4)
{
    // Store the values for fast access and easy
    // equations-to-code conversion
    double x1 = p1.x, x2 = p2.x, x3 = p3.x, x4 = p4.x;
    double y1 = p1.y, y2 = p2.y, y3 = p3.y, y4 = p4.y;
    
    double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    // If d is zero, there is no intersection
    if (d == 0) return Point2f();
    
    // Get the x and y
    double pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
    double x = (pre * (x3 - x4) - (x1 - x2) * post) / d;
    double y = (pre * (y3 - y4) - (y1 - y2) * post) / d;
    
    // Check if the x and y coordinates are within both lines
    //if (x < min(x1, x2) || x > max(x1, x2) ||
    //x < min(x3, x4) || x > max(x3, x4)) return Point2f();
    //if (y < min(y1, y2) || y > max(y1, y2) ||
    //y < min(y3, y4) || y > max(y3, y4)) return Point2f();
    
    // Return the point of intersection
    Point2f ret;
    ret.x = x;
    ret.y = y;
    return ret;
    
}
inline double CLC::GetDistance(const Point2f &p1, const Point2f &p2)
{
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}
bool CLC::SetOffCenteredQuad(vector<Point2f> &points)
{
    if (points.size() == 4 ) {
        //push pack points with initializing param.
        om = Point2f(0,0);
        
        quadOffCentered.points.clear();
        quadCentered.points.clear();
        quadOffCentered.points.push_back(points[0]);
        quadOffCentered.points.push_back(points[1]);
        quadOffCentered.points.push_back(points[2]);
        quadOffCentered.points.push_back(points[3]);
        if(!cv::isContourConvex(quadOffCentered.points))
            return false;
        quadOffCentered.SortPoints();
        return true;
    }
    
    return false;
}


bool CLC::FindProxyQuadrilateral()
{
    om = GetIntersectPoint(quadOffCentered.points[0], quadOffCentered.points[2], quadOffCentered.points[1], quadOffCentered.points[3]);
    
    w0 = GetIntersectPoint(quadOffCentered.points[0], quadOffCentered.points[1], quadOffCentered.points[2], quadOffCentered.points[3]);
    w1 = GetIntersectPoint(quadOffCentered.points[0], quadOffCentered.points[3], quadOffCentered.points[1], quadOffCentered.points[2]);
    
    wd0 = GetIntersectPoint(w0, w1, quadOffCentered.points[0], quadOffCentered.points[2]);
    wd1 = GetIntersectPoint(w0, w1, quadOffCentered.points[1], quadOffCentered.points[3]);
    
    wm = GetIntersectPoint(w0, w1, um, om);
    
    u0 = GetIntersectPoint(quadOffCentered.points[0], wm, um, wd0);
    u2 = GetIntersectPoint(quadOffCentered.points[2], wm, um, wd0);
    u1 = GetIntersectPoint(quadOffCentered.points[1], wm, um, wd1);
    u3 = GetIntersectPoint(quadOffCentered.points[3], wm, um, wd1);
    quadCentered.points.push_back(u0); quadCentered.points.push_back(u1);
    quadCentered.points.push_back(u2); quadCentered.points.push_back(u3);
    return true;
}
bool CLC::CalcCLC(Vector3d &trans, Quaternion<double> &q)
{
    //determinate that the projective quadrilateral is rectangle in real world
    double l0 = sqrt(pow((quadCentered.points[0].x - um.x), 2) + pow((quadCentered.points[0].y - um.y), 2));
    double l2 = sqrt(pow((quadCentered.points[2].x - um.x), 2) + pow((quadCentered.points[2].y - um.y), 2));
    double l1 = sqrt(pow((quadCentered.points[1].x - um.x), 2) + pow((quadCentered.points[1].y - um.y), 2));
    double l3 = sqrt(pow((quadCentered.points[3].x - um.x), 2) + pow((quadCentered.points[3].y - um.y), 2));
    
    double alpha0 = (l0 - l2) / (l0 + l2);
    double alpha1 = (l1 - l3) / (l1 + l3);
    
    double beta = l1 / l0;
    
    bool D = ( (beta >= ( (1-alpha0)/(1+alpha1) )) && (1 >= fabs(alpha1/alpha0) )) || ( (beta <= ( (1-alpha0)/(1+alpha1) )) && (1 <= fabs(alpha1/alpha0) ));
    //cout<<"determinant: "<<D<<endl;
    //    if(!D){
    //        std::cout<<"Couldn't pass the determinant."<<std::endl;
    //        return false;
    //    }
    
    double d = sqrt((pow((1 - alpha1)*beta, 2) - pow(1 - alpha0, 2)) / (pow((1 - alpha1)*alpha0*beta, 2) - pow((1 - alpha0)*alpha1, 2)));
    
    double theta0 = acos(d*alpha0);
    double theta1 = acos(d*alpha1);
    
    double x1 = (quadCentered.points[0].x - quadCentered.points[2].x);
    double x2 = (quadCentered.points[1].x - quadCentered.points[3].x);
    double y1 = (quadCentered.points[0].y - quadCentered.points[2].y);
    double y2 = (quadCentered.points[1].y - quadCentered.points[3].y);
    double rho = acos((x1*x2 +y1*y2) / (sqrt(x1*x1 + y1*y1)*sqrt(x2*x2+y1*y2)));
    
    if (rho < 0)
        rho = M_PI * 2 + rho;
    double phi = acos(cos(theta0)*cos(theta1) + sin(theta0)*sin(theta1)*cos(rho));
    if(std::isnan(phi))
    {
        std::cerr << "Crossing Angle is NaN" << std::endl;
    }
    cout<<"Crossing Angle : "<<phi<<endl;
    d=1;
    Point3d pc(d *cos(theta0)*sin(phi) / sin(phi), -d *cos(theta0)*cos(phi) + cos(theta1) / sin(phi), d *sin(theta0)*sin(theta1)*sin(rho) / sin(phi));
    //cout << "Principle point :\n" << pc << endl;
    
#if 1
    Point2f inputQuad[4];
    // Output Quadilateral or World plane coordinates
    Point2f outputQuad[4];
    
    // Lambda Matrix
    Mat H;
    
    // The 4 points that select quadilateral on the input , from top-left in clockwise order
    // These four pts are the sides of the rect box used as input
    double m = d*l0/((fx+fy)*0.5*sin(theta0)+l0*cos(theta0)) /*, phi = 0.61546*2*/;
    
    inputQuad[0] = Point2f( m, 0 );
    inputQuad[1] = Point2f( m*cos(phi), m*sin(phi) );
    inputQuad[2] = Point2f( -m, 0 );
    inputQuad[3] = Point2f( -m*cos(phi), -m*sin(phi) );
    
    outputQuad[0] = Point2f( quadOffCentered.points[0].x, quadOffCentered.points[0].y );
    outputQuad[1] = Point2f( quadOffCentered.points[1].x, quadOffCentered.points[1].y );
    outputQuad[2] = Point2f( quadOffCentered.points[2].x, quadOffCentered.points[2].y );
    outputQuad[3] = Point2f( quadOffCentered.points[3].x, quadOffCentered.points[3].y );
    
    // Get the Perspective Transform Matrix i.e. lambda
    H = getPerspectiveTransform( inputQuad, outputQuad );
    
    // Intrinsic
    cv::Mat K = (cv::Mat_<double>(3, 3) <<
                 fx,  0, cx,
                 0, fy, cy,
                 0,  0,  1
                 );
    H = K.inv()*H;
    cv::Mat pose;
    pose = Mat::eye(3, 4, CV_32FC1);      // 3x4 matrix, the camera pose
    float norm1 = (float)norm(H.col(0));
    float norm2 = (float)norm(H.col(1));
    float tnorm = (norm1 + norm2) / 2.0f; // Normalization value
    
    
    cv::normalize(H.col(0), pose.col(0));   // Normalize the rotation, and copies the column to pose
    
    cv::normalize(H.col(1),  pose.col(1));   // Normalize the rotation and copies the column to pose
    
    
    Mat p3 = pose.col(0).cross(pose.col(1));   // Computes the cross-product of p1 and p2
    p3.copyTo(pose.col(2));       // Third column is the crossproduct of columns one and two
    
    H.col(2).copyTo(pose.col(3)) /*/ tnorm*/;  //vector t [R|t] is the last column of pose
    pose.col(3) = pose.col(3)/tnorm;
    
    // Map the OpenCV matrix with Eigen:
    
    Eigen::Matrix3d rot;
    rot<< pose.at<float>(0,0), pose.at<float>(0,1), pose.at<float>(0,2),
    pose.at<float>(1,0), pose.at<float>(1,1), pose.at<float>(1,2),
    pose.at<float>(2,0), pose.at<float>(2,1), pose.at<float>(2,2);
    q = rot;
    trans << pose.at<float>(0,3), pose.at<float>(1,3), pose.at<float>(2,3) ;
#endif
#if 0
    ///Perspective-to-Euclidean transformation
    Point2f vecTranslate;
    double t0, t1, s0, s1;
    Point2f us0, us1;
    us0 = GetIntersectPoint(u0, om, um, wd1);
    us1 = GetIntersectPoint(u1, om, um, wd0);
    
    double m=0.1818175;
    
    s0 = GetDistance(us0, um) / l0;
    s1 = GetDistance(us1, um) / l1;
    
    t0 = s0*m*(l0 + l2) / (s0*m*l0 + ((1 - s0)*m+m)*l2);
    t1 = s1*m*(l1 + l3) / (s1*m*l1 + ((1 - s1)*m+m)*l3);
    
    vecTranslate.x = t0 + t1*cos(phi);
    vecTranslate.y = (t1)*sin(phi);
    
    trans[0]=pc.x;
    trans[1]=pc.y;
    trans[2]=pc.z;
    
    Vector3d N1(-pc.x, -pc.y, -pc.z);
    Vector3d N2(0, 0, 1);
    N1.normalize(), N2.normalize();
    Vector3d u = N1.cross(N2);
    u.normalize();
    double alpha = acos(N1.dot(N2)/(N1.norm()*N2.norm()));
    q = AngleAxis<double>(cos(0.5*alpha), sin(0.5*alpha)*u);
    
    //q = Quaternion<double>::FromTwoVectors(N1, N2);
    //(angle, axis.x, axis.y, axis.z);
    q.normalize();
    rot = q.matrix();
    //std::cout << "q: " << q.x() << ", " << q.y() << ", " << q.z() << ", " << q.w() << std::endl;
    std::cout << q.matrix()<<std::endl;
    std::cout << "vector translate:\n" << vecTranslate.x << ", " << vecTranslate.y << std::endl;
    //std::cout << lambda <<std::endl;
#endif
    return true;
}

void CLC::Visualization(cv::Mat &out)
{
    line(out, quadOffCentered.points[0], quadOffCentered.points[1], Scalar(0,255,0));
    line(out, quadOffCentered.points[1], quadOffCentered.points[2], Scalar(0,255,0));
    line(out, quadOffCentered.points[2], quadOffCentered.points[3], Scalar(0,255,0));
    line(out, quadOffCentered.points[3], quadOffCentered.points[0], Scalar(0,255,0));
    line(out, quadOffCentered.points[0], quadOffCentered.points[2], Scalar(0,255,0));
    line(out, quadOffCentered.points[1], quadOffCentered.points[3], Scalar(0,255,0));
    circle(out, quadOffCentered.points[0], 5, cv::Scalar(0, 255, 0), 4, 5);
    circle(out, quadOffCentered.points[1], 5, cv::Scalar(0, 255, 0), 4, 5);
    circle(out, quadOffCentered.points[2], 5, cv::Scalar(0, 255, 0), 4, 5);
    circle(out, quadOffCentered.points[3], 5, cv::Scalar(0, 255, 0), 4, 5);
    circle(out, um, 5, cv::Scalar(255, 0, 255), 4, 5);
    
    line(out, quadCentered.points[0], quadCentered.points[1], Scalar(0,255,255));
    line(out, quadCentered.points[1], quadCentered.points[2], Scalar(0,255,255));
    line(out, quadCentered.points[2], quadCentered.points[3], Scalar(0,255,255));
    line(out, quadCentered.points[3], quadCentered.points[0], Scalar(0,255,255));
    line(out, quadCentered.points[0], quadCentered.points[2], Scalar(0,255,255));
    line(out, quadCentered.points[1], quadCentered.points[3], Scalar(0,255,255));
    circle(out, quadCentered.points[0], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, quadCentered.points[1], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, quadCentered.points[2], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, quadCentered.points[3], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, om, 5, cv::Scalar(255, 255, 0), 4, 5);
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define DEG2RAD 0.017453293f
int thresh = 50, N = 1;
const char* wndname = "Square Detection Demo";
static double angle( Point pt1, Point pt2, Point pt0 )
{
    double dx1 = pt1.x - pt0.x;
    double dy1 = pt1.y - pt0.y;
    double dx2 = pt2.x - pt0.x;
    double dy2 = pt2.y - pt0.y;
    return (dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10);
}
struct LinePolar
{
    float rho;
    float angle;
};


struct hough_cmp_gt
{
    hough_cmp_gt(const int* _aux) : aux(_aux) {}
    bool operator()(int l1, int l2) const
    {
        return aux[l1] > aux[l2] || (aux[l1] == aux[l2] && l1 < l2);
    }
    const int* aux;
};
static void
HoughLinesStandard( const Mat& img, float rho, float theta,
                   int threshold, std::vector<Vec2f>& lines, int linesMax,
                   double min_theta, double max_theta )
{
    int i, j;
    float irho = 1 / rho;
    
    CV_Assert( img.type() == CV_8UC1 );
    
    const uchar* image = img.ptr();
    int step = (int)img.step;
    int width = img.cols;
    int height = img.rows;
    
    if (max_theta < min_theta ) {
        CV_Error( CV_StsBadArg, "max_theta must be greater than min_theta" );
    }
    int numangle = cvRound((max_theta - min_theta) / theta);
    int numrho = cvRound(((width + height) * 2 + 1) / rho);
    
    AutoBuffer<int> _accum((numangle+2) * (numrho+2));
    std::vector<int> _sort_buf;
    AutoBuffer<float> _tabSin(numangle);
    AutoBuffer<float> _tabCos(numangle);
    int *accum = _accum;
    float *tabSin = _tabSin, *tabCos = _tabCos;
    
    memset( accum, 0, sizeof(accum[0]) * (numangle+2) * (numrho+2) );
    
    float ang = static_cast<float>(min_theta);
    for(int n = 0; n < numangle; ang += theta, n++ )
    {
        tabSin[n] = (float)(sin((double)ang) * irho);
        tabCos[n] = (float)(cos((double)ang) * irho);
    }
    
    // stage 1. fill accumulator
    for( i = 0; i < height; i++ )
        for( j = 0; j < width; j++ )
        {
            if( image[i * step + j] != 0 )
                for(int n = 0; n < numangle; n++ )
                {
                    int r = cvRound( j * tabCos[n] + i * tabSin[n] );
                    r += (numrho - 1) / 2;
                    accum[(n+1) * (numrho+2) + r+1]++;
                }
        }
    
    // stage 2. find local maximums
    for(int r = 0; r < numrho; r++ )
        for(int n = 0; n < numangle; n++ )
        {
            int base = (n+1) * (numrho+2) + r+1;
            if( accum[base] > threshold &&
               accum[base] > accum[base - 1] && accum[base] >= accum[base + 1] &&
               accum[base] > accum[base - numrho - 2] && accum[base] >= accum[base + numrho + 2] )
                _sort_buf.push_back(base);
        }
    
    // stage 3. sort the detected lines by accumulator value
    std::sort(_sort_buf.begin(), _sort_buf.end(), hough_cmp_gt(accum));
    
    // stage 4. store the first min(total,linesMax) lines to the output buffer
    linesMax = std::min(linesMax, (int)_sort_buf.size());
    double scale = 1./(numrho+2);
    for( i = 0; i < linesMax; i++ )
    {
        LinePolar line;
        int idx = _sort_buf[i];
        int n = cvFloor(idx*scale) - 1;
        int r = idx - (n+1)*(numrho+2) - 1;
        line.rho = (r - (numrho - 1)*0.5f) * rho;
        line.angle = static_cast<float>(min_theta) + n * theta;
        lines.push_back(Vec2f(line.rho, line.angle));
    }
}

void CLC::findSquares( const Mat& image, vector<vector<Point2f> >& squares ){
    squares.clear();
    
    Mat grayImg, binImg, edgeImg;
    Mat dst;
    image.copyTo(dst);
    cvtColor(image, grayImg, CV_BGR2GRAY);
    threshold(grayImg, binImg, 150, 255, THRESH_OTSU);
    Canny(binImg, edgeImg, 100, 100, 3);
    
    vector<Vec2f> lines;
    //HoughLines(edgeImg, lines, 1, CV_PI/360, 30, 0, 0 );
    
    HoughLinesStandard(edgeImg, 1, CV_PI/360, 30, lines, INT_MAX, 0, CV_PI);

    vector<Point2d> crossPoints;
    
    RNG rng(53242);
    crossPoints.clear();
    for( size_t i = 0; i < lines.size(); i++ )
        for(size_t j =0 ; j < lines.size(); j++)
        {
            float rho0 = lines[i][0], theta0 = lines[i][1];
            Point pt1, pt2, pt3, pt4;
            double a0 = cos(theta0), b0 = sin(theta0);
            double x0 = a0*rho0, y0 = b0*rho0;
            pt1.x = cvRound(x0 + 1000*(-b0));
            pt1.y = cvRound(y0 + 1000*(a0));
            pt2.x = cvRound(x0 - 1000*(-b0));
            pt2.y = cvRound(y0 - 1000*(a0));
            if(i!=j){
                float rho1 = lines[j][0], theta1 = lines[j][1];
                double a1 = cos(theta1), b1 = sin(theta1);
                double x1 = a1*rho1, y1 = b1*rho1;
                pt3.x = cvRound(x1 + 1000*(-b1));
                pt3.y = cvRound(y1 + 1000*(a1));
                pt4.x = cvRound(x1 - 1000*(-b1));
                pt4.y = cvRound(y1 - 1000*(a1));
                //crossPoints.push_back(GetIntersectPoint(Point2d(pt1.x, pt1.y), Point2d(pt2.x, pt2.y), Point2d(pt3.x, pt3.y), Point2d(pt4.x, pt4.y)));
                line( dst, pt3, pt4, Scalar(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255)), 1, CV_AA);
            }
            if(i==0){
                //line( dst, pt1, pt2, Scalar(0,0,255), 1, CV_AA);
            }
        }
    for(size_t i =0 ; i < crossPoints.size(); i++){
        Point tmp;
        tmp.x=cvRound(crossPoints[i].x);
        tmp.y=cvRound(crossPoints[i].y);
        //circle(dst, tmp, 5, cv::Scalar(0, 255, 0), 4, 5);
    }
    imshow("test",dst);
    waitKey(1);
    
}

void CLC::drawSquares( Mat& image, const vector<vector<Point2f> >& squares )
{
    for( size_t i = 0; i < squares.size(); i++ )
    {
        line(image, squares[i][0], squares[i][1], Scalar(0,0,255));
        line(image, squares[i][1], squares[i][2], Scalar(0,0,255));
        line(image, squares[i][2], squares[i][3], Scalar(0,0,255));
        line(image, squares[i][3], squares[i][0], Scalar(0,0,255));
        
        circle(image, squares[i][0], 5, cv::Scalar(0, 255, 0), 4, 5);
        circle(image, squares[i][1], 5, cv::Scalar(0, 255, 0), 4, 5);
        circle(image, squares[i][2], 5, cv::Scalar(0, 255, 0), 4, 5);
        circle(image, squares[i][3], 5, cv::Scalar(0, 255, 0), 4, 5);
    }
}
