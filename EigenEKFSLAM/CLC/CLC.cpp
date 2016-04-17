#include "CLC.h"
bool _compare_min_x(Point2f const &p1, Point2f const &p2) { return p1.x < p2.x; }
bool _compare_min_y(Point2f const &p1, Point2f const &p2) { return p1.y < p2.y; }
CLC::CLC()
{
    um =Point2d(512,384);
    thresh = 40;
    N = 1;
}
bool Quadrilateral::SortPoints()
{
    Point2d min_x = *std::min_element(points.begin(), points.end(), &_compare_min_x);
    Point2d min_y = *std::min_element(points.begin(), points.end(), &_compare_min_y);
    Point2d max_x = *std::max_element(points.begin(), points.end(), &_compare_min_x);
    Point2d max_y = *std::max_element(points.begin(), points.end(), &_compare_min_y);
    
    points.clear();
    points.push_back(min_y);
    points.push_back(min_x);
    points.push_back(max_y);
    points.push_back(max_x);
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
    cout<<"determinant: "<<D<<endl;
    if(!D){
        std::cout<<"Couldn't pass the determinant."<<std::endl;
        return false;
    }

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
    cout<<"Crossing Angle : "<<phi<<endl;
    
    Point3d pc(d *cos(theta0)*sin(phi) / sin(phi), -d *cos(theta0)*cos(phi) + cos(theta1) / sin(phi), d *sin(theta0)*sin(theta1)*sin(rho) / sin(phi));
    cout << "Principle point :\n" << pc << endl;
        
    Point2f vecTranslate;
    double t0, t1, s0, s1;
    Point2f us0, us1;
    us0 = GetIntersectPoint(u0, om, um, wd1);
    us1 = GetIntersectPoint(u1, om, um, wd0);
        
    s0 = GetDistance(us0, um) / l0;
    s1 = GetDistance(us1, um) / l1;
        
    t0 = s0*(l0 + l2) / (s0*l0 + (2 - s0)*l2);
    t1 = s1*(l0 + l3) / (s1*l1 + (2 - s1)*l3);
        
    vecTranslate.x = (t1 - t0)*cos(phi);
    vecTranslate.y = -(t1 + t0)*sin(phi);
    
    trans[0]=pc.x;
    trans[1]=pc.y;
    trans[2]=pc.z;
    
    Vector3d N1(pc.x, pc.y, pc.z);
    Vector3d N2(0, 0, 1);
    N1.normalize(), N2.normalize();
    q = Quaternion<double>::FromTwoVectors(N1, N2);
    //(angle, axis.x, axis.y, axis.z);
    q.normalize();
    
    std::cout << "q: " << q.x() << ", " << q.y() << ", " << q.z() << ", " << q.w() << std::endl;
    std::cout << "vector translate:\n" << vecTranslate.x << ", " << vecTranslate.y << std::endl;
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
/*
    line(out, quadCentered.points[0], quadCentered.points[1], Scalar(0,255,255));
    line(out, quadCentered.points[1], quadCentered.points[2], Scalar(0,255,255));
    line(out, quadCentered.points[2], quadCentered.points[3], Scalar(0,255,255));
    line(out, quadCentered.points[3], quadCentered.points[0], Scalar(0,255,255));
    line(out, quadCentered.points[0], quadCentered.points[2], Scalar(0,255,255));
    line(out, quadCentered.points[1], quadCentered.points[3], Scalar(0,255,255));
    circle(out, quadCentered.points[0], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, quadCentered.points[1], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, quadCentered.points[2], 5, cv::Scalar(0, 255, 255), 4, 5);
    circle(out, quadCentered.points[3], 5, cv::Scalar(0, 255, 255), 4, 5);*/
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


void HT_findSquares( const Mat& image, vector<vector<Point2f> >& squares )
{
    
    Mat res, points, labels, centers;
    int x, y, n, nPoints;
    
    Mat grayImg, binImg, edgeImg;
    cvtColor(image, grayImg, CV_BGR2GRAY);
    threshold(grayImg, binImg, 150, 255, THRESH_OTSU);
    Canny(binImg, edgeImg, 100, 100, 3);
    int w = edgeImg.cols;
    int h = edgeImg.rows;
    int threshold;
    
    //Create the accu
    double hough_h = ((sqrt(2.0) * (double)max(h,w)) / 2.0);
    int _accu_h = hough_h * 2.0; // [-r, +r]
    int _accu_w = 360;
    Mat accu(_accu_h, _accu_w, CV_32SC1);
    unsigned int *_accu = new unsigned int[_accu_h*_accu_w];
    
    double center_x = w/2;
    double center_y = h/2;
    
    
    for(int y=0;y<h;y++)
    {
        for(int x=0;x<w;x++)
        {
            if( edgeImg.data[ (y*w) + x] > 250 )
            {
                for(int t=0;t<360;t++)
                {
                    double r = ( ((double)x - center_x) * cos((double)t * DEG2RAD)) + (((double)y - center_y) * sin((double)t * DEG2RAD));
                    //accu.data[(int)((round(r + hough_h ) * 360)) + t]++;
                    _accu[ (int)((round(r + hough_h ) * 360)) + t]++;
                }
            }
        }
    }
    
    threshold = 10;
    if(threshold == 0)
        threshold = w>h?w/4:h/4;
    unsigned int *thres_accu = new unsigned int[_accu_h*_accu_w];
    cv::Mat img_res = image.clone();
    
    //Search the accumulator
    std::vector< std::pair< std::pair<int, int>, std::pair<int, int> > > lines;
    
    if(_accu == 0){
        cerr<<"the accumulator is null"<<endl;
        return;
    }
    
    for(int r=0;r<_accu_h;r++)
    {
        for(int t=0;t<_accu_w;t++)
        {
            if((int)_accu[(r*_accu_w) + t] >= threshold)
            {
                thres_accu[(r*_accu_w) + t] = _accu[(r*_accu_w) + t];
                //Is this point a local maxima (9x9)
                /*int max = _accu[(r*_accu_w) + t];
                 for(int ly=-4;ly<=4;ly++)
                 {
                 for(int lx=-4;lx<=4;lx++)
                 {
                 if( (ly+r>=0 && ly+r<_accu_h) && (lx+t>=0 && lx+t<_accu_w)  )
                 {
                 if( (int)_accu[( (r+ly)*_accu_w) + (t+lx)] > max )
                 {
                 max = _accu[( (r+ly)*_accu_w) + (t+lx)];
                 ly = lx = 5;
                 }
                 }
                 }
                 }
                 if(max > (int)_accu[(r*_accu_w) + t])
                 continue;
                 */
                
                int x1, y1, x2, y2;
                x1 = y1 = x2 = y2 = 0;
                
                if(t >= 45 && t <= 135)
                {
                    //y = (r - x cos(t)) / sin(t)
                    x1 = 0;
                    y1 = ((double)(r-(_accu_h/2)) - ((x1 - (w/2) ) * cos(t * DEG2RAD))) / sin(t * DEG2RAD) + (h / 2);
                    x2 = w - 0;
                    y2 = ((double)(r-(_accu_h/2)) - ((x2 - (w/2) ) * cos(t * DEG2RAD))) / sin(t * DEG2RAD) + (h / 2);
                }
                else
                {
                    //x = (r - y sin(t)) / cos(t);
                    y1 = 0;
                    x1 = ((double)(r-(_accu_h/2)) - ((y1 - (h/2) ) * sin(t * DEG2RAD))) / cos(t * DEG2RAD) + (w / 2);
                    y2 = h - 0;
                    x2 = ((double)(r-(_accu_h/2)) - ((y2 - (h/2) ) * sin(t * DEG2RAD))) / cos(t * DEG2RAD) + (w / 2);
                }
                
                lines.push_back(std::pair< std::pair<int, int>, std::pair<int, int> >(std::pair<int, int>(x1,y1), std::pair<int, int>(x2,y2)));
                
            }
        }
    }
    
    std::cout << "lines: " << lines.size() << " " << threshold << std::endl;
    
    //Draw the results
    std::vector< std::pair< std::pair<int, int>, std::pair<int, int> > >::iterator it;
    RNG rng(12345);
    for(it=lines.begin();it!=lines.end();it++)
    {
        cv::line(img_res, cv::Point(it->first.first, it->first.second), cv::Point(it->second.first, it->second.second), Scalar(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255)), 1, 8);
    }
    
    //Visualize all
    int maxa = 0;
    for(int p=0;p<(_accu_h*_accu_w);p++)
    {
        if((int)_accu[p] > maxa)
            maxa = thres_accu[p];
    }
    double contrast = 5.0;
    double coef = 255.0 / (double)maxa * contrast;
    
    cv::Mat img_accu(_accu_h, _accu_w, CV_8UC3);
    for(int p=0;p<(_accu_h*_accu_w);p++)
    {
        unsigned char c = (double)_accu[p] * coef < 255.0 ? (double)_accu[p] * coef : 255.0;
        img_accu.data[(p*3)+0] = 255;
        img_accu.data[(p*3)+1] = 255-c;
        img_accu.data[(p*3)+2] = 255;
    }


    cv::imshow("lines", img_res);
    //cv::imshow("edges", edgeImg);
    cv::imshow("accumulator", img_accu);
    //cv::imshow("mat_accu",accu);
}

void ST_findSquares( const Mat& image, vector<vector<Point2f> >& squares )
{
    squares.clear();
    RNG rng(12345);
    Mat grayImg, binImg;
    cvtColor(image, grayImg, CV_BGR2GRAY);
    threshold(grayImg, binImg, thresh, 255, THRESH_OTSU);
    
    /// Parameters for Shi-Tomasi algorithm
    vector<Point2f> corners;
    double qualityLevel = 0.1;
    double minDistance = 1;
    int blockSize = 5;
    bool useHarrisDetector = true;
    double k = 0.04;
    
    /// Copy the source image
    Mat copy;
    copy = image.clone();
    
    /// Apply corner detection
    goodFeaturesToTrack( binImg,
                        corners,
                        100,
                        qualityLevel,
                        minDistance,
                        Mat(),
                        blockSize,
                        useHarrisDetector,
                        k );
    
    
    /// Draw corners detected
    cout<<"** Number of corners detected: "<<corners.size()<<endl;
    int r = 4;
    for( int i = 0; i < corners.size(); i++ )
    { circle( copy, corners[i], r, Scalar(rng.uniform(0,255), rng.uniform(0,255),
                                          rng.uniform(0,255)), -1, 8, 0 ); }
    
    /// Show what you got
    namedWindow( "corner", CV_WINDOW_AUTOSIZE );
    imshow( "corner", copy );
    
    /// Set the neeed parameters to find the refined corners
    Size winSize = Size( 5, 5 );
    Size zeroZone = Size( -1, -1 );
    TermCriteria criteria = TermCriteria( CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 40, 0.001 );
    
    /// Calculate the refined corner locations
    cornerSubPix( binImg, corners, winSize, zeroZone, criteria );
    imshow("corner",copy);
    waitKey(0);
}
static void OpenCV_findSquares( const Mat& image, vector<vector<Point> >& squares ){
    Mat pyr, timg, gray0(image.size(), CV_8U), gray;
    
    // down-scale and upscale the image to filter out the noise
    pyrDown(image, pyr, Size(image.cols/2, image.rows/2));
    pyrUp(pyr, timg, image.size());
    vector<vector<Point> > contours;
    
    // find squares in every color plane of the image
    for( int c = 0; c < 3; c++ )
    {
        int ch[] = {c, 0};
        mixChannels(&timg, 1, &gray0, 1, ch, 1);
        
        // try several threshold levels
        for( int l = 0; l < N; l++ )
        {
            // hack: use Canny instead of zero threshold level.
            // Canny helps to catch squares with gradient shading
            if( l == 0 )
            {
                // apply Canny. Take the upper threshold from slider
                // and set the lower to 0 (which forces edges merging)
                Canny(gray0, gray, 0, thresh, 5);
                // dilate canny output to remove potential
                // holes between edge segments
                dilate(gray, gray, Mat(), Point(-1,-1));
            }
            else
            {
                // apply threshold if l!=0:
                //     tgray(x,y) = gray(x,y) < (l+1)*255/N ? 255 : 0
                gray = gray0 >= (l+1)*255/N;
            }
            
            // find contours and store them all as a list
            findContours(gray, contours, RETR_LIST, CHAIN_APPROX_SIMPLE);
            
            vector<Point> approx;
            
            // test each contour
            for( size_t i = 0; i < contours.size(); i++ )
            {
                // approximate contour with accuracy proportional
                // to the contour perimeter
                approxPolyDP(Mat(contours[i]), approx, arcLength(Mat(contours[i]), true)*0.02, true);
                
                // square contours should have 4 vertices after approximation
                // relatively large area (to filter out noisy contours)
                // and be convex.
                // Note: absolute value of an area is used because
                // area may be positive or negative - in accordance with the
                // contour orientation
                if( approx.size() == 4 &&
                   fabs(contourArea(Mat(approx))) > 1000 &&
                   isContourConvex(Mat(approx)) )
                {
                    double maxCosine = 0;
                    
                    for( int j = 2; j < 5; j++ )
                    {
                        // find the maximum cosine of the angle between joint edges
                        double cosine = fabs(angle(approx[j%4], approx[j-2], approx[j-1]));
                        maxCosine = MAX(maxCosine, cosine);
                    }
                    
                    // if cosines of all angles are small
                    // (all angles are ~90 degree) then write quandrange
                    // vertices to resultant sequence
                    if( maxCosine < 0.5 )
                        squares.push_back(approx);
                }
            }
        }
    }
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
