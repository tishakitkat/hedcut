#include "hedcut.h"
#include <time.h>



UnclusterGraph::UnclusterGraph(){
    int cluster_value = 250;
    //small distance you dont want too large a distance
    float cluster_dist = 0.0015f;
    
}


int UnclusterGraph::addEdgeToGraph(cv::Point2d pt, cv::Mat & img){
    Vertex vx = Vertex(pt);
    Vertex *v = &vx;
    this->vertices.push_back(v);
    int size = this->vertices.size();
    std::vector<Vertex*> reverstvertices;
    
  for (int i=0; i<size-1; i++)
         reverstvertices.push_back(this->vertices[i]);
    
    //cout<<"size is a "<<vertices.size()<<endl;
    if(this->vertices.empty())
        return 0;
    
    
    for(int i=0;i<this->vertices.size()-1;i++){
        //vertices.at(i)->pt.x

            float xback = this->vertices.at(i)->pt.x/img.size().width - cluster_dist;
            float xfor = this->vertices.at(i)->pt.x/img.size().width + cluster_dist;
            float ytop =this->vertices.at(i)->pt.y/img.size().height + cluster_dist;
            float ybott = this->vertices.at(i)->pt.y/img.size().height - cluster_dist;
            if(xback > pt.x/img.size().width || xfor > pt.x/img.size().width || ytop > pt.y/img.size().height || ybott > pt.y/img.size().height){
                //cout<<"clustered edges"<<vertices.at(i)->pt.x/img.size().width<<" > "<<cluster_value<<endl;               sleep(2);
                this->vertices.at(i)->pt.x /=img.size().width;
                this->vertices.at(i)->pt.y/=img.size().height;
                
                vx.edges.push_back(this->vertices.at(i)->pt);
                //vx.edges.push_back(this->vertices.at(i)->pt);
                this->vertices.at(i)->edges.push_back(pt);
                if(this->vertices.at(i)->edges.size() > cluster_value || vx.edges.size() > cluster_value ){
                    //this->vertices = reverstvertices;
                  for (int i=0; i<size-1; i++)
                         this->vertices.push_back(reverstvertices[i]);
                    return 1;
                }
                    
            
            }
            else{
                continue;
            }
        }
        

    return 0;   
            
}



//hedcut



Hedcut::Hedcut()
{
	//control flags
	disk_size = 1;        //if uniform_disk_size is true, all disks have radius=disk_size,
	                      //othewise, the largest disks will have their radii=disk_size 

	uniform_disk_size = false; //true if all disks have the same size. disk_size is used when uniform_disk_size is true.
	black_disk = false;        //true if all disks are black ONLY

	//cvt control flags
	cvt_iteration_limit = 100; //max number of iterations when building cvf
	max_site_displacement = 1.01f; //max tolerable site displacement in each iteration. 
	average_termination = false;
	gpu = false;
	subpixels = 1;

	debug = false;

    //for the graph
    bool usegraph = false;
    int cluster_value = 0;
    float cluster_dist = 0.0f;
    //for time
    int timeout = 0;
    //for q sort
    bool useqsort=false;


}



bool Hedcut::build(cv::Mat & input_image, int n)
{
	cv::Mat grayscale;
	cv::cvtColor(input_image, grayscale, CV_BGR2GRAY);

	//sample n points
	std::vector<cv::Point2d> pts;
	sample_initial_points(grayscale, n, pts);

	//initialize cvt
	CVT cvt;
	
	cvt.iteration_limit = this->cvt_iteration_limit;
	cvt.max_site_displacement = this->max_site_displacement;
	cvt.average_termination = this->average_termination;
	cvt.gpu = this->gpu;
	cvt.subpixels = this->subpixels;
	cvt.debug = this->debug;
    cvt.useqsort = this->useqsort;
   cvt.timeout = this->timeout;

	clock_t startTime, endTime;
	startTime = clock();
	
	//compute weighted centroidal voronoi tessellation
	if (cvt.gpu)
		cvt.compute_weighted_cvt_GPU(input_image, pts);
	else
		cvt.compute_weighted_cvt(grayscale, pts);	//*****

	endTime = clock();
	std::cout << "Total time: "<< ((double)(endTime - startTime)) / CLOCKS_PER_SEC << std::endl;

	if (debug) cv::waitKey();

	//create disks
	create_disks(input_image, cvt);

	return true;
}


void Hedcut::sample_initial_points(cv::Mat & img, int n, std::vector<cv::Point2d> & pts)
{


    UnclusterGraph g;
    if(this->cluster_value || this->cluster_dist > 0.0){
        this->usegraph = true;
        g.usegraph = this->usegraph;
        g.cluster_value = this->cluster_value > 0 ? this->cluster_value:g.cluster_value;
        g.cluster_dist = this->cluster_dist > 0.0 ? this->cluster_dist:g.cluster_dist;


    }


	//create n points that spread evenly that are in areas of black points...
	int count = 0;

	cv::RNG rng_uniform(time(NULL));
	cv::RNG rng_gaussian(time(NULL));
	cv::Mat visited(img.size(), CV_8U, cv::Scalar::all(0)); //all unvisited

	while (count < n)
	{
		//generate a random point
		int c = (int)floor(img.size().width*rng_uniform.uniform(0.f, 1.f));
		int r = (int)floor(img.size().height*rng_uniform.uniform(0.f, 1.f));

		//decide to keep basic on a probability (black has higher probability)
		float value = img.at<uchar>(r, c)*1.0/255; //black:0, white:1
		float gr = fabs(rng_gaussian.gaussian(0.8));
		if ( value < gr && visited.at<uchar>(r, c) ==0) //keep
		{
            if(g.usegraph){
                int check = g.addEdgeToGraph(cv::Point(r, c),img);
                    if(check == 0){
                
                    count++;
                    pts.push_back(cv::Point(r, c));
                    visited.at<uchar>(r,c)=1;

                    }
                    else{
                        continue;   
                    }
                    
                }
            }

            else{
			count++;
			pts.push_back(cv::Point(r, c));
			visited.at<uchar>(r,c)=1;
            }
	}

	if (debug)
	{
		cv::Mat tmp = img.clone();
		for (auto& c : pts)
		{
			cv::circle(tmp, cv::Point(c.y, c.x), 2, CV_RGB(0, 0, 255), -1);
		}
		cv::imshow("samples", tmp);
		cv::waitKey();
	}
}

void Hedcut::create_disks(cv::Mat & img, CVT & cvt)
{
	cv::Mat grayscale;
	cv::cvtColor(img, grayscale, CV_BGR2GRAY);

	disks.clear();

	//create disks from cvt
	for (auto& cell : cvt.getCells())
	{
		//compute avg intensity
		unsigned int total = 0;
		unsigned int r = 0, g = 0, b = 0;
		for (auto & resizedPix : cell.coverage)
		{
			cv::Point pix(resizedPix.x / subpixels, resizedPix.y / subpixels);
			total += grayscale.at<uchar>(pix.x, pix.y);
			r += img.at<cv::Vec3b>(pix.x, pix.y)[2];
			g += img.at<cv::Vec3b>(pix.x, pix.y)[1];
			b += img.at<cv::Vec3b>(pix.x, pix.y)[0];
		}
		float avg_v = floor(total * 1.0f/ cell.coverage.size());
		r = floor(r / cell.coverage.size());
		g = floor(g / cell.coverage.size());
		b = floor(b / cell.coverage.size());

		//create a disk
		HedcutDisk disk;
		disk.center.x = cell.site.y; //x = col
		disk.center.y = cell.site.x; //y = row
		disk.color = (black_disk) ? cv::Scalar::all(0) : cv::Scalar(r, g, b, 0.0);
		disk.radius = (uniform_disk_size) ? disk_size : (100 * disk_size / (avg_v + 100));

		//remember
		this->disks.push_back(disk);

	}//end for cell

	//done
}