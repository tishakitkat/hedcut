#include "wcvt.h"
#include <iostream>
#include <unistd.h>
bool compareCell(const std::pair<float, cv::Point>& p1, const std::pair<float, cv::Point>& p2)
{
	if (p1.first == p2.first)
	{
		if (p1.second.x == p2.second.x) return p1.second.y > p2.second.y;
		return p1.second.x > p2.second.x;
	}
	return p1.first > p2.first;
}

 


void make_equal(std::pair<float, cv::Point> p1, std::pair<float, cv::Point> p2){
	p1.first = p2.first;
	p1.second.x = p2.second.x;
	p1.second.y = p2.second.y;
	//std::cout<<"in make equal "<<" first "<<p1.first<<" second "<<p2.first<<std::endl;
	//std::cout<<"in make equal "<<" 1y"<<p1.second.x<<" 2y"<<p2.second.x<<std::endl;
	//std::cout<<"in make equal "<<" 1y"<<p1.second.y<<" 2y"<<p2.second.y<<std::endl;
}



void insertionSort(std::vector< std::pair<float, cv::Point>> point){
	//std::cout<<"insertion sort"<<std::endl;
  int i;
	const std::pair<float, cv::Point> key;
	int j;
	
    for (i = 1; i < point.size(); i++)
    {
        make_equal(key,point[i]);
        j = i-1;
 
        /* Move elements of arr[0..i-1], that are
           greater than key, to one position ahead
           of their current position */
        while (j >= 0 && compareCell(point[j], key))
        {
					//std::cout<<"there are there "<<std::endl;
				make_equal(point[j+1] ,point[j]);
           // point[j+1] = point[j];
            j = j-1;
        }
				make_equal(point[j+1] ,key);
        //point[j+1] = key;
    }
}



//buil the VOR once
void CVT::vor(cv::Mat &  img)
{
	//Generate virtual high resolution image
	cv::Size res(img.size().width * subpixels, img.size().height * subpixels);
	cv::Mat resizedImg(res.width, res.height, CV_LOAD_IMAGE_GRAYSCALE);
	cv::resize(img, resizedImg, res, 0, 0, CV_INTER_LINEAR);

	cv::Mat dist(res, CV_32F, cv::Scalar::all(FLT_MAX));
	cv::Mat root(res, CV_16U, cv::Scalar::all(USHRT_MAX));
	cv::Mat visited(res, CV_8U, cv::Scalar::all(0));
	

	//init
	std::vector< std::pair<float, cv::Point> > open;
	ushort site_id = 0;
	//clock_t startTime, endTime;
	
	for (auto& c : this->cells)
	{
		//std::cout <<"before (x" << c.site.x <<" , y "<<c.site.y<<")"<<" subpixels " <<subpixels<<std::endl;
		//Redirect input site to the position in resized Image
		c.site.x *= subpixels;
		c.site.y *= subpixels;

		if (debug)
		{
			if (c.site.x<0 || c.site.x>res.height)
				std::cout << "! Warning: c.site.x=" << c.site.x << std::endl;

			if (c.site.y<0 || c.site.y>res.width)
				std::cout << "! Warning: c.site.y=" << c.site.y << std::endl;
		}

		
		
		cv::Point pix((int)c.site.x, (int)c.site.y);
		float d = color2dist(resizedImg, pix);
		dist.at<float>(pix.x, pix.y) = d;
		root.at<ushort>(pix.x, pix.y) = site_id++;
		
		
		//startTime = clock();
		if(useqsort){
			//std::cout<<" begin "<<std::endl;
			
			//make for less swapping
			open.push_back(std::make_pair(d, pix));
			if(open.size() > 1){
			
			insertionSort(open);
		}
			//sleep(5);
			
		}
		else{
		open.push_back(std::make_pair(d, pix));
		}
		//std::cout << "the open ("<<d<<","<<" (x" << c.site.x <<" , y "<<c.site.y<<"))"<< std::endl;
		c.coverage.clear();
	}


	std::make_heap(open.begin(), open.end(), compareCell);
	
	
	
	
	//endTime = clock();
	//std::cout<<"total sort time "<<(double)endTime-startTime<<std::endl;
	//sleep(5);
	
	
	
	//propagate
	while (open.empty() == false)
	{
		//
		
		std::pop_heap(open.begin(), open.end(), compareCell);
		auto cell = open.back();
		auto& cpos = cell.second;
		open.pop_back();

		//check if the distance from this cell is already updated
		if (cell.first > dist.at<float>(cpos.x, cpos.y)) continue;
		if (visited.at<uchar>(cpos.x, cpos.y) > 0) continue; //visited
		visited.at<uchar>(cpos.x, cpos.y) = 1;

		//check the neighbors
		for (int dx =-1; dx <= 1; dx++) //x is row
		{
			int x = cpos.x + dx;
			if (x < 0 || x >= res.height) continue;
			for (int dy = -1; dy <= 1; dy++) //y is column
			{
				if (dx == 0 && dy == 0) continue; //itself...

				int y = cpos.y + dy;
				if (y < 0 || y >= res.width) continue;
                cv::Point pt_tmp(x, y);
				float newd = dist.at<float>(cpos.x, cpos.y) + color2dist(resizedImg, pt_tmp);
				float oldd = dist.at<float>(x, y);

				if (newd < oldd)
				{
					if(debug)
					dist.at<float>(x, y)=newd;
					root.at<ushort>(x, y) = root.at<ushort>(cpos.x, cpos.y);
					open.push_back(std::make_pair(newd, cv::Point(x, y)));
					std::push_heap(open.begin(), open.end(), compareCell);
				}
			}//end for dy
		}//end for dx
	}//end while

	//collect cells
	for (int x = 0; x < res.height; x++)
	{
		for (int y = 0; y < res.width; y++)
		{
			ushort rootid = root.at<ushort>(x, y);
			this->cells[rootid].coverage.push_back(cv::Point(x,y));
		}//end y
	}//end x
	
	//remove empty cells...
	int cvt_size = this->cells.size();
	for (int i = 0; i < cvt_size; i++)
	{
		if (this->cells[i].coverage.empty())
		{
			this->cells[i] = this->cells.back();
			this->cells.pop_back();
			i--;
			cvt_size--;
		}
	}//end for i

	if (debug)
	{
		//this shows the progress...
		double min;
		double max;
		cv::minMaxIdx(dist, &min, &max);
		cv::Mat adjMap;
		cv::convertScaleAbs(dist, adjMap, 255 / max);
		//cv::applyColorMap(adjMap, adjMap, cv::COLORMAP_JET);

		for (auto& c : this->cells)
		{
			cv::circle(adjMap, cv::Point( c.site.y, c.site.x ), 2, CV_RGB(0, 0, 255), -1);
		}

		cv::imshow("CVT", adjMap);
		cv::waitKey(5);
	}
}

void CVT::compute_weighted_cvt(cv::Mat &  img, std::vector<cv::Point2d> & sites)
{
	//inint 
	int site_size = sites.size();
	this->cells.resize(site_size);
	for (int i = 0; i < site_size; i++)
	{
		this->cells[i].site = sites[i];
	}

	float max_dist_moved = FLT_MAX;

	int iteration = 0;
	//llyods algorithms
	do
	{
		vor(img); //compute voronoi
		max_dist_moved = move_sites(img);

		//if (debug) std::cout << "[" << iteration << "] max dist moved = " << max_dist_moved << std::endl;
		iteration++;
		//sleep(3);
	} while (max_dist_moved>max_site_displacement && iteration < this->iteration_limit);

	//if (debug) cv::waitKey();
}
