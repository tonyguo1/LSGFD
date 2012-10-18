//Using the data in DataFiles directory (which consists of points randomly distributed in the unit cube [0,1]x[0,1]x[0,1])
// For running the code, ./<exectuable> ./Datafiles/<.mat> <treedepth> 
#include "octree.h"
#include "particles.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <string>

void bruteforcesphere(const double *xp,
                      const double *yp,
                      const double *zp,
		               int  numpoints,
                      const double& search_x,
                      const double& search_y,
                      const double& search_z,
                      const double& radius,
                      std::deque<int>& nbrlist );

void printbruteforcesphere(const double *xp,
                           const double *yp,
                           const double *zp,
                           const double& search_x,
                           const double& search_y,
                           const double& search_z,
                           const double& search_radius,
                           const std::deque<int>& nbrlist) ;

void printoctreesearchsphere(const std::string& str,
                             const double *xp,
                             const double *yp,
                             const double *zp,
                             const double& search_x,
                             const double& search_y,
                             const double& search_z,
                             const double& search_radius,
                             const std::deque<int>& nbrlist);


int main(int argc, char *argv[]) 
{

  clock_t start,stop;
  int treedepth= atoi(argv[2]); // Get tree depth from the command line.
  std::ifstream	readin(argv[1]);// argv[1]=data file containing the coordinates of the particles.
  int           N              ;// Number of particles. Specified on the first line of the data file
  readin>>N                    ;

  //Now we initiate the Particle class we will use throughout the program.
  //Data will be fed into this particle class.
  //A simple particle just has the attribute of position and nothing else. 
 
 
  double *xp = new double[N];
  double *yp = new double[N];
  double *zp = new double[N];

  start=clock();
  for(int i=0 ; i < N; ++i)  
  {
    readin >>  xp[i]  >>  yp[i]  >>  zp[i]    ;
  } 
 
  readin.close();
  stop=clock();
  std::cout << " Time taken for reading the data file is " <<  (stop - start) / (double) CLOCKS_PER_SEC << std::endl;

 

  std::cout << " Number of particles " << N          << std::endl;
  std::cout << " Depth               " << treedepth  << std::endl;
  // Create the particles data structure and feed the information on the  arrays xp, yp and zp into this data structure. 
  // You don't need to delete xp, yp, and zp. The destructor of DefaultCPUParticles will delete them automatically when it perishes.
  DefaultCPUParticles some_particles(xp, yp, zp , N, treedepth);
 
   start = clock();
   some_particles.buildParticleSearchDataStructure();
   stop = clock();
   std::cout << " Time taken for building the Octree is " << (stop - start) / (double) CLOCKS_PER_SEC  << std::endl;


   //**************************************************************************************************** 
   // Print the neighbours of a sample search point to a file, one using brute force and one using octree
   // Done for verification of correctness of octree vs brute force. 
   // Note that the ORDER in which the neighbours are printed is different 
   

  double sample_search_pt_x = 0.840187717155;
  double sample_search_pt_y = 0.482490656656;
  double sample_search_pt_z = 0.979434129307;
  double search_radius      = 0.09;

  std::cout<< "\n\nFor verifying correctness we print the neighbours of " << "[  " << sample_search_pt_x <<"   " << sample_search_pt_y <<"   " << sample_search_pt_z <<"  ]  to 2 files "<< std::endl;
  std::cout<< "Radius chosen for the range search is "                         << search_radius                                                                 <<std::endl;
   // Searching via octree
  std::deque<int> *octree_search_result2 = new std::deque<int>;
  some_particles.searchNeighbor(sample_search_pt_x,
                                sample_search_pt_y,
                                sample_search_pt_z ,
                                search_radius,
                                (void **)(&octree_search_result2) );   

  std::cout<<"\n\nWriting octree search based  neighbours to text file octree_search_neighbours.txt\n\n"<<std::endl;
  printoctreesearchsphere("octree_search_neighbours.txt", xp, yp, zp, sample_search_pt_x, sample_search_pt_y, sample_search_pt_z, search_radius, *octree_search_result2);

  // Searching via brute force
  std::deque<int> nbrlist_bruteforce ;
  bruteforcesphere( xp , 
                    yp, 
                    zp ,  
                    N , 
                    sample_search_pt_x , 
                    sample_search_pt_y ,
                    sample_search_pt_z , 
                    search_radius , 
                    nbrlist_bruteforce);
  std::cout<<"\nWriting brute force search based neighbours to text file brute_force_neighbours.txt"<<std::endl;
  printbruteforcesphere(xp, yp, zp, sample_search_pt_x, sample_search_pt_y, sample_search_pt_z,   search_radius,  nbrlist_bruteforce);
  return 0;
}


















/////////////////////////////////////////////////////////////////////////////////////////////////
// Extra functions for printing out on the screen and files.
/////////////////////////////////////////////////////////////////////////////////////////////////

inline double distsq(const double& xp, 
                     const double& yp, 
                     const double& zp, 
                     const double& search_x,
                     const double& search_y,
                     const double& search_z)
{

  //Return the square of the distance between search point and the considered point.
return (search_x-xp)*(search_x-xp) +	\
       (search_y-yp)*(search_y-yp) +    \
       (search_z-zp)*(search_z-zp) ;
}




void printoctreesearchsphere(const std::string& str,
                             const double *xp,
                             const double *yp,
                             const double *zp,
                             const double& search_x,
                             const double& search_y,
                             const double& search_z,
                             const double& search_radius,
                             const std::deque<int>& nbrlist) 
{


 std::ofstream octreesearch_file;
 octreesearch_file.open(str.c_str());

  octreesearch_file<<"Number of neighbours detected are "<<nbrlist.size()                      <<std::endl;

  octreesearch_file<<"\n\nSearch point is "                  <<search_x<<"\t"
                                                         <<search_y<<"\t"
	                                                 <<search_z                           <<std::endl;
  octreesearch_file<<"Search Radius is "                 <<search_radius                      <<std::endl;


  octreesearch_file<<"\nIn the table below POSN indicates the position of the neighbours in the XYZ arrays used while building the octree.\n"<<std::endl;
  octreesearch_file<<"\n-------------------------------------------------------------------- " << std::endl;
  octreesearch_file<<"              OCTREE SEARCH NEIGHBOURS                     "               <<std::endl;
  octreesearch_file<<"\n------------------------------------------------------------------"    << std::endl;
  octreesearch_file<<"#\tPOSN\tXCD\tYCD\tZCD\t\tDISTANCE TO SEARCHPT"                       << std::endl;
  octreesearch_file<<"-------------------------------------------------------------------- "   << std::endl;      

  for (unsigned int i = 0; i < nbrlist.size(); ++i)
    {
          octreesearch_file << std::setprecision(3)
                << i <<"\t"
                << nbrlist[i]     <<"\t"
                << xp[nbrlist[i]] <<"\t"
                << yp[nbrlist[i]] <<"\t"
		<< zp[nbrlist[i]] <<"\t\t"
                << sqrt( distsq( xp[nbrlist[i]], yp[nbrlist[i]], zp[nbrlist[i]], search_x, search_y, search_z) )<<std::endl;
    }//end for
  octreesearch_file<<"-------------------------------------------------------------------- " << std::endl;

  octreesearch_file.close();

} 



//Section related to the brute-force.

void printbruteforcesphere(const double *xp,
                           const double *yp,
                           const double *zp,
                           const double& search_x,
                           const double& search_y,
                           const double& search_z,
                           const double& search_radius,
                           const std::deque<int>& nbrlist) 

{

  std::ofstream bruteforce_file;
  bruteforce_file.open("bruteforce_neighbours.txt");
  bruteforce_file<<"Number of neighbours detected are "<<nbrlist.size()                      <<std::endl;
 
  bruteforce_file<<"\n\nSearch point is "                  <<search_x<<"\t"
                                                         <<search_y<<"\t"
	                                                 <<search_z                           <<std::endl;
  bruteforce_file<<"Search Radius is "                 <<search_radius                      <<std::endl;


  bruteforce_file<<"\nIn the table below POSN indicates the position of the neighbours in the XYZ arrays used while building the octree.\n"<<std::endl;
  bruteforce_file<<"\n-------------------------------------------------------------------- " << std::endl;
  bruteforce_file<<"              BRUTE FORCE NEIGHBOURS                     "               <<std::endl;
  bruteforce_file<<"\n------------------------------------------------------------------"    << std::endl;
  bruteforce_file<<"NBR#\tPOSN\tXCD\tYCD\tZCD\t\tDISTANCE TO SEARCHPT"                       << std::endl;
  bruteforce_file<<"-------------------------------------------------------------------- "   << std::endl;      

  for (unsigned int i = 0; i < nbrlist.size(); ++i)
    {
          bruteforce_file << std::setprecision(3)
                << i <<"\t"
                << nbrlist[i]     <<"\t"
                << xp[nbrlist[i]] <<"\t"
                << yp[nbrlist[i]] <<"\t"
		<< zp[nbrlist[i]] <<"\t\t"
                << sqrt( distsq( xp[nbrlist[i]], yp[nbrlist[i]], zp[nbrlist[i]], search_x, search_y, search_z) )<<std::endl;
    }//end for
  bruteforce_file<<"-------------------------------------------------------------------- " << std::endl;

  bruteforce_file.close();

}//end function printbruteforcesphere



void bruteforcesphere(const double *xp,
                      const double *yp,
                      const double *zp,
		               int  numpoints,
                      const double& search_x,
                      const double& search_y,
                      const double& search_z,
                      const double& radius,
                      std::deque<int>& nbrlist )
{

for (int i = 0; i <  numpoints; ++i)
  {
    double distance_squared_to_search_pt =  distsq( xp[i] , yp[i], zp[i], search_x, search_y, search_z);
   //Check if point lies within the sphere. Ensure that the point does not count itself as its own neighbour!
    if (    (distance_squared_to_search_pt < radius*radius) && (distance_squared_to_search_pt > 0)   )
      {
     	  	nbrlist.push_back(i);
      } 
    
  }//end for


}//end function bruteforce sphere.

