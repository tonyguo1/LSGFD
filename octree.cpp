/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \author Gaurish Telang <gaurish108@gaurish108> 
 * \date   Mon Jan. 30 2012
 * 
 * \brief Implemntations of the methods and helper functions mentioned in the Octree class.
 */
#include "octree.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <numeric> 
#include <fstream>
#include <vector>
inline bool is_node_intersect_search_region(const double min_x,
                                            const double max_x,
                                            const double min_y,
                                            const double max_y,
                                            const double min_z,
                                            const double max_z,
                                            const double& search_x, 
                                            const double& search_y, 
                                            const double& search_z, const double& radius) 
{
  //The following calculation has been fashioned such that:
  //if final value of squared_dmin == 0 then the point lies inside the node. 
  //if final value of squared_dmin !=0  then the point lies outside the node, AND 
  //             tells the SQUARE of the minimum distance of the point to the points on the node boundary(surface).
  float squared_dmin=0;
  float temp; //Used as a temporary variable to store some intermediate results.
  //Process the x cooridinates
  if( search_x < min_x ) 
    { 
       temp = search_x - min_x; 
       squared_dmin += temp*temp;
    }  
   else if( search_x > max_x )	
    { 
       temp         = search_x - max_x ; 
       squared_dmin += temp*temp                  ;
    }   
 

  //Process the Y-coorindtaes
  if( search_y < min_y )  
  { 
    temp = search_y - min_y ; 
    squared_dmin += temp*temp;
  } 
   else if( search_y > max_y ) 
   { 
    temp = search_y - max_y ; 
    squared_dmin += temp*temp;
   }   



  //Process the Z-coorindtaes
  if( search_z < min_z ) 
  { 
    temp          = search_z - min_z; 
    squared_dmin += temp*temp;
  } 
  
  else if( search_z > max_z ) 
   { 
    temp          = search_z - max_z; 
    squared_dmin += temp*temp;
   }   


  // Based on the calculated value give a  YES/NO answer to the question
  if (squared_dmin <= radius*radius) 
    {
       return true;
    }
  else 
    { 
       return false;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
Octree::Octree(const vector<double> xp,
               const vector<double> yp,
               const vector<double> zp,
               int treedepth, 
               int numOfParticles) :m_vCoordX(xp), m_vCoordY(yp), m_vCoordZ(zp), m_iTotalNumberOfParticles(numOfParticles), m_iMaxDepth(treedepth) {}


int Octree::buildOctree() {

  //Calculate the bounding box limits
   m_dBoundingBox_min_x = *std::min_element(m_vCoordX.begin(), m_vCoordX.end());
   m_dBoundingBox_max_x = *std::max_element(m_vCoordX.begin(), m_vCoordX.end());

   m_dBoundingBox_min_y = *std::min_element(m_vCoordY.begin(), m_vCoordY.end());
   m_dBoundingBox_max_y = *std::max_element(m_vCoordY.begin(), m_vCoordY.end());

   m_dBoundingBox_min_z = *std::min_element(m_vCoordZ.begin(), m_vCoordZ.end());
   m_dBoundingBox_max_z = *std::max_element(m_vCoordZ.begin(), m_vCoordZ.end());

   //Create an array which will hold the morton key and index of each particle.
   //After we fill this array, it will be sorted by the "key" field in the "KeyIndex" struct.
   //The operator < has been overloaded for this purpose
   m_vParticleKeyIndex = new KeyIndex [ m_iTotalNumberOfParticles ];

  for (int i = 0; i < m_iTotalNumberOfParticles; ++i)
    {
      m_vParticleKeyIndex[i].key   = computeKey( m_vCoordX[i], m_vCoordY[i], m_vCoordZ[i] );//function call which computes the morton key from the coordinate arrays.defined below.
      m_vParticleKeyIndex[i].index = i;
    }
  //Now sort the key-index array
  std::sort(m_vParticleKeyIndex, m_vParticleKeyIndex + m_iTotalNumberOfParticles);

  //Now that we have the sorted key-index array we can begin the tree construction level by level
  //For that we need to calculate three helper vectors. bitmasks, baseaddresses and numnodes. 
  //numnodes[i]     is the number of nodes at level i of the octree. zero level correposnds to root nodes. hence numnodes[0] = 1
  //basaddresses[i] is the position in the octree arrays of first node at level i
  //bitmasks[i]     is the Bit-Mask for level i. used in the computation of the numnodes vector 
  std::vector<int> bitmasks(m_iMaxDepth+1);
  std::vector<int> baseaddresses(m_iMaxDepth+1);
  std::vector<int> numnodes(m_iMaxDepth+1);

  //CALCULATE THE BITMASKS VECTOR
   bitmasks[m_iMaxDepth] = ( 1 << (3 * m_iMaxDepth) ) - 1 ;//Deepest level bitmask 
   bitmasks[    0      ] = 0                       ;//Root    level bitmask 
   //Compute remaining level bitmasks
   for (int k = m_iMaxDepth -1 ; k >= 1 ; --k)
    {
      int shiftbits = 3 * ( m_iMaxDepth - k ) ;
      bitmasks[k] = ( bitmasks[ m_iMaxDepth ] >> shiftbits ) << shiftbits;  
    } 

  //COMPUTE THE NUMNODES VECTOR
   numnodes[ 0 ] =  1 ;//Root level has only 1 node i.e. the root node itself.
  for (int k = 1; k < m_iMaxDepth + 1; ++k)
   {
    int count=1;
    
     for ( int i = 1; i < m_iTotalNumberOfParticles ; ++i)
     {
       if ( (m_vParticleKeyIndex[i].key & bitmasks[k]) != (m_vParticleKeyIndex[i-1].key & bitmasks[k] ))
	 {
	   ++count;
	 }
     }
     numnodes[k]=count;   
   }

  //COMPUTE THE BASEADDRESSES VECTOR.
  baseaddresses[0]=0; 
  baseaddresses[1]=1;
  for (unsigned int k = 2; k < baseaddresses.size(); ++k)
   {
     baseaddresses[k]=baseaddresses[k-1]+numnodes[k-1];
   }   

  //With the bitmasks, numnodes and baseaddresses vectors we can now begin the tree computation.
  
  //This is the length of the arrays used in descirbing octrees. 
 m_iTreeLength = std::accumulate(numnodes.begin(),numnodes.end(),0);

  //Now we allocate the arrays using the new operator and m_iTreeLength
  m_vNodeKey			= new int[m_iTreeLength];
  m_vDepth			= new int[m_iTreeLength];
  m_vFirstParticleIndex		= new int[m_iTreeLength];
  m_vNumberOfContainedParticles = new int[m_iTreeLength];
  m_vFirstChildIndex		= new int[m_iTreeLength];
  m_vNumberOfChildren		= new int[m_iTreeLength];  

  m_vLowerLimitOfX		= new double[m_iTreeLength];     
  m_vUpperLimitOfX		= new double[m_iTreeLength]; 

  m_vLowerLimitOfY		= new double[m_iTreeLength]; 
  m_vUpperLimitOfY		= new double[m_iTreeLength]; 

  m_vLowerLimitOfZ		= new double[m_iTreeLength]; 
  m_vUpperLimitOfZ		= new double[m_iTreeLength]; 
  //Start building the last level
  
  //Construct the deepest level of the tree first.
  m_vNodeKey[baseaddresses[m_iMaxDepth]]   = m_vParticleKeyIndex[0].key;
  m_vDepth[baseaddresses[m_iMaxDepth]]     = m_iMaxDepth;

  m_vFirstParticleIndex[baseaddresses[m_iMaxDepth]]	= 0;
 //tree[baseaddressed[m_iMaxDepth]].pnum calculated in the for loop below.

  m_vNumberOfChildren[baseaddresses[m_iMaxDepth]]	= 0;
  m_vFirstChildIndex[baseaddresses[m_iMaxDepth]]	= (-1);
  

  int j   = 1; //This counter changes on encountering a new key.  
  int fix = 0; //Used in the pnum calculation. Records the index where the key change last took place.  
               //pnum = i-fix where i and fix are positions where successive key changes take place.  

   for (int i = 1; i < m_iTotalNumberOfParticles ; ++i)
    {
       if(m_vParticleKeyIndex[i].key!=m_vParticleKeyIndex[i-1].key)
  	{
	     m_vNodeKey[baseaddresses[m_iMaxDepth]+j]				=m_vParticleKeyIndex[i].key;  
     
	     m_vDepth[baseaddresses[m_iMaxDepth]+j]				=m_iMaxDepth; 

	     m_vFirstParticleIndex[baseaddresses[m_iMaxDepth]+j]		=i         ; 
  
       	     m_vNumberOfChildren[baseaddresses[m_iMaxDepth]+j]  		=0         ; //No child nodes. Deepest level.

	     m_vFirstChildIndex[baseaddresses[m_iMaxDepth]+j]   		=(-1);//Special sentinel value. No children since leaf node.    
                                                       
	     m_vNumberOfContainedParticles[baseaddresses[m_iMaxDepth]+j-1]	=i-fix; //Justified above!
             ++j;  

             fix								=i;//change origin of the measure tape. 
  	}//END computation on encountering a new label

    }
   m_vNumberOfContainedParticles[baseaddresses[m_iMaxDepth]+j-1]  = m_iTotalNumberOfParticles - fix;//This cannot be tackled within above for loop.
  //END DEEPEST LEVEL NODE CONSTRUCTION.

  //Now that the deepest level has been constructed lets us start filling in the rest of the octree nodes.
  //Fill the other levels.
   for (int level=m_iMaxDepth-1; level>=0; --level)
    {
      j=1;
      fix=0;    //used in cnum calculation. Records the index where key change LAST took place. 
                //cnum=i-fix. 'i' and 'fix' are positions where successive MASKED key changes take place
            
      //LEVEL CONSTRUCTION BEGINS. WE MUST SCAN LEVEL+1 SECTION TO INSERT APPROPRIATE VALUES IN THE LEVEL SECTION.
      m_vNodeKey[baseaddresses[level]]			= ( m_vNodeKey[baseaddresses[level+1]] & bitmasks[level])     ;
      m_vDepth[baseaddresses[level]]			= level;
      m_vFirstParticleIndex[baseaddresses[level]]	= 0                         ;
      m_vFirstChildIndex[baseaddresses[level]]		= baseaddresses[level+1]     ;
 
      int sum=m_vNumberOfContainedParticles[baseaddresses[level+1]];//This variable used in pnum calculations.
      m_vNumberOfContainedParticles[baseaddresses[level]]	= sum  ; 

      for (int i = 1; i < numnodes[level+1] ; ++i)
  	{
     	  if( ( m_vNodeKey[baseaddresses[level+1]+i] & bitmasks[level] ) != ( m_vNodeKey[baseaddresses[level+1] + i-1] & bitmasks[level] )  )//Relate i and i-1
  	    {
  	       m_vNodeKey[baseaddresses[level]+j]			= ( m_vNodeKey[baseaddresses[level+1]+i] & bitmasks[level] )  ; //Give it the bitmasked key.  
	       m_vDepth[baseaddresses[level]+j]				= level                       ;//Assign the depth
	       m_vFirstParticleIndex[baseaddresses[level]+j]		= m_vFirstParticleIndex[baseaddresses[level+1]+i]                              ;
	       m_vNumberOfContainedParticles[baseaddresses[level]+j-1]	= sum;//Calculate the difference where succesive changes take place in keys masked with a level mask.
	       m_vFirstChildIndex[baseaddresses[level]+j]		= baseaddresses[level+1]+i        ;//Since we encountered a new masked level key. we record the place where this is encountered/
	       m_vNumberOfChildren[baseaddresses[level]+j-1]		= i-fix                           ;//Will have to keep some kind of running sum
	       
  	      ++j  ;
  	      fix=i;
	      sum = m_vNumberOfContainedParticles[baseaddresses[level+1]+i];//initializing sum to pnum of the new key encountered.              
            }

      	  else sum += m_vNumberOfContainedParticles[baseaddresses[level+1]+i];//This is executed only if baseaddresses[level+1]+i and baseaddresses[level+1]+i-1 share the same parent . 
	                                                  //Hence we increment the sum counter which is keeping track of the number of particles within the parent node.
     	}//end inner for 

       m_vNumberOfChildren[baseaddresses[level]+j-1]                = numnodes[level+1]-fix   ;
       m_vNumberOfContainedParticles[baseaddresses[level]+j-1]      = sum                     ;
    }//end outer for

//END LEVEL CONSTRUCTION.

  //Finally we fill in the vertex information required by all the nodes.
  //This is a separate computatin because of its length.
  //Loop over all the tree vector from left to right
  //for center & vertex computation. The only information necessary 
  //for the entire computation here is the DEPTH and the KEY of all the octree nodes.
  //These have already been computed above.

  //Calculate widths of the bounding box along each dimension.
    double xwidth =  m_dBoundingBox_max_x - m_dBoundingBox_min_x   ;
    double ywidth =  m_dBoundingBox_max_y - m_dBoundingBox_min_y   ;
    double zwidth =  m_dBoundingBox_max_z - m_dBoundingBox_min_z  ;

  for (int i = 0; i < m_iTreeLength; ++i)
  {
      //X-y-Z coordinates of the curreent-node center. 
    double center_x; 
    double center_y;
    double center_z;
    //Parse the node-key from left to right by extracting bit-triplets 
    //to calculate the node center. Set the INITIAL values of the center 
    //coordinates to the midpoint coodinates of the bounding box.
    center_x = (  m_dBoundingBox_max_x  + m_dBoundingBox_min_x   ) / 2.0 ;
    center_y = (  m_dBoundingBox_max_y  + m_dBoundingBox_min_y   ) / 2.0 ;
    center_z = (  m_dBoundingBox_max_z  + m_dBoundingBox_min_z   ) / 2.0 ;
    //Following for loop will NOT be executed for the root-node i.e when (i == 0) ! 
    //Will fail the for loop condition right away.( which is good :D )
    for (int j = m_iMaxDepth - 1  ; j >= (m_iMaxDepth - m_vDepth[i]) ; --j)
      {
        double num = (double)(1 << ( m_iMaxDepth - j + 1 ));
        //Perform a right-shift and bit-mask operation (Bitwise AND) with bit-mask 0x7 == (...000111)_2  to extract lowest 3 bits
	unsigned int bit_triplet_extracted = ( m_vNodeKey[i] >> (3*j) ) & 0x7;   
        //Cases below are mutually exclusive.
      	//Case 000
        if (bit_triplet_extracted ==       0)
	  {
            center_x  -=  xwidth/num ;  center_y  -=  ywidth/num  ;  center_z  -=  zwidth/num ; 
    	  }
        //Case 001
        else if (bit_triplet_extracted  == 1)
	  {
	    center_x -=  xwidth/num;    center_y -=  ywidth/num;     center_z +=  zwidth/num; 
	  }
         //Case 010
        else if (bit_triplet_extracted  == 2)
	  {
	    center_x -=  xwidth/num;    center_y +=  ywidth/num;     center_z -=  zwidth/num; 
	  }
         //Case 011
        else if (bit_triplet_extracted  == 3)
	  {
	    center_x -=  xwidth/num;     center_y +=  ywidth/num;     center_z +=  zwidth/num; 
	  }
          //Case 100
        else if (bit_triplet_extracted ==  4)
	  {
	    center_x +=  xwidth/num;     center_y -=  ywidth/num;      center_z -=  zwidth/num; 
	  } 
        //Case 101
        else if (bit_triplet_extracted  == 5)
	  {
	    center_x +=  xwidth/num;      center_y -=  ywidth/num;     center_z +=  zwidth/num; 
	  }
          //Case 110 
        else if (bit_triplet_extracted  == 6)
	  {
	    center_x +=  xwidth/num;       center_y +=  ywidth/num;     center_z -=  zwidth/num; 
	  }
         //Case 111
        else if (bit_triplet_extracted  == 7)
	  {
	    center_x +=  xwidth/num;       center_y +=  ywidth/num;     center_z +=  zwidth/num; 
	  }

      }//END key-parsing for loop. This finishes the computation of the center of node tree[i]
  
    //Calculate vmin and vmax from the center just calculated above.
     m_vLowerLimitOfX[i] = center_x -  ( xwidth / (1 << (m_vDepth[i] + 1 )) ); 
     m_vLowerLimitOfY[i] = center_y -  ( ywidth / (1 << (m_vDepth[i] + 1 )) );
     m_vLowerLimitOfZ[i] = center_z -  ( zwidth / (1 << (m_vDepth[i] + 1 )) ); 

     m_vUpperLimitOfX[i] = center_x + ( xwidth / (1 << (m_vDepth[i]  + 1  )) ); 
     m_vUpperLimitOfY[i] = center_y + ( ywidth / (1 << (m_vDepth[i]  + 1  )) );
     m_vUpperLimitOfZ[i] = center_z + ( zwidth / (1 << (m_vDepth[i]  + 1  )) ); 

  }//END tree-parsing for loop


  //Print the octree information to a file. 

  // std::ofstream octree_file;
  // octree_file.open("printed_octree_vector.txt")			;
  // octree_file << "BOUNDING_BOX INFORMATION :" << std::endl	;
  // octree_file << "Min_x:  "  << m_dBoundingBox_min_x   << std::endl			;
  // octree_file << "Max_x:  "  << m_dBoundingBox_max_x<< std::endl			;

  // octree_file << "\nMin_y:  " << m_dBoundingBox_min_y<< std::endl			;
  // octree_file << "Max_y:  "   << m_dBoundingBox_max_y<< std::endl			;

  // octree_file << "\nMin_z:  "  << m_dBoundingBox_min_z<< std::endl			;
  // octree_file << "Max_z:  "   << m_dBoundingBox_max_z<< std::endl			;
  
  // octree_file << "\nLEVEL\t"<<"BITMASK\t"<<"BASEADDRESSES\t"<<"NUMNODES"<<std::endl;
  // for (int i = 0; i < m_iMaxDepth+1; ++i)
  //   {
  //     octree_file << i << "\t" << bitmasks[i] << "\t" << baseaddresses[i] << "\t\t" << numnodes[i] << std::endl;
  //   }
  //  octree_file << "Length of the octree vector is m_iTreeLength = " << m_iTreeLength << std::endl;  
  // for (int level = 0; level <= m_iMaxDepth; ++level)
  // {
  //   octree_file<<"============================================================================================================================================================"<<std::endl;
  //   octree_file<<"RANK\t"<<"DEPTH\t"<< "KEY\t"<<"PIDX\t"<<"PNUM\t # "<<"CIDX\t"<<"CNUM\t"<<"GLOBAL POSN\t\t"<<"VMIN[0]\tVMIN[1]\tVMIN[2]\t\t\tVMAX[0]\tVMAX[1]\tVMAX[2]"<<std::endl;
  //   octree_file<<"============================================================================================================================================================"<<std::endl; 
  //    for (int i = 0; i < numnodes[level]; ++i)
  // 	{
  // 	  octree_file<<std::setprecision(5)<< i      <<"\t"
  //             << m_vDepth[baseaddresses[level]+i]    <<"\t"  
  // 	      << m_vNodeKey[baseaddresses[level]+i]      <<"\t"
  // 	      << m_vFirstParticleIndex[baseaddresses[level]+i] <<"\t"
  // 	      << m_vNumberOfContainedParticles[baseaddresses[level]+i] <<"\t # "
  // 	      << m_vFirstChildIndex[baseaddresses[level]+i] <<"\t"
  // 	      << m_vNumberOfChildren[baseaddresses[level]+i] <<"\t\t"
  // 		     <<        baseaddresses[level]+i     <<"\t\t" 

  // 		   << m_vLowerLimitOfX[baseaddresses[level]+i] << "\t"
  // 		   << m_vLowerLimitOfY[baseaddresses[level]+i] << "\t"
  // 		   << m_vLowerLimitOfZ[baseaddresses[level]+i] << "\t\t\t"

  // 		   << m_vUpperLimitOfX[baseaddresses[level]+i] << "\t"
  // 		   << m_vUpperLimitOfY[baseaddresses[level]+i] << "\t"
  // 		   << m_vUpperLimitOfZ[baseaddresses[level]+i] << "\t" <<std::endl;

  //      	}
  //     octree_file<<'\n';
  // }



  return 0;

}//End of the buildOctree method
 
////////////////////////////////////////////////////////////////////////////////////////////////////

inline uint32_t Octree::computeKey(const double& x, const double& y, const double& z)
{
  uint32_t    	key=0;//Key which will be returned. 
                      //Excluding this bit, the key can containg upto 30 bits 
                      //corresponding to a max-tree depth of 10. 

  double	left_x=m_dBoundingBox_min_x ,right_x=m_dBoundingBox_max_x ;
  double	left_y=m_dBoundingBox_min_y ,right_y=m_dBoundingBox_max_y ;
  double	left_z=m_dBoundingBox_min_z ,right_z=m_dBoundingBox_max_z ;

  //Midpoint of the various intervals in the following for loop;
  double	midpt_x;
  double        midpt_y;
  double        midpt_z;

   for (int i = 1; i <= m_iMaxDepth; ++i)
    {
      //Compute midpoints.
      midpt_x=(left_x+right_x)/2.0;
      midpt_y=(left_y+right_y)/2.0;
      midpt_z=(left_z+right_z)/2.0;

      //Now we consider 8 cases. EXACTLY ONE of the following will be true within this for loop. 
      //Hence we place it within an if-elsef-if condtruct which guarantees mutual exclusion
      
      //left_x, left_y et al. will get modified within these if-else constructs 
      //since the boundary of the containing interval is continously shrinking.

      //Case 1.
     if (x<midpt_x && \
         y<midpt_y && \
         z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc000 
 	 key=(key << 3);//insert three zeros at the end
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y=left_y;
         right_y=midpt_y;

         left_z=left_z;
         right_z=midpt_z;
       }
   

      //Case 2.
     else if (x>=midpt_x &&   \
              y<midpt_y &&        \
              z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc100 
 	 key=(key << 3) + (1<<2);
 
         left_x=midpt_x;
         right_x=right_x;
       
         left_y=left_y;
         right_y=midpt_y;

         left_z=left_z;
         right_z=midpt_z;
       }

     //Case 3
        else if (x<midpt_x &&   \
                 y>=midpt_y &&        \
                 z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc010 
 	 key=(key << 3) + (1<<1);
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y =midpt_y;
         right_y=right_y;

         left_z=left_z;
         right_z=midpt_z;
       }

     
     //Case 4
      else if ( x<midpt_x &&   \
                y<midpt_y &&        \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc001 
 	 key=(key << 3) + 1;
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y=left_y;
         right_y=midpt_y;

         left_z=midpt_z;
         right_z=right_z;
       }

     //Case 5
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc110 
 	 key=(key << 3) + (1<<2) + (1<<1) ;
 
         left_x=midpt_x;
         right_x=right_x;
       
         left_y=midpt_y;
         right_y=right_y;

         left_z=left_z;
         right_z=midpt_z;
       }

     //Case 6
      else if ( x>=midpt_x &&   \
                y<midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc101 
 	 key=(key << 3) + (1<<2) + 1 ;
 
         left_x =midpt_x;
         right_x=right_x;
       
         left_y =left_y;
         right_y=midpt_y;

         left_z=midpt_z;
         right_z=right_z;
       }

     //Case 7
      else if ( x<midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc011 
 	 key=(key << 3) + (1<<1) + 1 ;
 
         left_x=left_x;
         right_x=midpt_x;
       
         left_y =midpt_y;
         right_y=right_y;

         left_z =midpt_z;
         right_z=right_z;
       }

     //Case 8
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc111 
 	 key=(key << 3) + (1<<2) + (1<<1) +1 ;
 
         left_x=midpt_x;
         right_x=right_x;
       
         left_y =midpt_y;
         right_y=right_y;

         left_z=midpt_z;
         right_z=right_z;
       }

    }//End for loop
  return key ;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
int Octree::searchNeighbor(const double search_x, 
                           const double search_y, 
                           const double search_z, 
                           const double radius, 
			   std::vector<int> *result
			   std::vector<double> *distance) 
{
   std::deque<int> nodelist;
   nodelist.push_back(0); // Push back the index of the root node, which is at position 0 in the octree vector.

   // In this while loop, we find all the leaf nodes of the octree which intersect the search region.  
   // If, at the end of the while loop, nodelist queue is empty, that means no leaf node intersects the search region and hence the search point has no particles in it neighbourhood.
   // In this we don't check whether an octree node is a subset of the search region. 
   // Paradoxically, this seems to be performing better than an algorithm which does. This needs some more investigating. 
   // Probably for deeper octrees and smaller search radii. that algorithm will endup outperforming this one. Also that means the subset check would need to be improved upon, since it could be the bottleneck. 
   while( 1 ) 
   {

     //std::cout << nodelist.size() << std::endl;
     // If the node-list is empty stop the execution of he current function, and return 0, since 0 neighbours have been founf. 
     // Also as a result, you see that nothing is in the deque pointed to be result.
    if( nodelist.empty() ) 
      { 
        return 0;
      }
    // Get the node number at the front of the queue.
    // The popping will be done at the end of the while loop.
    int inode       = nodelist[0];
    int first_index = m_vFirstChildIndex[inode];

    // If the head node of the deque is a leaf that means all the remaining nodes in the deque are also leaves. 
    // So we can now exit the while loop and process look through all the particles within the leaf nodes. 
    if( first_index == -1 ) 
      {
        break;
      }

    int last_index = first_index + m_vNumberOfChildren[ inode ];
    //  Find all the child nodes of the CURRENT node which intersect the sphere having radius "radius".
    //  After finding push_back their index onto the stack.
    for(int i=first_index ; i < last_index ; i++) 
     {
       if(   is_node_intersect_search_region( m_vLowerLimitOfX[i],\
                                              m_vUpperLimitOfX[i],\
                                              m_vLowerLimitOfY[i],\
                                              m_vUpperLimitOfY[i],\
                                              m_vLowerLimitOfZ[i],\
                                              m_vUpperLimitOfZ[i],\
                                              search_x,\
                                              search_y,\
                                              search_z, radius) )
	 {
             nodelist.push_back( i );
	 }
    }
    nodelist.pop_front( );
  }
  // When the  above while-loop ends normally, the queue consists only of leaf nodes intersecting the search region.
  // We now look at every particle in each leaf node looking for neighbouring particles. Probably as an optimization step,
  // it would be interesting to first mark all the leaf nodes which are a subset of the search region , so all the particles within
  // such leaf nodes can be directly added. 
  while( !nodelist.empty() ) 
   {
    
     // Number at the front of the queue
    int index = nodelist[0];
     // Positions of the first and the last particles within the node pointed to by index. 
    int first = m_vFirstParticleIndex[index];
    int end   = first + m_vNumberOfContainedParticles[index];
    
    // Cycle through all the particles within the current leaf node and look for neighbours.
    for (int i = first; i < end; i++) 
     {
      double dx = search_x - m_vCoordX[ m_vParticleKeyIndex[i].index ];
      double dy = search_y - m_vCoordY[ m_vParticleKeyIndex[i].index ];
      double dz = search_z - m_vCoordZ[ m_vParticleKeyIndex[i].index ];
      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( squared_distance < radius*radius) 
        {
          result->push_back(m_vParticleKeyIndex[i].index);
	  distance->push_back(sqrt(squared_distance));
        }
    }
    // Since all the information needed from the queue head has been processes we have no need for it any more. So pop it.
    nodelist.pop_front();
  }
    // Return the number of neighbours found.
    // The actual neighbours are pointed to by the result pointer. 
    return result->size();
}                
