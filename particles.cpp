#include "particles.h"

namespace std{
////////////////////////////////////////////////////////////////////////////////////////////////////
DefaultParticles::DefaultParticles(vector<double> xp,
                                   vector<double> yp,
                                   vector<double> zp ,
                                   intnum_particles)
{
  m_vCoordX = xp;
  m_vCoordY = yp;
  m_vCoordZ = zp;
  m_iNumOfParticles = num_particles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
DefaultParticles::~DefaultParticles() 
{
  delete m_vCoordX;
  delete m_vCoordY;
  delete m_vCoordZ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
intDefaultParticles::buildParticleSearchDataStructure() {return 0;}
////////////////////////////////////////////////////////////////////////////////////////////////////
intDefaultParticles::searchNeighbor(const double x,
                                     const double y, 
                                     const double z, 
                                     const double radius, 
                                     void **result)  {return 0;}

intDefaultParticles::saveParticleData(ofstream *o) {return 0;}
////////////////////////////////////////////////////////////////////////////////////////////////////
DefaultCPUParticles::DefaultCPUParticles(vector<double> xp,
                                         vector<double> yp,
                                         vector<double> zp , intnum_particles,
                                         inttreedepth)  : DefaultParticles(xp, yp, zp,num_particles), m_Octree(NULL), m_iTreeDepth(treedepth) {}

DefaultCPUParticles::~DefaultCPUParticles() 
{
  if(m_Octree != NULL) 
  {
    delete m_Octree;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////
intDefaultCPUParticles::buildParticleSearchDataStructure()
{
  //Invoke the OCtree building routine
  if(m_Octree != NULL) 
  {
    delete m_Octree;
  }
  m_Octree = new Octree(m_vCoordX, m_vCoordY, m_vCoordZ, m_iTreeDepth, m_iNumOfParticles);
  return m_Octree->buildOctree();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
intDefaultCPUParticles::searchNeighbor(const double x,
                                        const double y, 
                                        const double z, 
                                        const double radius, void **result, void **distance) 
{
  return m_Octree->searchNeighbor( x, y, z, radius, *((vector<int>**)(result)), *((vector<double>**)(distance));
}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
