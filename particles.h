/*!
 * \brief SPH Particle Data Structure. 
 */

#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#include "octree.h"
#include <iostream>

namespace std{
/*!
\brief Partilces interface. Every particles has to implement this interface. DefaultParticles is a implementation class. You can inherit from DefaultParticles.
\see DefaultParticles
*/
class Particles {
public:
  virtual void setLocation(int index, double x, double y, double z) = 0;
  virtual int getLocation(int index, vector<double> location) = 0;
  virtual void setVelocity(int index, double x, double y, double z) = 0;
  virtual int getVelocity(int index, vector<double> velocity) = 0;
  virtual void setDensity(int index, double value) = 0;
  virtual double getDensity(int index) = 0;
  virtual void setRadius(int index, double value) = 0;
  virtual double getRadius(int index) = 0;
  virtual void setMass(int index, double value) = 0;
  virtual double getMass(int index) = 0;
  virtual void setType(int index, int value) = 0;
  virtual int getType(int index) = 0;

  virtual int getNumberOfParticle() = 0;
  virtual int buildParticleSearchDataStructure() = 0;
  virtual int searchNeighbor(const double x, const double y, const double z, const double radius, void **result) = 0;
  virtual int saveParticleData(ofstream *o) = 0;
};


/*!
\brief A default implementation of Particles interface. If you need Particles class for GPU, it it better inheriting from DefaultGPUParticles. If you need Particles class for CPU, you can inherit form DefaultCPUParticles.
\see DefaultCPUParticles DefaultGPUParticles
*/
class DefaultParticles : public Particles {
  
protected:
  /*!
  \brief The total number of particles. This number must be coincide with the length of array member fields.
  */
  int m_iNumOfParticles;

  /*!
  \brief array of x-coordinates of particles.
  */
  vector<double> m_vCoordX;
  /*!
  \brief array of y-coordinates of particles.
  */
  vector<double> m_vCoordY;
  /*!
  \brief array of z-coordinates of particles.
  */
  vector<double> m_vCoordZ;
  /*!
  \brief array of x-direction velocity of particles.
  */
  vector<double> m_vVelocityX;
  /*!
  \brief array of y-direction velocity of particles.
  */
  vector<double> m_vVelocityY;
  /*!
  \brief array of z-direction velocity of particles.
  */
  vector<double> m_vVelocityZ;
  /*!
  \brief array of density of particles.
  */
  vector<double> m_vDensity;
  /*!
  \brief array of radius of particles.
  */
  vector<double> m_vRadius;
  /*!
  \brief array of mass velocity of particles.
  */
  vector<double> m_vMass;
  /*!
  \brief array of type of particles. this field is reserved for particles have different properties.
  */
  vector<int> m_vType;

public:
  DefaultParticles(vector<double> xp, vector<double> yp, vector<double> zp , int num_particles);
  ~DefaultParticles();

public:
  void setLocation(int index, double x, double y, double z) {m_vCoordX[index] = x;
                                                             m_vCoordY[index] = y; 
                                                             m_vCoordZ[index] = z;}

  int getLocation(int index, vector<double> location) {location[0] = m_vCoordX[index];
                                                   location[1] = m_vCoordY[index]; 
                                                   location[2] = m_vCoordZ[index]; return 3;}
 
  void setVelocity(int index, double x, double y, double z) {m_vVelocityX[index] = x;
                                                             m_vVelocityY[index] = y; 
                                                             m_vVelocityZ[index] = z;}


  int getVelocity(int index, vector<double> velocity) {velocity[0] = m_vVelocityX[index];
                                                   velocity[1] = m_vVelocityY[index]; 
                                                   velocity[2] = m_vVelocityZ[index]; return 3;}


  void   setDensity(int index, double value) {  m_vDensity[index] = value;  }
  double getDensity(int index) {  return m_vDensity[index];  }


  void   setRadius(int index, double value) {m_vRadius[index] = value;}
  double getRadius(int index) {return m_vRadius[index];}


  void   setMass(int index, double value) {m_vMass[index] = value;}
  double getMass(int index) {return m_vMass[index];}


  int getNumberOfParticle() {return m_iNumOfParticles;}

  void setType(int index, int value) {m_vType[index] = value;}
  int  getType(int index) {return m_vType[index];}

  virtual int buildParticleSearchDataStructure();

  virtual int searchNeighbor(const double x,
                             const double y, 
                             const double z, 
                             const double radius, void **result);
  int saveParticleData(ofstream *o);
};



/*!
\brief A default implementation of Particles interface for CPU code.
*/
class DefaultCPUParticles: public DefaultParticles 
{
    protected:
       Octree *m_Octree;
       int m_iTreeDepth;
    public:
       DefaultCPUParticles(vector<double> xp,
                           vector<double> yp,
                           vector<double> zp ,
                           int num_particles,
                           int treedepth);

       ~DefaultCPUParticles();
       //void setTreeDepth(int depth);
       virtual int buildParticleSearchDataStructure();
       virtual int searchNeighbor(const double x,
                                  const double y, 
                                  const double z, 
                                  const double radius, 
                                  void **result
                                  void **distance);
};
}

#endif // __PARTICLES_H__
