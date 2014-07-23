// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <cstdlib>
#include <math.h>
#include <map>
#include <stdio.h>
#include <iostream>
#undef max
#undef min
#include <fstream>

//#define _MAC_GRID_DEBUG


// Globals:
MACGrid target;

//Added by Hao:
int runTimes = 200;
bool firstRun = true;
const double vf = 0.1;
const double vd = 0.2;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 


#define FOR_EACH_FACE_X \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X] + 1; i++) 

#define FOR_EACH_FACE_Y \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y] + 1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_FACE_Z \
   for(int k = 0; k < theDim[MACGrid::Z] + 1; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 



MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();		//the default value for each unit is 0;
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);


/* Hao's Data Structure init (BEGIN) */

/*
   mFilter.initialize();
   mDensityBlur.initialize();
   mTargetDensity.initialize();
   mTargetDensityBlur.initialize();
   mDrivingForce.initialize();

   mFilterWidth = theDim[MACGrid::X];
   FilterGenGaussian2D(mFilter, mFilterWidth, 0.7);
   cout << "********** FILTER info *********" << endl;
   for (int i = 0; i < mFilterWidth; i++) {
	   for (int j = 0; j < mFilterWidth; j++) {
		   //cout << "i: " << i << " j: " << j << "  " << mFilter(i, j, 0) << endl;
		   cout << mFilter(i, j, 0) << "   ";
		   if (j == mFilterWidth - 1)
			   cout << endl;
	   }
   }
*/

/* Hao's Data Structure init (END) */

   std::cout << "get called 1" << endl;
   setUpAMatrix();
}

void MACGrid::initialize()
{
   reset();
   std::cout << "Z max: " <<  theDim[MACGrid::Z] << std::endl;
   std::cout << "X max: " <<  theDim[MACGrid::X] << std::endl;
   std::cout << "Y max: " <<  theDim[MACGrid::Y] << std::endl;
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, and velocity.
	//target.mV(0, 1, 0) = 1;
	//target.mV(1, 1, 0) = 1;
	
	/*
	if (firstRun) {
		mTargetDensity(5, 5, 0) = 10;
		mTargetDensity(5, 4, 0) = 10;
		mTargetDensity(5, 3, 0) = 10;
		mTargetDensity(5, 2, 0) = 10;
		mTargetDensity(5, 1, 0) = 10;

		mTargetDensity(4, 3, 0) = 10;
		mTargetDensity(3, 3, 0) = 10;
		mTargetDensity(2, 3, 0) = 10;
		

		mTargetDensity(1, 5, 0) = 10;
		mTargetDensity(1, 4, 0) = 10;
		mTargetDensity(1, 3, 0) = 10;
		mTargetDensity(1, 2, 0) = 10;
		mTargetDensity(1, 1, 0) = 10;
		firstRun = false;
		Convolve(mTargetDensityBlur, mTargetDensity, mFilter, mFilterWidth);
	}
	*/
	if (runTimes > 0) {
		mV(5, 2, 0) = 4;
		mV(5, 3, 0) = 4;
		mV(5, 4, 0) = 4;
		mV(5, 5, 0) = 4;
		mV(10, 6, 0) = 4;
		mV(10, 7, 0) = 4;
		mV(10, 8, 0) = 4;
		mV(10, 9, 0) = 4;
		mV(10, 10, 0) = 4;
		mV(10, 11, 0) = 4;
		mV(10, 12, 0) = 4;
		mV(10, 13, 0) = 4;
		mV(10, 14, 0) = 4;
		mV(10, 15, 0) = 4;
		mV(10, 16, 0) = 4;
		mV(10, 17, 0) = 4;
		mV(10, 18, 0) = 4;
		mV(10, 19, 0) = 4;
		//mV(10, 20, 0) = 4;
		//mV(10, 21, 0) = 4;
		//mV(10, 22, 0) = 4;
		//mV(10, 23, 0) = 4;
		//mU(1, 0, 0) = 1;
		
		
		//mD(5, 0, 0) = 2;
		//mV(5, 1, 0) = 1;
		--runTimes;
		mD(10, 0, 0) = 2;
		mV(10, 1, 0) = 10;
	}
}


void MACGrid::FilterGenGaussian2D(GridData &filter, int halfCellNum, float sigma)
{
	const uint w_2 = halfCellNum >> 1;
	const float fw_2 = float(halfCellNum) / 2.0f * theCellSize;
		
	if (sigma == 0.0f) {
		sigma = w_2 / 2.0f;
	}
	const float c1 = 1.0f / (sqrt(2.0 * M_PI) * sigma);
	const float c2 = 2.0 * sigma * sigma;

	for (int gi = 0; gi < halfCellNum; gi++) {
		for (int gj = 0; gj < halfCellNum; gj++) {
			const float cx = ((float)gi * theCellSize + 0.5 * theCellSize) - fw_2;	//position in the gaussian
			const float cy = ((float)gj * theCellSize + 0.5 * theCellSize) - fw_2;

			filter(gi, gj, 0) = c1 * c1 * exp(-cx * cx / c2) * exp(-cy * cy / c2);
		}
	}
}

void MACGrid::Convolve(GridData &out, const GridData &in, GridData &filter, int filterWidth, int absVal)
{
	cout << "Convolve once" << endl;
	// find center position of kernel (half of kernel size)
	int kCenterX = filterWidth / 2;
	int kCenterY = filterWidth / 2;
	int rows = theDim[MACGrid::X];
	int cols = theDim[MACGrid::Y];
	int kRows, kCols;
	kRows = kCols = filterWidth;

	for(int i = 0; i < rows; ++i)              // rows
	{
		for(int j = 0; j < cols; ++j)          // columns
		{
			int sum = 0;                     // init to 0 before sum

			for(int m = 0; m < kRows; ++m)     // kernel rows
			{
				int mm = kRows - 1 - m;      // row index of flipped kernel

				for(int n = 0; n < kCols; ++n) // kernel columns
				{
					int nn = kCols - 1 - n;  // column index of flipped kernel

					// index of input signal, used for checking boundary
					int ii = i + (m - kCenterY);
					int jj = j + (n - kCenterX);

					// ignore input samples which are out of bound
					if( ii >= 0 && ii < rows && jj >= 0 && jj < cols )
					out(i, j, 0) += in(ii, jj, 0) * filter(mm, nn, 0);
				}
			}
		}
	}
/*
	int ii, jj, kk;
	int lowOrder = (n - 1) / 2;
	int highOrder = n / 2;
	int minx, miny, maxx, maxy;
	int xndx, yndx, ndx;
	int width = theDim[MACGrid::X];
	int height = theDim[MACGrid::Y];

	float *kptr;
	float div;
	float val;
	
	for (int ii = 0; ii < height; i++) {
		miny = ii - lowOrder;
		maxy = ii + highOrder;
		for (int jj = 0; jj < width; jj++) {
			minx = jj - lowOrder;
			maxx = jj + highOrder;
			kptr = filter;
			div = 1.0;
			val = 0.0;
			for (kk = miny; kk <= maxy; k++) {
				yndx = kk * width;
				for (xndx = minx; xndx <= maxx; xndx++) {
					if (kk < 0 || kk >= height || xndx < 0 || xndx >= width)
					{
						div -= *kptr;		//modify
					} else {
						ndx = xndx + yndx;
						val +=
					}
				}
			}
		}
	}
*/
}

void MACGrid::advectVX(const vec3 &index, double dt)
{
	int i = index.n[VX]; 
	int j = index.n[VY];
	int k = index.n[VZ];

	vec3 pt(i * theCellSize,
			j * theCellSize + 0.5 * theCellSize,
			k * theCellSize + 0.5 * theCellSize);
	
	vec3 vel = getVelocity(pt);

	pt.n[VX] -= vel.n[VX] * dt;
	pt.n[VY] -= vel.n[VY] * dt;
	pt.n[VZ] -= vel.n[VZ] * dt;
	target.mU(i, j, k) = getVelocityX(pt);
}

void MACGrid::advectVY(const vec3 &index, double dt)
{
	int i = index.n[VX]; 
	int j = index.n[VY];
	int k = index.n[VZ];

	vec3 pt(i * theCellSize + 0.5 * theCellSize,
			j * theCellSize,
			k * theCellSize + 0.5 * theCellSize);
	
	vec3 vel = getVelocity(pt);

	pt.n[VX] -= vel.n[VX] * dt;
	pt.n[VY] -= vel.n[VY] * dt;
	pt.n[VZ] -= vel.n[VZ] * dt;
	target.mV(i, j, k) = getVelocityY(pt);
}

void MACGrid::advectVZ(const vec3 &index, double dt)
{
	int i = index.n[VX]; 
	int j = index.n[VY];
	int k = index.n[VZ];

	vec3 pt(i * theCellSize + 0.5 * theCellSize,
			j * theCellSize + 0.5 * theCellSize,
			k * theCellSize);

	vec3 vel = getVelocity(pt);

	pt.n[VX] -= vel.n[VX] * dt;
	pt.n[VY] -= vel.n[VY] * dt;
	pt.n[VZ] -= vel.n[VZ] * dt;
	target.mW(i, j, k) = getVelocityZ(pt);
}


void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
/*	
	vec3 vel1, vel2, pt;
	for (int k = 0; k < theDim[MACGrid::Z]; k++) {
		for (int i = 0; i < theDim[MACGrid::X]; i++) {
			for (int j = 0; j < theDim[MACGrid::Y]; j++) {
				if (i == 0)
					continue;
				pt.n[VX] = i * theCellSize;
				pt.n[VY] = j + 0.5 * theCellSize;
				pt.n[VZ] = k + 0.5 * theCellSize;
				
				vel1 = getVelocity(pt);

				pt.n[VX] -= vel1.n[VX] * dt;
				pt.n[VY] -= vel1.n[VY] * dt;
				pt.n[VZ] -= vel1.n[VZ] * dt;

				target.mU(i, j, k) = getVelocityX(pt);
			}
		}
	}
*/

/********************* PRINT(BEGIN) ***************************/
#ifdef _MAC_GRID_DEBUG
	cout << "********** BEFORE ************" << endl;
	

	FOR_EACH_FACE_X {
		cout << "mU i: " << i << " j: " << j << " k: " << k  \
				<< " --> " << this->mU(i, j, k) << endl;
	}
	
	cout << endl;
	FOR_EACH_FACE_Y {
		cout << "mV i: " << i << " j: " << j << " k: " << k  \
				<< " --> " << this->mV(i, j, k) << endl;
	}

	cout << endl;
	FOR_EACH_FACE_Z {
		cout << "mW i: " << i << " j: " << j << " k: " << k  \
				<< " --> " << this->mW(i, j, k) << endl;
	}
#endif

/********************* PRINT(END) ***************************/
	
	FOR_EACH_CELL {
		vec3 index(i, j, k);

		if (i == 0 && j == 0 && k == 0) {
			continue;
		} else if (i == 0 && j == 0 && k != 0) { 		//only update the face on Z direction
			advectVZ(index, dt);
		} else if (i == 0 && j != 0 && k == 0) {	//only update the face on Y direction 
			advectVY(index, dt);
		} else if (i != 0 && j == 0 && k == 0) {
			advectVX(index, dt);
		} else if (i == 0 && j != 0 && k != 0) {
			advectVY(index, dt);
			advectVZ(index, dt);
		} else if (i != 0 && j == 0 && k != 0) {
			advectVX(index, dt);
			advectVZ(index, dt);
		} else if (i != 0 && j != 0 && k == 0) {
			advectVX(index, dt);
			advectVY(index, dt);
		} else if (i != 0 && j != 0 && k != 0) {
			advectVX(index, dt);
			advectVY(index, dt);
			advectVZ(index, dt);
		}
	}

//
//	for (int k = 0; k < theDim[MACGrid::Z]; k++)
//		for (int j = 0; j < theDim[MACGrid::Y] + 1; j++)
//			for (int i = 0; i < theDim[MACGrid::X]; i++) {
//				cout << "i: " << i << " j: " << j << " k: " << k  \
//						<< " --> " << target.mV(i, j, k) << endl;
//			}

	
	//target.mU = mU;
    //target.mV = mV;
    //target.mW = mW;


    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;

/********************* PRINT(BEGIN) ***************************/
#ifdef _MAC_GRID_DEBUG
	cout << "********** AFTER ************" << endl;
	

	FOR_EACH_FACE_X {
		cout << "mU i: " << i << " j: " << j << " k: " << k  \
				<< " --> " << this->mU(i, j, k) << endl;
	}
	
	cout << endl;
	FOR_EACH_FACE_Y {
		cout << "mV i: " << i << " j: " << j << " k: " << k  \
				<< " --> " << this->mV(i, j, k) << endl;
	}

	cout << endl;
	FOR_EACH_FACE_Z {
		cout << "mW i: " << i << " j: " << j << " k: " << k  \
				<< " --> " << this->mW(i, j, k) << endl;
	}


	cout << "#################################" << endl;
#endif
/********************* PRINT(END) ***************************/
}


void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.

	FOR_EACH_CELL {
		vec3 pt(i * theCellSize + theCellSize * 0.5,
				 j * theCellSize + theCellSize * 0.5,
				 k * theCellSize + theCellSize * 0.5);
		vec3 vel = getVelocity(pt);
		pt.n[VX] -= vel[VX] * dt;
		pt.n[VY] -= vel[VY] * dt;
		pt.n[VZ] -= vel[VZ] * dt;
		target.mT(i, j, k) = getTemperature(pt);
	}
	//target.mT = mT;
    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.
	target.mD = mD;
    // Then save the result to our object.


	FOR_EACH_CELL {
		vec3 pt(i * theCellSize + theCellSize * 0.5,
				 j * theCellSize + theCellSize * 0.5,
				 k * theCellSize + theCellSize * 0.5);

		vec3 vel = getVelocity(pt);
		pt.n[VX] -= vel[VX] * dt;
		pt.n[VY] -= vel[VY] * dt;
		pt.n[VZ] -= vel[VZ] * dt;
		target.mD(i, j, k) = getDensity(pt);
	}

    mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.
	//target.mV = mV;
	
	GridData fbuoy;
	fbuoy.initialize();
	double alpha = 1, beta = 1;

	FOR_EACH_CELL {
		if (j < theDim[MACGrid::Y] - 1) {
			fbuoy(i, j + 1, k) = -alpha * (mD(i, j + 1, k) + mD(i, j, k)) / 2 + beta * (mT(i, j + 1, k) + mT(i, j, k));
			target.mV(i, j + 1, k) += fbuoy(i, j + 1, k) * dt;
		}
	}

   // Then save the result to our object.
   mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.


	GridData mUCenter, mVCenter, mWCenter;
	mUCenter.initialize();
	mVCenter.initialize();
	mWCenter.initialize();

	GridData omegaX, omegaY, omegaZ;
	omegaX.initialize();
	omegaY.initialize();
	omegaZ.initialize();

	GridData lapOmegaX, lapOmegaY, lapOmegaZ;
	lapOmegaX.initialize();
	lapOmegaY.initialize();
	lapOmegaZ.initialize();

	GridData nX, nY, nZ;
	nX.initialize();
	nY.initialize();
	nZ.initialize();

	GridData fX, fY, fZ;
	fX.initialize();
	fY.initialize();
	fZ.initialize();

	FOR_EACH_CELL {
		
		vec3 pt(i * theCellSize + 0.5 * theCellSize,
				j * theCellSize + 0.5 * theCellSize,
				k * theCellSize + 0.5 * theCellSize);

		vec3 vel = getVelocity(pt);
		mUCenter(i, j, k) = vel.n[0];
		mVCenter(i, j, k) = vel.n[1];
		mWCenter(i, j, k) = vel.n[2];
	}

	FOR_EACH_CELL {
		omegaX(i, j, k) = (mWCenter(i, j + 1, k) - mWCenter(i, j - 1, k) - mVCenter(i, j, k + 1) + mVCenter(i, j, k - 1)) / (2 * theCellSize);
		omegaY(i, j, k) = (mUCenter(i, j, k + 1) - mUCenter(i, j, k - 1) - mWCenter(i + 1, j, k) + mWCenter(i - 1, j, k)) / (2 * theCellSize);
		omegaZ(i, j, k) = (mVCenter(i + 1, j, k) - mVCenter(i - 1, j, k) - mUCenter(i, j + 1, k) + mUCenter(i, j - 1, k)) / (2 * theCellSize);
	}

	FOR_EACH_CELL {
		lapOmegaX(i, j, k) = (omegaX(i + 1, j, k) - omegaX(i - 1, j, k)) / (2 * theCellSize);
		lapOmegaY(i, j, k) = (omegaY(i, j + 1, k) - omegaY(i, j - 1, k)) / (2 * theCellSize);
		lapOmegaZ(i, j, k) = (omegaZ(i, j, k + 1) - omegaZ(i, j, k - 1)) / (2 * theCellSize);
	}

	FOR_EACH_CELL {
		double norm = sqrt(lapOmegaX(i, j, k) * lapOmegaX(i, j, k) + lapOmegaY(i, j, k) * lapOmegaY(i, j, k) + lapOmegaZ(i, j, k) * lapOmegaZ(i, j, k)) + 1 / pow(10, -20);
		nX(i, j, k) = lapOmegaX(i, j, k) / norm;
		nY(i, j, k) = lapOmegaY(i, j, k) / norm;
		nZ(i, j, k) = lapOmegaZ(i, j, k) / norm;
	}


	/* calculate the fconf: */

	double epsilon = 5;

	FOR_EACH_CELL {
		/* calculate the cross product and then multiply by theCellSize */
		fX(i, j, k) = nY(i, j, k) * omegaZ(i, j, k) - nZ(i, j, k) * omegaY(i, j, k);
		fY(i, j, k) = nZ(i, j, k) * omegaX(i, j, k) - nX(i, j, k) * omegaZ(i, j, k);
		fZ(i, j, k) = nX(i, j, k) * omegaY(i, j, k) - nY(i, j, k) * omegaX(i, j, k);
		fX(i, j, k) *= theCellSize * epsilon;
		fY(i, j, k) *= theCellSize * epsilon;
		fZ(i, j, k) *= theCellSize * epsilon;
	}
	
	FOR_EACH_CELL {
		
		if (i != theDim[MACGrid::X] - 1 && j == theDim[MACGrid::Y] - 1
								 && k == theDim[MACGrid::Z] - 1) {

			target.mU(i + 1, j, k) += (fX(i + 1, j, k) + fX(i, j, k)) / 2 * dt;

		} else if (i == theDim[MACGrid::X] - 1 && j != theDim[MACGrid::Y] - 1
								 && k == theDim[MACGrid::Z] - 1) {

			target.mV(i, j + 1, k) += (fY(i, j + 1, k) + fY(i, j, k)) / 2 * dt;

		} else if (i == theDim[MACGrid::X] - 1 && j == theDim[MACGrid::Y] - 1
								 && k != theDim[MACGrid::Z] - 1) {

			target.mW(i, j, k + 1) += (fZ(i, j, k + 1) + fZ(i, j, k)) / 2 * dt;

		} else if (i != theDim[MACGrid::X] - 1 && j != theDim[MACGrid::Y] - 1
								 && k == theDim[MACGrid::Z] - 1) {
			
			target.mU(i + 1, j, k) += (fX(i + 1, j, k) + fX(i, j, k)) / 2 * dt;
			target.mV(i, j + 1, k) += (fY(i, j + 1, k) + fY(i, j, k)) / 2 * dt;

		} else if (i == theDim[MACGrid::X] - 1 && j != theDim[MACGrid::Y] - 1
								 && k != theDim[MACGrid::Z] - 1) {

			target.mV(i, j + 1, k) += (fY(i, j + 1, k) + fY(i, j, k)) / 2 * dt;
			target.mW(i, j, k + 1) += (fZ(i, j, k + 1) + fZ(i, j, k)) / 2 * dt;

		} else if (i != theDim[MACGrid::X] - 1 && j == theDim[MACGrid::Y] - 1
								 && k != theDim[MACGrid::Z] - 1) {

			target.mU(i + 1, j, k) += (fX(i + 1, j, k) + fX(i, j, k)) / 2 * dt;
			target.mW(i, j, k + 1) += (fZ(i, j, k + 1) + fZ(i, j, k)) / 2 * dt;

		} else if (i != theDim[MACGrid::X] - 1 && j != theDim[MACGrid::Y] - 1
								 && k != theDim[MACGrid::Z] - 1) {

			target.mU(i + 1, j, k) += (fX(i + 1, j, k) + fX(i, j, k)) / 2 * dt;
			target.mV(i, j + 1, k) += (fY(i, j + 1, k) + fY(i, j, k)) / 2 * dt;
			target.mW(i, j, k + 1) += (fZ(i, j, k + 1) + fZ(i, j, k)) / 2 * dt;

		}
	}

	//target.mU = mU;
	//target.mV = mV;
	//target.mW = mW;
	
	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}


void MACGrid::computeDrivingForce(double dt)
{
	/** IMPORTANT !! */
	/* if you wanna see the driving force event, please make sure the
	 * dimension is 2-D not 3-D. which means theDim[MACGrid::Z] = 1;
	 */
	Convolve(mDensityBlur, mD, mFilter, mFilterWidth);
	//Convolve(mTargetDensityBlur, mTargetDensity, mFilter, mFilterWidth);
	cout << "======================================================" << endl;
	cout << "************ mDensityBlur  ***********" << endl;
	FOR_EACH_CELL {
		cout << "i: " << i << " j: " << j << " k: " << k << "   " << mDensityBlur(i, j, 0) << endl;
	}
	cout << "************ mTargetDensityBlur **********" << endl;
	FOR_EACH_CELL {
		if (i == 0 && j == 0) {
			continue;
		} else if (i != 0 && j == 0) { 		//only update the face on Z direction

			double fx = (mDensityBlur(i, j, 0) + mDensityBlur(i - 1, j, 0)) * 0.5 \
					* ( (mTargetDensityBlur(i, j, 0) - mTargetDensityBlur(i - 1, j, 0)) / theCellSize ) \
					/ ( ( mTargetDensityBlur(i, j, 0) + mTargetDensityBlur(i - 1, j, 0) ) * 0.5 );
					
			target.mU(i, j, 0) += vf * fx * dt;

		} else if (i == 0 && j != 0) {	//only update the face on Y direction 

			double fy = (mDensityBlur(i, j, 0) + mDensityBlur(i, j - 1, 0)) * 0.5 \
					* ( (mTargetDensityBlur(i, j, 0) - mTargetDensityBlur(i, j - 1, 0)) / theCellSize ) \
					/ ( ( mTargetDensityBlur(i, j, 0) + mTargetDensityBlur(i, j - 1, 0) ) * 0.5 );
			
			target.mV(i, j, 0) += vf * (mTargetDensityBlur(i, j, 0) - mTargetDensityBlur(i, j - 1, 0)) * dt;

		} else if (i != 0 && j != 0) {

			double fx = (mDensityBlur(i, j, 0) + mDensityBlur(i - 1, j, 0)) * 0.5 \
					* ( (mTargetDensityBlur(i, j, 0) - mTargetDensityBlur(i - 1, j, 0)) / theCellSize ) \
					/ ( ( mTargetDensityBlur(i, j, 0) + mTargetDensityBlur(i - 1, j, 0) ) * 0.5 );

			double fy = (mDensityBlur(i, j, 0) + mDensityBlur(i, j - 1, 0)) * 0.5 \
					* ( (mTargetDensityBlur(i, j, 0) - mTargetDensityBlur(i, j - 1, 0)) / theCellSize ) \
					/ ( ( mTargetDensityBlur(i, j, 0) + mTargetDensityBlur(i, j - 1, 0) ) * 0.5 );

			target.mU(i, j, 0) += vf * fx * dt;
			target.mV(i, j, 0) += vf * fy * dt;
		}
		cout << "i: " << i << " j: " << j << " k: " << k << "   " << mTargetDensityBlur(i, j, 0) << endl;
	}
	mU = target.mU;
	mV = target.mV;
	
	FOR_EACH_CELL {
		if (i == 0 && j == 0) {
			continue;
		} else if (i != 0 && j == 0) {
			target.mU(i, j, k) -= target.mU(i, j, k) * vd;
		} else if (i == 0 && j != 0) {
			target.mV(i, j, k) -= target.mV(i, j, k) * vd;
		} else if (i != 0 && j != 0) {
			target.mU(i, j, k) -= target.mU(i, j, k) * vd;
			target.mV(i, j, k) -= target.mV(i, j, k) * vd;
		}
	}

	mU = target.mU;
	mV = target.mV;
}

void MACGrid::addExternalForces(double dt)
{
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
   //if (runTimes <= 0)
   //		computeDrivingForce(dt);
}


int MACGrid::getNonSolid(int i, int j, int k)
{
	int count = 0;
	count = i == 0 ? count + 1 : count;
	count = i == theDim[X] ? count + 1 : count;

	count = j == 0 ? count + 1 : count;
	count = j == theDim[Y] ? count + 1 : count;

	count = k == 0 ? count + 1 : count;
	count = k == theDim[Z] ? count + 1 : count;

	return count;
}

bool MACGrid::checkDivergence()
{
	bool ret = true;
	double threshold = 1e-4;
	double divergence = 0;
	FOR_EACH_CELL {

		vec3 pt = getCenter(i, j, k);
		pt.Print("The position");
		
		cout << "====================================== " << endl;
		cout << "i: " << i << " j: " << j << " k: " << k << endl;
		double diffX = getVelocityX(vec3(pt.n[0] + theCellSize * 0.5, pt.n[1], pt.n[2]))
						- getVelocityX(vec3(pt.n[0] - theCellSize * 0.5, pt.n[1], pt.n[2]));
		cout << "diffX: " << diffX << endl;

		double diffY = getVelocityY(vec3(pt.n[0], pt.n[1] + theCellSize * 0.5 , pt.n[2]))
						- getVelocityY(vec3(pt.n[0], pt.n[1] - theCellSize * 0.5, pt.n[2]));

		cout << "diffY: " << diffY << endl;

		double diffZ = getVelocityZ(vec3(pt.n[0], pt.n[1], pt.n[2] + theCellSize * 0.5))
						- getVelocityZ(vec3(pt.n[0], pt.n[1], pt.n[2] - theCellSize * 0.5));
		cout << "diffZ: " << diffZ << endl;

		cout << "diffX + diffY + diffZ: " << diffX + diffY + diffZ << endl;
		divergence = diffX + diffY + diffZ;

		if (divergence  > threshold) {
			//cout << "DIVERGENCE is larger than threshold in cell i: " << i << " j: "
			//			<< j  << " k: " << k << endl;
			cout << "Larger than threshold" << endl;
			cout << threshold << endl;
			cout << "Divergence = " << divergence << endl;
		}
		
	}
	return ret;
	
}

void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	// 1. Construct d
	// 2. Construct A
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target.
	GridData d, p;
	d.initialize();
	p.initialize();
	double rho = 1;

	FOR_EACH_CELL {
		//double rho = getDensity(vec3(i * theCellSize + 0.5 * theCellSize,
		//							 j * theCellSize + 0.5 * theCellSize,
		//							 k * theCellSize + 0.5 * theCellSize));

		vec3 ptX(i * theCellSize,
				 j * theCellSize + 0.5 * theCellSize,
				 k * theCellSize + 0.5 * theCellSize);
		vec3 ptXplusOne((i + 1) * theCellSize,
						 j * theCellSize + 0.5 * theCellSize,
						 k * theCellSize + 0.5 * theCellSize);

		vec3 ptY(i * theCellSize + 0.5 * theCellSize,
				 j * theCellSize,
				 k * theCellSize + 0.5 * theCellSize);
		vec3 ptYplusOne(i * theCellSize + 0.5 * theCellSize,
						 (j + 1) * theCellSize,
						 k * theCellSize + 0.5 * theCellSize);
		
		vec3 ptZ(i * theCellSize + 0.5 * theCellSize,
				 j * theCellSize + 0.5 * theCellSize,
				 k * theCellSize);
		vec3 ptZplusOne(i * theCellSize + 0.5 * theCellSize,
						j * theCellSize + 0.5 * theCellSize,
						(k + 1) * theCellSize);
		
		double u_right = getVelocityX(ptX);
		double u_left = getVelocityX(ptXplusOne);
		double v_down = getVelocityY(ptY);
		double v_up = getVelocityY(ptYplusOne);
		double w_front = getVelocityZ(ptZ);
		double w_back = getVelocityZ(ptZplusOne);


/*********************** PRINT (BEGIN) ********************************/
#ifdef _MAC_GRID_DEBUG
		cout << "----------" << i << " " << j << " " << k << "-------------" << endl;
		cout << "u_left = " << u_left << endl;
		cout << "u_right = " << u_right << endl;
		cout << "v_up = " << v_up << endl;
		cout << "v_down = " << v_down << endl;
		cout << "w_back = " << w_back << endl;
		cout << "w_fron = " << w_front << endl;
#endif
/*********************** PRINT (END) ********************************/

		d(i, j, k) = - (u_left - u_right + v_up - v_down + w_back - w_front) * rho * theCellSize / dt;
	}

/*********************** PRINT (BEGIN) ********************************/
#ifdef _MAC_GRID_DEBUG
	cout << "*********************" << endl;
	FOR_EACH_CELL {
		cout << "d(" << i << ", " << j << ", " << k << ") = " << d(i, j, k) << endl;
	}
#endif
/*********************** PRINT (END) ********************************/

	conjugateGradient(this->AMatrix, p, d, 100, 1 / 1000000);
	

/*********************** PRINT (BEGIN) ********************************/

#ifdef _MAC_GRID_DEBUG
	cout << endl;
	FOR_EACH_CELL {
		cout << "p(" << i << ", " << j << ", " << k << ") = " << p(i, j, k) << endl;
	}

	cout << "#####################" << endl;
	cout << endl;
#endif
	
/*********************** PRINT (END) ********************************/

	FOR_EACH_CELL {
		
		if (i == 0 && j == 0 && k == 0) {
			continue;
		} else if (i == 0 && j == 0 && k != 0) { 		//only update the face on Z direction

			target.mW(i, j, k) -= dt * (p(i, j, k) - p(i, j, k - 1)) / theCellSize;

		} else if (i == 0 && j != 0 && k == 0) {	//only update the face on Y direction 

			target.mV(i, j, k) -= dt * (p(i, j, k) - p(i, j - 1, k)) / theCellSize;

		} else if (i != 0 && j == 0 && k == 0) {

			target.mU(i, j, k) -= dt * (p(i, j, k) - p(i - 1, j, k)) / theCellSize;

		} else if (i == 0 && j != 0 && k != 0) {

			target.mV(i, j, k) -= dt * (p(i, j, k) - p(i, j - 1, k)) / theCellSize;
			target.mW(i, j, k) -= dt * (p(i, j, k) - p(i, j, k - 1)) / theCellSize;

		} else if (i != 0 && j == 0 && k != 0) {

			target.mU(i, j, k) -= dt * (p(i, j, k) - p(i - 1, j, k)) / theCellSize;
			target.mW(i, j, k) -= dt * (p(i, j, k) - p(i, j, k - 1)) / theCellSize;

		} else if (i != 0 && j != 0 && k == 0) {

			target.mU(i, j, k) -= dt * (p(i, j, k) - p(i - 1, j, k)) / theCellSize;
			target.mV(i, j, k) -= dt * (p(i, j, k) - p(i, j - 1, k)) / theCellSize;

		} else if (i != 0 && j != 0 && k != 0) {
			
			target.mU(i, j, k) -= dt * (p(i, j, k) - p(i - 1, j, k)) / theCellSize;
			target.mV(i, j, k) -= dt * (p(i, j, k) - p(i, j - 1, k)) / theCellSize;
			target.mW(i, j, k) -= dt * (p(i, j, k) - p(i, j, k - 1)) / theCellSize;
		}
		
/*********************** PRINT (BEGIN) ********************************/
#ifdef _MAC_GRID_DEBUG
		cout << "--------" << i << " " << j << " " << k << "-----------" << endl;
		cout << "mU = " << target.mU(i, j, k) << endl;
		cout << "mV = " << target.mV(i, j, k) << endl;
		cout << "mW = " << target.mW(i, j, k) << endl;
		cout << "------------------------------------------------------" << endl;
#endif
/*********************** PRINT (END) ********************************/
	}



//	FOR_EACH_CELL {
//		target(i, j, k) = p(i, j, k);
//	}
	
	//target.mP = mP;
	//target.mU = mU;
	//target.mV = mV;
	//target.mW = mW;
	

	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
	//checkDivergence();
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}







/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:
	z = r;

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true;
		}

		// TODO: Apply a preconditioner here.
		// For now, just bypass the preconditioner:
		z = r;

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	//PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();
           vel *= theCellSize/2.0;
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(22.0, 2.0, 1.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(255.0, 0.0, 0.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
