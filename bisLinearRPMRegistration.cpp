/*  License
 
 _This file is Copyright 2018 by the Image Processing and Analysis Group (BioImage Suite Team). Dept. of Radiology & Biomedical Imaging, Yale School of Medicine._ It is released under the terms of the GPL v2.
 
 ----
     
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
   See also  http: www.gnu.org/licenses/gpl.html
   
   If this software is modified please retain this statement and add a notice
   that it had been modified (and by whom).  
 
 Endlicense */


#include "bisLinearRPMRegistration.h"
#include "bisUtil.h"

bisLinearRPMRegistration::bisLinearRPMRegistration(std::string n) : bisRPMCorrespondenceFinder(n) {
  this->Output=new bisMatrixTransformation();
}

bisLinearRPMRegistration::~bisLinearRPMRegistration()
{
  if (this->Output)
    delete this->Output;
}


  // Description:
  // get Output Transformation  
bisSimpleMatrix<float>* bisLinearRPMRegistration::getOutputMatrix() {
  return this->Output->getSimpleMatrix("lin_reg_matrix");
}

/** run transformation
 * @param transformMode (0=rigid,1=similarity,2=affine)
 * @param initialTemperature initial temperature
 * @param finalTemperature final temperature
 * @param annealRate anneal rate for RPM
 * @param useCentroids if true center points first
 * @param initialTransformaiton use this to initialize the mappings
 * @param debug if true print extra messages
 * @returns 1 if success, 0 if failed
 */
int bisLinearRPMRegistration::run(int in_transformMode,
                                  int in_correspondenceMode,
                                  float in_initialTemperature,
                                  float in_finalTemperature,
                                  float in_annealRate,
                                  int in_useCentroids,
                                  bisMatrixTransformation* in_initialTransformation,
                                  int in_debug) {


  if (this->locator==NULL) {
    std::cerr << "___ bisLinearRPMRegistration not initialized" << std::endl;
    return 0;
  }
  
  int TransformMode=bisUtil::irange(in_transformMode,0,2);
  int CorrespondenceMode=bisUtil::irange(in_correspondenceMode,0,2);
  float FinalTemperature=bisUtil::frange(in_finalTemperature,0.01,1000.0);
  float InitialTemperature=bisUtil::frange(in_initialTemperature,FinalTemperature,1000.0);
  float AnnealRate=bisUtil::frange(in_annealRate,0.5,0.999);
  int debug=in_debug;

  // Initial Mapping including centroid shift!
  bisUtil::mat44 m;
  if (!in_initialTransformation) {
    bisUtil::makeIdentityMatrix(m);
    if (in_useCentroids) {
      float cx[3],cy[3];
      bisPointRegistrationUtils::computeCentroid(this->SampledReferencePoints.get(),cx,debug);
      bisPointRegistrationUtils::computeCentroid(this->SampledTargetPoints.get(),cy,debug);
      if (debug)
        std::cout << "___ using initial centroid alignment :";
      for (int ia=0;ia<=2;ia++) {
        m[ia][3]=cy[ia]-cx[ia];
        if (debug)
          std::cout << m[ia][3] << " ";
      }
      
      if (debug)
        std::cout << std::endl;
    } else {
      if (debug)
        std::cout << "Not using centroids " << std::endl;
    }
    this->Output->setMatrix(m);
  }  else {
    if (debug)
      std::cout << "Using Initial Transformation" << std::endl;
    in_initialTransformation->getMatrix(m);
    this->Output->setMatrix(m);
  }
    bisPointRegistrationUtils::printMatrix(this->Output,"Initial Mapping");

  bisSimpleMatrix<float>* OutputRefLandmarks=new bisSimpleMatrix<float>();
  bisSimpleMatrix<float>* OutputTargetLandmarks=new bisSimpleMatrix<float>();
  bisSimpleVector<float>* OutputWeights=new bisSimpleVector<float>();
                                       
  
  float Temperature=InitialTemperature;
  int numpoints=this->SampledReferencePoints->getNumRows();
  int numpoints2=this->SampledTargetPoints->getNumRows();
    
  if (debug) {
    std::cout << "___ Beginning  Linear RPM Registration tmode=" << TransformMode << " cmode=" << CorrespondenceMode << std::endl;
    std::cout << "___            NumPoints= " << numpoints << " vs " << numpoints2 << " temperatures=" << InitialTemperature << ":" << AnnealRate << ":" << FinalTemperature << std::endl;
  }

  if (debug) {
    bisPointRegistrationUtils::printTwoPoints(this->SampledReferencePoints.get(),"ref");
    bisPointRegistrationUtils::printTwoPoints(this->SampledTargetPoints.get(),"target");
  }

  
  
  int iteration=0;
  while (Temperature > FinalTemperature) {
    iteration=iteration+1;
    if (debug) 
      std::cout << "Beginning iteration " << iteration << " t=" << Temperature << "  TransformMode=" << TransformMode << std::endl;


    this->estimateCorrespondence(this->Output,
                                 Temperature,
                                 CorrespondenceMode,
                                 OutputRefLandmarks,
                                 OutputTargetLandmarks,
                                 OutputWeights,
                                 debug);

    if (debug) {
      bisPointRegistrationUtils::printJointPoints(OutputRefLandmarks,OutputTargetLandmarks,OutputWeights,"out_ref->targ",117);
      bisPointRegistrationUtils::printTwoPoints(OutputRefLandmarks,"ref");
      bisPointRegistrationUtils::printTwoPoints(OutputTargetLandmarks,"targ");
    }
    bisPointRegistrationUtils::computeLandmarkTransformation(OutputRefLandmarks,
                                                             OutputTargetLandmarks,
                                                             TransformMode,
                                                             this->Output,
                                                             OutputWeights,
                                                             debug);

    if (debug)
      bisPointRegistrationUtils::printMatrix(this->Output,"End");
    
    Temperature*=AnnealRate;

  }


  delete OutputRefLandmarks;
  delete OutputTargetLandmarks;
  delete OutputWeights;

  return iteration;
}


