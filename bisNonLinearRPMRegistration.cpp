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


#include "bisNonLinearRPMRegistration.h"
#include "bisApproximateLandmarkDisplacementsWithGridTransform.h"
#include "bisJSONParameterList.h"
#include "bisUtil.h"

bisNonLinearRPMRegistration::bisNonLinearRPMRegistration(std::string n) : bisRPMCorrespondenceFinder(n) {

  std::shared_ptr<bisComboTransformation> tmp(new bisComboTransformation());
  this->Output=std::move(tmp);
  
  std::shared_ptr<bisGridTransformation> tmpgrid(new bisGridTransformation());
  this->Grid=std::move(tmpgrid);
  this->Output->addTransformation(this->Grid);
  
  tmpgrid=0;
  tmp=0;
}

bisNonLinearRPMRegistration::~bisNonLinearRPMRegistration()
{
  this->Output=0;
  this->Grid=0;
}


  // Description:
  // get Output Transformation  
std::shared_ptr<bisComboTransformation> bisNonLinearRPMRegistration::getOutput() {
  return this->Output;
}

// Description:
int bisNonLinearRPMRegistration::initializeWithLinear(bisMatrixTransformation* in_initialTransformation,
                                                      bisSimpleMatrix<float>* source,
                                                      bisSimpleMatrix<float>* target,
                                                      int maxnumlandmarks,
                                                      int samplingweight, // 1 = equal sampling, if  > 1 get extra points from points whose label > 0
                                                      bisSimpleVector<int>* source_labels,
                                                      bisSimpleVector<int>* target_labels,
                                                      int debug) {


  // Two parts
  // 1. Set the linear component of Output to in_initialTransformation
  // 2. Map target points by Inverse of in_initialTransformation and use these to initialize the locator
  this->Output->setInitialTransformation(in_initialTransformation);

  // Compute inverse
  std::unique_ptr<bisMatrixTransformation> tmp(new bisMatrixTransformation());
  bisUtil::mat44 m; in_initialTransformation->getMatrix(m);

  Eigen::MatrixXf eigen=Eigen::MatrixXf::Zero(4,4);
  for (int ia=0;ia<=3;ia++)
    for (int ib=0;ib<=3;ib++)
      eigen(ia,ib)=m[ia][ib];

  Eigen::MatrixXf inverse=eigen.inverse();
  for (int ia=0;ia<=3;ia++)
    for (int ib=0;ib<=3;ib++)
      m[ia][ib]=inverse(ia,ib);
  tmp->setMatrix(m);

  bisSimpleMatrix<float>* mappedPoints=bisPointRegistrationUtils::transformPoints(source,tmp.get(),debug);
  return bisRPMCorrespondenceFinder::initialize(source,mappedPoints,maxnumlandmarks,samplingweight,source_labels,target_labels,debug);
  
  
}

int bisNonLinearRPMRegistration::initializeGrid(float spacing,float cps_begin,float cps_end,bisSimpleMatrix<float>* RefLandmarks) {

  float current_spa[3]; this->Grid->getGridSpacing(current_spa);
  int current_dim[3]; this->Grid->getGridDimensions(current_dim);

  if (current_dim[0]>0) {
    // Do we need to reshape the grid, may be it is close enough
    float gap=0.1*fabs(cps_begin-cps_end);
    float dg=fabs(spacing-current_spa[0]);
    if (dg<gap)
      return 0;
  }

  // Find bounds
  float minc[3],maxc[3];
  float* pts=RefLandmarks->getData();
  int numpoints=RefLandmarks->getNumRows();

  // Find min and max and create bounds
  for (int ia=0;ia<=2;ia++) {
    minc[ia]=pts[ia];
    maxc[ia]=pts[ia];
  }

  for (int index=1;index<numpoints;index++) {
    int offset=index*3;
    for (int ia=0;ia<=2;ia++) {
      float p=pts[offset+ia];
      if (minc[ia]>p) minc[ia]=p;
      if (maxc[ia]<p) maxc[ia]=p;
    }
  }

  // Compute Grid Dimensions, origin and spacing
  int dim[3];
  float ori[3],spa[3];
  
  for (int ia=0;ia<=2;ia++) {
    float l=maxc[ia]-minc[ia];
    dim[ia]=int(l/spacing)+1;
    float l2=(dim[ia]-1)*spacing;
    ori[ia]=minc[ia]-(l2-l)*0.5;
    spa[ia]=spacing;
  }

  // Initialize the grid
  Grid->initializeGrid(dim,spa,ori,1);
  return 1;

}

int bisNonLinearRPMRegistration::approximateGrid(bisSimpleMatrix<float>* sourcePts,
                                                 bisSimpleMatrix<float>* targetPts,
                                                 bisSimpleVector<float>* weights,
                                                 float Smoothness,
                                                 int accurate,int debug) {


  int iterations=20;
  int steps=1;
  float stepsize=1.0;
  float tolerance=0.2;
  
  if (accurate) {
    steps=2;
    tolerance=0.05;
    stepsize=0.5;
  }

  std::unique_ptr<bisJSONParameterList> plist(new bisJSONParameterList());
  plist->setFloatValue("lambda",Smoothness);
  plist->setFloatValue("stepsize",stepsize);
  plist->setFloatValue("tolerance",tolerance);
  plist->setIntValue("steps",steps);
  plist->setIntValue("iterations",iterations);
  
  std::unique_ptr<bisApproximateLandmarkDisplacementsWithGridTransform> landmarkFit(new bisApproximateLandmarkDisplacementsWithGridTransform());
  return landmarkFit->run(sourcePts,targetPts,weights,this->Grid.get(),plist.get(),debug);
}
  
int bisNonLinearRPMRegistration::run(int in_correspondenceMode,
                                     float in_cps_begin,
                                     float in_cps_end,
                                     float in_smoothness_begin,
                                     float in_smoothness_end,
                                     float in_initialTemperature,
                                     float in_finalTemperature,
                                     int   in_iterationPerTemperature,
                                     float in_annealRate,
                                     int in_debug) {


  if (this->locator==NULL) {
    std::cerr << "___ bisNonLinearRPMRegistration not initialized" << std::endl;
    return 0;
  }

  float CPSBegin=bisUtil::frange(in_cps_begin,1.0,100.0);
  float CPSEnd=bisUtil::frange(in_cps_end,0.2,CPSBegin);
  float SmoothnessBegin=bisUtil::frange(in_smoothness_begin,0.001,10.0);
  float SmoothnessEnd=bisUtil::frange(in_smoothness_end,0.0001,SmoothnessBegin);

  
  int IterationPerTemperature = bisUtil::irange(in_iterationPerTemperature,1,10);
  int CorrespondenceMode=bisUtil::irange(in_correspondenceMode,0,2);
  float FinalTemperature=bisUtil::frange(in_finalTemperature,0.01,1000.0);
  float InitialTemperature=bisUtil::frange(in_initialTemperature,FinalTemperature,1000.0);
  float AnnealRate=bisUtil::frange(in_annealRate,0.5,0.999);
  int debug=in_debug;
  
  bisSimpleMatrix<float>* OutputRefLandmarks=new bisSimpleMatrix<float>();
  bisSimpleMatrix<float>* OutputTargetLandmarks=new bisSimpleMatrix<float>();
  bisSimpleVector<float>* OutputWeights=new bisSimpleVector<float>();
  
  
  float Temperature=InitialTemperature;
  float CPS=CPSBegin;
  float Smoothness=SmoothnessBegin;
  int numpoints=this->SampledReferencePoints->getNumRows();
  int numpoints2=this->SampledTargetPoints->getNumRows();

  int totaliter=int(fabs(log(InitialTemperature/FinalTemperature)/log(AnnealRate))+1.0);
  float cps_annealrate=exp(-fabs(log(CPSBegin/CPSEnd)/float(totaliter-1.0)));
  float smooth_annealrate=exp(-fabs(log(SmoothnessBegin/SmoothnessEnd)/float(totaliter-1.0)));
  
  if (debug) {
    std::cout << "___ Beginning  NonLinear RPM Registration tmode=" << " cmode=" << CorrespondenceMode << std::endl;
    std::cout << "___            CPS= " << CPSBegin << ":" << CPSEnd << " , smoothnesss=" << SmoothnessBegin << ":" << SmoothnessEnd << std::endl;
    std::cout << "___            NumPoints= " << numpoints << " vs " << numpoints2 << " temperatures=" << InitialTemperature << ":" << AnnealRate << ":" << FinalTemperature << std::endl;
  }

  if (debug) {
    bisPointRegistrationUtils::printTwoPoints(this->SampledReferencePoints.get(),"ref");
    bisPointRegistrationUtils::printTwoPoints(this->SampledTargetPoints.get(),"target");
  }

  
  
  int iteration_count=0;
  totaliter*=IterationPerTemperature;
  while (Temperature > FinalTemperature)
    {
      for (int it=1;it<=IterationPerTemperature;it++) {
        iteration_count=iteration_count+1;
        if (debug) 
          std::cout << "Beginning iteration " << iteration_count << "/" << totaliter << ". Temp=" << Temperature << "  CPS=" << CPS << " Smooth=" << Smoothness << std::endl;
        
        
        this->estimateCorrespondence(this->Output.get(),
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

        this->initializeGrid(CPS,CPSBegin,CPSEnd,OutputRefLandmarks);
        int accurate=0;
        // Do accurate at last step
        if (it==IterationPerTemperature && Temperature*AnnealRate< FinalTemperature)
          accurate=1;
        this->approximateGrid(OutputRefLandmarks,
                              OutputTargetLandmarks,
                              OutputWeights,
                              Smoothness,
                              accurate,debug);
        
      }
      Temperature*=AnnealRate;
      CPS*=cps_annealrate;
      Smoothness*=smooth_annealrate;
    }
  
  delete OutputRefLandmarks;
  delete OutputTargetLandmarks;
  delete OutputWeights;

  return iteration_count;
}


