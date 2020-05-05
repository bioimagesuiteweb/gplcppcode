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
#include "bisGPLTestFunctions.h"
#include "bisSimpleDataStructures.h"
#include "bisJSONParameterList.h"
#include "bisMatrixTransformation.h"
#include "bisGridTransformation.h"
#include "bisOptimizer.h"
#include "bisApproximateLandmarkDisplacementsWithGridTransform.h"
#include "bisRPMCorrespondenceFinder.h"
#include "bisPointLocator.h"
#include <memory>



class bisTestOptimizable : public bisOptimizableAlgorithm {

public:
  bisTestOptimizable() : bisOptimizableAlgorithm("test") { };
  
  float computeGradient(std::vector<float>& params,std::vector<float>& grad) {

    int num=params.size();
    if (num==2) {
      
      float x=params[0],y=params[1];
      float dx=2*(x-9), dy=2*y;
      float s=float(sqrt(dx*dx+dy*dy)+0.00001);
      grad[0]=dx/s; grad[1]=dy/s;
      return s;
    }
    
    double x=params[0], dx=2.0f*(x-9.0f);
    grad[0]=float(dx/fabs(dx));
    return float(fabs(dx));
  };
  
  float computeValue(std::vector<float>& params) {

    int num=params.size();
    float x=params[0];
    float v=0.0,y=0.0;
    if (num==2)
      y=params[1];
    v=(x-9)*(x-9)+y*y;
    return v;
  };
  
  float comparePos(std::vector<float> p,std::vector<float> tp) {
    float sum=0.0;
    for (unsigned int ia=0;ia<p.size();ia++)
      sum+=powf(p[ia]-tp[ia],2.0f);
    return sum;
  }
};


int  test_optimizer(int numparam) {

  int numdof=numparam;
  int numfail=0;
  
  for (int mode=0;mode<=2;mode++)
    {
      std::unique_ptr<bisTestOptimizable> test_optimizable(new bisTestOptimizable());
      std::unique_ptr<bisOptimizer> optimizer(new bisOptimizer(test_optimizable.get()));

      std::vector<float> position(numdof);
      position[0]=15;
      if (numdof>1)
	position[1]=5;
      
      std::vector<float> truepos(numdof);
      truepos[0]=9;
      if (numdof>1)
	truepos[1]=0;


      if (position.size()==2)
	std::cout << std::endl << "________________ mode=" << mode << " pos="  << position[0] << ", " << position[1] << std::endl;
      else
	std::cout << std::endl << "________________ mode=" << mode << " pos="  << position[0] << ", " << std::endl;
      
      if (mode==0) 
	optimizer->computeSlowClimb(position,0.5,25);
      else if (mode==1)
	optimizer->computeGradientDescent(position,25,0.01f);
      else
	optimizer->computeConjugateGradient(position,25,0.01f);
      
      float d=test_optimizable->comparePos(position,truepos);

      if (d>0.1)
	numfail++;

      if (position.size()==2)
	std::cout << "\t\t Final " << position[0] << ", " << position[1] << ". Diff2=" << d << " numfail=" << numfail << std::endl;
      else
      	std::cout << "\t\t Final " << position[0] << ". Diff2=" << d << " numfail=" << numfail << std::endl;
      
      if (d>0.1)
	return numfail;
      

    }
  return numfail;
}

// ------------------------------------------------------------------------

unsigned char*  test_landmarkApproximationWASM(unsigned char* in_source_ptr,
                                               unsigned char* in_target_ptr,
                                               const char* jsonstring,
                                               int debug)
{

  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  float spacing=params->getFloatValue("spacing",20.0);
  if (debug)
    std::cout << "___ spacing=" << spacing << std::endl;

  std::unique_ptr<bisSimpleMatrix<float> > source(new bisSimpleMatrix<float>("source_points_json"));
  if (!source->linkIntoPointer(in_source_ptr))
    return 0;

  if (debug)
    std::cout << "___ Ref Allocated = " << source->getNumRows() << "*" << source->getNumCols() << std::endl;
  
  std::unique_ptr<bisSimpleMatrix<float> > target(new bisSimpleMatrix<float>("target_points_json"));
  if (!target->linkIntoPointer(in_target_ptr))
    return 0;

  if (debug) 
    std::cout << "___ Target Allocated = " << target->getNumRows() << "*" << target->getNumCols() << std::endl;

  int numpoints=source->getNumRows();
  std::unique_ptr<bisSimpleVector<float> > weights(new bisSimpleVector<float>("weights_points_json"));
  weights->allocate(numpoints);
  weights->fill(1.0);

  if (debug) 
    std::cout << "___ Weights Allocated = " << weights->getLength() << std::endl;

  
  float minc[3],maxc[3];
  float* pts=source->getData();

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

  int dim[3];
  float ori[3],spa[3];
  
  for (int ia=0;ia<=2;ia++) {
    float l=maxc[ia]-minc[ia];
    dim[ia]=int(l/spacing)+1;
    float l2=(dim[ia]-1)*spacing;
    ori[ia]=minc[ia]-(l2-l)*0.5;
    spa[ia]=spacing;
  }

  if (debug) {
    std::cout << "__ Creating Grid " << std::endl;
    std::cout << "__    min bounds " << minc[0] << "," << minc[1] << "," << minc[2] << std::endl;
    std::cout << "__    max bounds " << maxc[0] << "," << maxc[1] << "," << maxc[2] << std::endl;
    std::cout << "__    grid dim   " << dim[0] << "," << dim[1] << "," << dim[2] << std::endl;
    std::cout << "__    grid spa   " << spa[0] << "," << spa[1] << "," << spa[2] << std::endl;
    std::cout << "__    grid ori   " << ori[0] << "," << ori[1] << "," << ori[2] << std::endl;
  }
  std::unique_ptr<bisGridTransformation> grid(new bisGridTransformation());
  grid->initializeGrid(dim,spa,ori,1);
  
  std::unique_ptr<bisApproximateLandmarkDisplacementsWithGridTransform> landmarkFit(new bisApproximateLandmarkDisplacementsWithGridTransform());
  landmarkFit->run(source.get(),target.get(),weights.get(),grid.get(),params.get(),debug);

  if (debug)
    std::cout << "___ Transforming points " << std::endl;
  
  std::unique_ptr<bisSimpleMatrix<float> > output(new bisSimpleMatrix<float>("warped_points_json"));
  output->zero(source->getNumRows(),3);

  float* srcPts=source->getData();
  float* outPts=output->getData();
  for (int i=0;i<source->getNumRows();i++) {
    float x1[3],x2[3];
    for (int ia=0;ia<=2;ia++) 
      x1[ia]=srcPts[i*3+ia];
    grid->transformPoint(x1,x2);
    for (int ia=0;ia<=2;ia++)
      outPts[i*3+ia]=x2[ia];
  }
    
      
  return output->releaseAndReturnRawArray();

}
// ------------------------------------------------------------------------

unsigned char*  test_rpmCorrespondenceEstimatorWASM(unsigned char* in_source,
                                                    unsigned char* in_target,
                                                    const char* jsonstring,
                                                    int debug)
{
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  float temperature=params->getFloatValue("temperature",10.0);
  int mode=params->getIntValue("mode",2);
  int numlandmarks=params->getIntValue("numpoints",1000);
    
  if (debug)
    std::cout << "___ temperature=" << temperature << " mode=" << mode << " numpoints=" << numlandmarks << std::endl;

  std::unique_ptr<bisSimpleMatrix<float> > source(new bisSimpleMatrix<float>("source_points_json"));
  if (!source->linkIntoPointer(in_source))
    return 0;

  if (debug)
    std::cout << "___ Ref Allocated = " << source->getNumRows() << "*" << source->getNumCols() << std::endl;
  
  std::shared_ptr<bisSimpleMatrix<float> > target(new bisSimpleMatrix<float>("target_points_json"));
  if (!target->linkIntoPointer(in_target))
    return 0;

  if (debug) 
    std::cout << "___ Target Allocated = " << target->getNumRows() << "*" << target->getNumCols() << std::endl;

  std::unique_ptr<bisPointLocator> locator(new bisPointLocator());
  locator->initialize(target,0.2,0);

  std::unique_ptr<bisMatrixTransformation> matrix(new bisMatrixTransformation());
  matrix->identity();

  int numref=source->getNumRows();
  int numtarget=target->getNumRows();  
  std::unique_ptr<bisSimpleMatrix<float> > OutputRefLandmarks(new bisSimpleMatrix<float>());
  std::unique_ptr<bisSimpleMatrix<float> > OutputTargetLandmarks(new bisSimpleMatrix<float>());
  std::unique_ptr<bisSimpleVector<float> > OutputWeights(new bisSimpleVector<float>());
  OutputRefLandmarks->zero(numref,3);
  OutputTargetLandmarks->zero(numref,3);
  OutputWeights->zero(numref);


  std::unique_ptr<bisSimpleVector<int> > RefLabels(new bisSimpleVector<int>());
  RefLabels->zero(numref);

  std::unique_ptr<bisSimpleVector<int> > TargLabels(new bisSimpleVector<int>());
  TargLabels->zero(numtarget);
  

  if (mode == 0) {
    bisRPMCorrespondenceFinder::computeCorrespodnencesICP(matrix.get(),locator.get(),
                                                          source->getData(),
                                                          OutputRefLandmarks->getData(),
                                                          OutputTargetLandmarks->getData(),
                                                          OutputWeights->getData(),
                                                          numref,debug);

    return OutputTargetLandmarks->releaseAndReturnRawArray();
  }


  Eigen::SparseMatrix<float,Eigen::RowMajor> M(numref,numtarget);
  bisRPMCorrespondenceFinder::computeCorrespondencesRPM(matrix.get(),locator.get(),M,
                                                        mode,
                                                        source->getData(),RefLabels->getData(),
                                                        target->getData(),TargLabels->getData(),
                                                        OutputRefLandmarks->getData(),
                                                        OutputTargetLandmarks->getData(),
                                                        OutputWeights->getData(),
                                                        temperature,numref,numtarget,debug);

  //  std::cout << "M=" << M << std::endl;
  
  return OutputTargetLandmarks->releaseAndReturnRawArray();
}

unsigned char*  test_rpmSamplingWASM(unsigned char* in_points,
                                     unsigned char* in_labels,
                                     const char* jsonstring,
                                     int debug)
{

  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  int prefsampling=params->getIntValue("prefsampling",2);
  int numlandmarks=params->getIntValue("numpoints",1000);

    
  if (debug)
    std::cout << " prefsampling=" << prefsampling << " numpoints=" << numlandmarks << std::endl;

  std::unique_ptr<bisSimpleMatrix<float> > points(new bisSimpleMatrix<float>("points_points_json"));
  if (!points->linkIntoPointer(in_points))
    return 0;

  if (debug)
    std::cout << "___ Ref Allocated = " << points->getNumRows() << "*" << points->getNumCols() << std::endl;
  
  std::unique_ptr<bisSimpleVector<int> > labels(new bisSimpleVector<int>("labels_points_json"));
  if (!labels->linkIntoPointer(in_labels))
    return 0;

  if (debug) 
    std::cout << "___ Labels Allocated = " << labels->getLength()  << std::endl;


  std::unique_ptr<bisSimpleMatrix<float> > out_points(new bisSimpleMatrix<float>("out_points_points_json"));
  std::unique_ptr<bisSimpleVector<int> > out_labels(new bisSimpleVector<int>("out_labels_points_json"));
  
  bisRPMCorrespondenceFinder::samplePoints(points.get(),
                                           labels.get(),
                                           numlandmarks,
                                           prefsampling,
                                           out_points.get(),
                                           out_labels.get(),
                                           debug);

  return out_points->releaseAndReturnRawArray();
}
