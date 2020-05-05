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

#include "bisRPMCorrespondenceFinder.h"
#include "bisPointLocator.h"
#include "bisPointRegistrationUtils.h"


bisRPMCorrespondenceFinder::bisRPMCorrespondenceFinder(std::string n) : bisObject(n) {

  this->locator=NULL;
  this->SampledReferencePoints=0;
  this->SampledTargetPoints=0;
  this->SampledReferenceLabels=NULL;
  this->SampledTargetLabels=NULL;
  this->class_name="bisRPMCorrespondenceFinder";
}

bisRPMCorrespondenceFinder::~bisRPMCorrespondenceFinder() {
  this->cleanup();
}

void bisRPMCorrespondenceFinder::cleanup() {
  
  if (this->locator)
    delete this->locator;

  this->SampledReferencePoints=0;
  this->SampledTargetPoints=0;
  
  if (this->SampledReferenceLabels)
    delete this->SampledReferenceLabels;
  
  if (this->SampledTargetLabels)
    delete this->SampledTargetLabels;

  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------
int bisRPMCorrespondenceFinder::samplePoints(bisSimpleMatrix<float>* inputPoints,
                                             bisSimpleVector<int>* inputLabels,
                                             int maxnumpoints,
                                             int samplingweight,
                                             bisSimpleMatrix<float>* outputPoints,
                                             bisSimpleVector<int>* outputLabels,
                                             int debug)
{
  samplingweight=bisUtil::irange(samplingweight,1,5);
  int uselabels=0;
  if (inputLabels)
    uselabels=1;

  int numrows=inputPoints->getNumRows();
  int numcols=inputPoints->getNumCols();
  int numrows2=numrows;
  if (uselabels) {
    numrows2=inputLabels->getLength();
  }

  if (numrows<4 || numcols!=3 || numrows2 != numrows) {
    std::cerr << "___ Failed to sample " << numrows << "*" << numcols << std::endl;
    if (uselabels)
      std::cerr << "___      uselabels=" << uselabels << " numrows=" << numrows2 << std::endl;
  }
  
  if (!uselabels || samplingweight<=1) {
    
    int step=1;
    if (maxnumpoints < numrows) {
      step=int(numrows/maxnumpoints);
    }
    int actualnumpoints=numrows/step;
    outputPoints->zero(actualnumpoints,3);
    outputLabels->zero(actualnumpoints);
    float* outpts=outputPoints->getData();

    if (debug) {
      std::cout << "___ Sampling from " << numrows << " points to " << actualnumpoints << " (step=" << step << ")" << std::endl;
      std::cout << "___ Actual Num Points=" << actualnumpoints << std::endl;
    }
    float* pts=inputPoints->getData();
    
    for (int i=0;i<actualnumpoints;i++) {
      int index=i*step*3;
      for (int ia=0;ia<=2;ia++)
        outpts[i*3+ia]=pts[index+ia];
    }
    return actualnumpoints;
  }

  // Preferential sampling ... harder
  int count[2]={0,0};
  float* pts=inputPoints->getData();
  int* labels=inputLabels->getData();
  
  for (int i=0;i<numrows;i++) {
    if (labels[i]>0)
      count[1]=count[1]+1;
  }
  count[0]=numrows-count[1];

  int nominaltotal=count[0]+count[1]*samplingweight;
  int step[2]={1,1};
  if (maxnumpoints < nominaltotal) {
    float fraction=float(count[0])/float(nominaltotal);
    step[0]=int(count[0]/(fraction*maxnumpoints));
    step[1]=int(count[1]/((1.0-fraction)*maxnumpoints));
  }

  
  int actual[2]= { count[0]/step[0], count[1]/step[1]};
  
  int actualnumpoints=actual[0]+actual[1];
  if (debug) {
    std::cout << "___ Pref Sampling from " << numrows << "(" << count[0] << "," << count[1] << ") points to " << actualnumpoints;
    std::cout << "(" << actual[0] << "," << actual[1] << ") ( step= " << step[0] << "," << step[1] << ")" << std::endl;
  }
  

  outputPoints->zero(actualnumpoints,3);
  outputLabels->zero(actualnumpoints);
  float* outpts=outputPoints->getData();
  int* outlabels=outputLabels->getData();

  int out_index=0;

  int added[2] = { 0,0 };
  
  for (int pass=0;pass<=1;pass++) {
    int count=0;
    for (int i=0;i<numrows;i++) {
      int add=0;
      if ( count==0 && ( ( labels[i]==0 && pass==0 ) || (labels[i]>0 && pass==1)))
        add=1;
      
      if (add) {
        for (int ia=0;ia<=2;ia++)
          outpts[out_index*3+ia]=pts[i*3+ia];
        outlabels[out_index]=labels[i];
        out_index=out_index+1;
        added[pass]=added[pass]+1;
        if (added[pass] >=  actual[pass])
          i=numrows;
      }
      if ( (labels[i]==0 && pass==0) || (labels[i]>0 && pass==1)) {
        count=count+1;
        if (count==step[pass])
          count=0;
      }
    }
  }

  return actualnumpoints;
}
                 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------
// Description:
// Specify the ref and target data sets.
int bisRPMCorrespondenceFinder::initialize(bisSimpleMatrix<float>* ref,
                                           bisSimpleMatrix<float>* target,
                                           int maxnumlandmarks,
                                           int samplingweight, // 1 = equal sampling, if  > 1 get extra points from points whose label > 0
                                           bisSimpleVector<int>* ref_labels,
                                           bisSimpleVector<int>* target_labels,
                                           int debug) {

  this->cleanup();
  int ok[2];
  
  std::shared_ptr<bisSimpleMatrix<float> > reft(new bisSimpleMatrix<float>());
  this->SampledReferencePoints=std::move(reft);

  std::shared_ptr<bisSimpleMatrix<float> > targt(new bisSimpleMatrix<float>());
  this->SampledTargetPoints=std::move(targt);
  
  this->SampledTargetLabels=new bisSimpleVector<int>();
  this->SampledReferenceLabels=new bisSimpleVector<int>();
  
  ok[0]=this->samplePoints(ref,ref_labels,maxnumlandmarks,samplingweight,this->SampledReferencePoints.get(),this->SampledReferenceLabels,debug);
  ok[1]=this->samplePoints(target,target_labels,int(maxnumlandmarks*1.5),samplingweight,this->SampledTargetPoints.get(),this->SampledTargetLabels,debug);

  if (ok[0] == 0 || ok[1] ==0) {
    std::cerr << "Failed to sample points " << std::endl;
    this->cleanup();
    return 0;
  }

  this->locator=new bisPointLocator();
  this->locator->initialize(this->SampledTargetPoints,0.1,0);
  return 1;
   
  
}
// ----------------------------------------------------------------------------------------------------------------------------------------------------------

int bisRPMCorrespondenceFinder::estimateCorrespondence(bisAbstractTransformation* Transformation,
                                                       float temperature,
                                                       int mode,
                                                       bisSimpleMatrix<float>* OutputRefLandmarks,
                                                       bisSimpleMatrix<float>* OutputTargetLandmarks,
                                                       bisSimpleVector<float>* OutputWeights,
                                                       int debug)
{
  int numref=this->SampledReferencePoints->getNumRows();
  int numtarget=this->SampledTargetPoints->getNumRows();
  if (OutputRefLandmarks->getNumRows()!=numref || OutputRefLandmarks->getNumCols()!=3) 
    OutputRefLandmarks->zero(numref,3);
  
  if (OutputTargetLandmarks->getNumRows()!=numref || OutputTargetLandmarks->getNumCols()!=3) 
    OutputTargetLandmarks->zero(numref,3);
  
  if (OutputWeights->getLength()!=numref)
    OutputWeights->zero(numref);
  
  float* reference_pts=this->SampledReferencePoints->getData();
  
  float* out_ref=OutputRefLandmarks->getData();
  float* out_target=OutputTargetLandmarks->getData();
  float* out_weights=OutputWeights->getData();
  
  if (mode==0) {
    return bisRPMCorrespondenceFinder::computeCorrespodnencesICP(Transformation,
                                                                 locator,
                                                                 reference_pts,
                                                                 out_ref,
                                                                 out_target,
                                                                 out_weights,
                                                                 numref,
                                                                 debug);
  }
  
  float* target_pts=this->SampledTargetPoints->getData();
  int* reference_labels=this->SampledReferenceLabels->getData();
  int* target_labels=this->SampledTargetLabels->getData();
  
  bisRPMCorrespondenceFinder::computeCorrespondencesRPM(Transformation,locator,
                                                        mode,
                                                        reference_pts,reference_labels,
                                                        target_pts,target_labels,
                                                        out_ref,out_target,out_weights,
                                                        temperature,numref,numtarget,debug);

  
  


  return 1;
}
// ----------------------------------------------------------------------------------------------------------------------------------------------------------

int bisRPMCorrespondenceFinder::computeCorrespodnencesICP(bisAbstractTransformation* Transformation,
                                                          bisPointLocator* locator,
                                                          float* reference_pts,
                                                          float* out_ref,
                                                          float* out_target,
                                                          float* out_weights,
                                                          int numref,
                                                          int debug) {

  if (debug)
    std::cout << "___ Computing ICP correspondences for " << numref << " points" << std::endl;

  int half=int(numref/2);
  
  for (int i=0;i<numref;i++) {
    float x[3],y[3],tx[3];
    for (int ia=0;ia<=2;ia++)
      x[ia]=reference_pts[i*3+ia];
    Transformation->transformPoint(x,tx);
    
    locator->getNearestPoint(tx,y,0);
    for (int ia=0;ia<=2;ia++) {
      out_ref[i*3+ia]=x[ia];
      out_target[i*3+ia]=y[ia];
      out_weights[i]=1.0;
    }

    if  ( (i==0 || i==half) && debug>0) {
        std::cout << "Point = " << i << " x=" << x[0] << "," << x[1] << " " << x[2] << " --> ";
        std::cout << "tx=" << tx[0] << "," << tx[1] << " " << tx[2] << " --> ";
        std::cout << "y=" << y[0] << "," << y[1] << " " << y[2] << std::endl;
    }
    
  }
  return 1;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------
// RPM
// ----------------------------------------------------------------------------------------------------------------------------------------------------------

int bisRPMCorrespondenceFinder::normalizeMatrixMixture(Eigen::SparseMatrix<float,Eigen::RowMajor> M) {

  // Compute Sums
  for (int k=0; k<M.outerSize(); ++k) {

    float sum=0.01;
    
    for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it ; ++it)
      sum+=it.value();

    //if (k==468 || k==117)
    //std::cout << "Sum for row=" << k << "=" << sum << std::endl;
    
    for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it ; ++it) 
      M.coeffRef(it.row(),it.col())=it.value()/sum;
  }

  return 1;
}

int bisRPMCorrespondenceFinder::normalizeMatrixRPM(Eigen::SparseMatrix<float,Eigen::RowMajor> M) {

  int numrefpts=M.rows();
  std::vector<float> outlier_ref(numrefpts);
  for (int i=0;i<numrefpts;i++) 
    outlier_ref[i]=0.01;

  int numtargetpts=M.cols();
  std::vector<float> outlier_targ(numtargetpts);
  for (int i=0;i<numtargetpts;i++)
    outlier_targ[i]=0.01;

  std::vector<float> sumcols(numtargetpts);
  for (int sinkolhm_iteration=0;sinkolhm_iteration<=4;sinkolhm_iteration++)
    {
      // Sum Rows
      // Compute Sums
      for (int k=0; k<M.outerSize(); ++k) {
        float sum=outlier_ref[k];
        
        for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it ; ++it)
          sum+=it.value();

        for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it; ++it)
          M.coeffRef(it.row(),it.col())=it.value()/sum;
        outlier_ref[k]/=sum;
      }

      // Sum Columns now (a little messier as this is the inner loop
      for (int column=0;column<numtargetpts;column++)
        sumcols[column]=outlier_targ[column];
      
      // Now Columns Sum
      for (int k=0; k<M.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it ; ++it)  {
          sumcols[it.col()]+=it.value();
        }
      }
      
      for (int column=0;column<numtargetpts;column++)
        outlier_targ[column]/=sumcols[column];

      for (int k=0; k<M.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it ; ++it) {
          M.coeffRef(it.row(),it.col())=it.value()/sumcols[it.col()];
        }
      }
    }
  return 1;
}
// ----------------------------------------------------------------------------------------------------------------------------------------------------------
int bisRPMCorrespondenceFinder::computeCorrespondencesRPM(bisAbstractTransformation* Transformation,bisPointLocator* locator,
                                                          int mode,
                                                          float* reference_pts,int* reference_labels,
                                                          float* target_pts,int* target_labels,
                                                          float* out_ref, float* out_target,float* out_weights,
                                                          float temperature,int numrefpts,int numtargetpts,int debug) {

  // Create extra columns and rows
  Eigen::SparseMatrix<float,Eigen::RowMajor> M(numrefpts,numtargetpts);
  Eigen::VectorXf outlier_ref=Eigen::VectorXf::Zero(numrefpts);
  Eigen::VectorXf sum_ref=Eigen::VectorXf::Zero(numrefpts);
  Eigen::VectorXf outlier_targ=Eigen::VectorXf::Zero(numtargetpts);
  Eigen::VectorXf sum_targ=Eigen::VectorXf::Zero(numtargetpts);

  // Initialize outliers
  for (int i=0;i<numrefpts;i++)
    outlier_ref(i)=0.01;
  for (int i=0;i<numtargetpts;i++)
    outlier_targ(i)=0.01;
  
  std::cout << "_____ Computing RPM Correspondences mode=" << mode << " Points=" << numrefpts << "*" << numtargetpts << " temperature=" << temperature << std::endl;
  float threshold=3*temperature;
  float T2=2.0*temperature*temperature;
  std::vector<int> pointlist;
  
  //std::cout << "___ nums=" << numrefpts << "*" << numtargetpts << " M=" <<  M.outerSize() << "*" << M.innerSize() << std::endl;
    
  for (int row=0;row<numrefpts;row++)
    {
      float x[3],tx[3];
      for (int ia=0;ia<=2;ia++) {
        x[ia]=reference_pts[row*3+ia];
        out_ref[row*3+ia]=x[ia];
      }

      Transformation->transformPoint(x,tx);
      int nump=locator->getPointsWithinRadius(tx,threshold,pointlist,0);
      if (nump<3) {
        pointlist.clear();
        locator->getPointsWithinRadius(tx,threshold*2.0,pointlist,0);
      }

      float sum=0.0;

      //if (row==468 || row==117)  {
      //        std::cout << "___ row=" << row << " nump=" << nump << std::endl;
      //      }
      
      for (int i=0;i<nump;i++) {
        int col=pointlist[i];
        if (reference_labels[row] == target_labels[col]) {
          float y[3];
          for (int ia=0;ia<=2;ia++) 
            y[ia]=target_pts[col*3+ia];

          float dist2=bisPointRegistrationUtils::distance2(tx,y);
          float d=exp(-dist2/T2);
          sum+=d;
          M.coeffRef(row,col)=d;
          /*if (row==468 || row==117)  {
            std::cout << "row=" << row << "," << col << " (nump=" << nump << ")=" << d << "(dist2=" << dist2 << "," << T2 << ") sum=" << sum << " x=" << x[0] << "," << x[1] << "," << x[2];
            std::cout << " -> tx=" << tx[0] << "," << tx[1] << "," << tx[2];
            std::cout << " -> y=" << y[0] << "," << y[1] << "," << y[2] << std::endl;
            }*/
          
        }
      }

      // If total sum is too low < 0.001, then stablize using closest point and set its weight to 0.001
      if (sum<0.1) {
        float y[3];
        int nearest=locator->getNearestPoint(tx,y,0);
        M.coeffRef(row,nearest)=0.001;
        sum=0.001;
        /*if (row==468 || row==117)  {
          std::cout << "SMALL SUM row=" << row << "," << nearest << " sum=" << sum << " x=" << x[0] << "," << x[1] << "," << x[2];
          std::cout << " -> tx=" << tx[0] << "," << tx[1] << "," << tx[2];
          std::cout << " -> y=" << y[0] << "," << y[1] << "," << y[2] << std::endl;
          }*/
      }
    }

  if (mode==1)  {
    if (debug)
      std::cout << "__ Normalizing Matrix as Mixture" << std::endl;
    bisRPMCorrespondenceFinder::normalizeMatrixMixture(M);
  } else {
    std::cout << "__ Normalizing Matrix as Full RPM" << std::endl;
    bisRPMCorrespondenceFinder::normalizeMatrixRPM(M);
  }
  
  for (int i=0;i<numrefpts;i++) {
    out_weights[i]=0.0;
    out_target[i*3]=0.0;
    out_target[i*3+1]=0.0;
    out_target[i*3+2]=0.0;
  }
    
  for (int k=0; k<M.outerSize(); ++k) {
    for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it; ++it) {
      float v=it.value();
      int row=it.row();   // row index
      int col=it.col();   // col index (here it is equal to k)
      out_weights[row]+=v;
      for (int ia=0;ia<=2;ia++) 
        out_target[row*3+ia]+=v*target_pts[col*3+ia];
    }
  }

  /*if (debug) {
    int row=117;
    std::cout << "weighted_out_target=  for row=" << row << "= " 
              <<  out_target[row*3] << "," 
              <<  out_target[row*3+1] << "," 
              <<  out_target[row*3+2] 
              <<  " wgt=" << out_weights[row] << std::endl;
    row=468;
    std::cout << "weighted_out_target=  for row=" << row << "= " 
              <<  out_target[row*3] << "," 
              <<  out_target[row*3+1] << "," 
              <<  out_target[row*3+2] 
              <<  " wgt=" << out_weights[row] << std::endl;
              }*/
      
  for (int i=0;i<numrefpts;i++) {
    for (int ia=0;ia<=2;ia++)  {
      out_target[i*3+ia]/=out_weights[i];
    }
  }

  /*if (debug) {
    int row=117;
    std::cout << "out_target=  for row=" << row << "= "
              <<  out_target[row*3] << "," 
              <<  out_target[row*3+1] << "," 
              <<  out_target[row*3+2] 
              <<  " wgt=" << out_weights[row] << std::endl;
    row=468;
    std::cout << "out_target=  for row=" << row << "= "
              <<  out_target[row*3] << "," 
              <<  out_target[row*3+1] << "," 
              <<  out_target[row*3+2] 
              <<  " wgt=" << out_weights[row] << std::endl;
              }*/

  
  return 1;
}

                 
