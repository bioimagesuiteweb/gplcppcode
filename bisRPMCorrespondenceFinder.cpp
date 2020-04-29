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

bisRPMCorrespondenceFinder::bisRPMCorrespondenceFinder(std::string n) : bisObject(n) {

  this->locator=NULL;
  this->SampledTargetPoints=NULL;
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

  if (this->SampledTargetPoints)
    delete this->SampledTargetPoints;
  
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
  int count[2];
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
  
  this->SampledReferenceLabels=new bisSimpleVector<int>();
  this->SampledTargetPoints=new bisSimpleMatrix<float>();
  this->SampledTargetLabels=new bisSimpleVector<int>();
  
  ok[0]=this->samplePoints(ref,ref_labels,maxnumlandmarks,samplingweight,this->SampledReferencePoints.get(),this->SampledReferenceLabels,debug);
  ok[1]=this->samplePoints(target,target_labels,maxnumlandmarks,samplingweight,this->SampledTargetPoints,this->SampledTargetLabels,debug);

  if (ok[0] == 0 || ok[1] ==0) {
    std::cerr << "Failed to sample points " << std::endl;
    this->cleanup();
    return 0;
  }

  this->locator=new bisPointLocator();
  this->locator->initialize(this->SampledReferencePoints,0.1,0);
  return 1;
   
  
}
// ----------------------------------------------------------------------------------------------------------------------------------------------------------

int bisRPMCorrespondenceFinder::estimateCorrespondence(bisAbstractTransformation* Transformation,
                                                       float temperature,
                                                       int mode,
                                                       bisSimpleMatrix<float>* OutputRefLandmarks,
                                                       bisSimpleMatrix<float>* OutputTargetLandmarks,
                                                       bisSimpleMatrix<float>* OutputWeights,
                                                       int debug)
{
  std::cout << "Debug = " << debug << std::endl;
  
  int numref=this->SampledReferencePoints->getNumRows();
  int numtarget=this->SampledTargetPoints->getNumRows();
  if (OutputRefLandmarks->getNumRows()!=numref || OutputRefLandmarks->getNumCols()!=3) 
    OutputRefLandmarks->zero(numref,3);
  else
    OutputRefLandmarks->fill(0.0);
  
  if (OutputTargetLandmarks->getNumRows()!=numref || OutputTargetLandmarks->getNumCols()!=3) 
    OutputTargetLandmarks->zero(numref,3);
  else
    OutputTargetLandmarks->fill(0.0);
  
  if (OutputWeights->getNumRows()!=numref)
    OutputWeights->zero(numref,1);
  else
    OutputTargetLandmarks->fill(0.0);

  float* reference_pts=this->SampledReferencePoints->getData();

  float* out_ref=OutputRefLandmarks->getData();
  float* out_target=OutputTargetLandmarks->getData();
  float* out_weights=OutputWeights->getData();
  
  if (mode==0) 
    return bisRPMCorrespondenceFinder::computeCorrespodnencesICP(Transformation,locator,reference_pts,out_ref,out_target,out_weights,numref,debug);

  float* target_pts=this->SampledTargetPoints->getData();
  int* reference_labels=this->SampledReferenceLabels->getData();
  int* target_labels=this->SampledTargetLabels->getData();
  
  Eigen::SparseMatrix<float,Eigen::RowMajor> M(numref,numtarget);
  bisRPMCorrespondenceFinder::computeCorrespondencesRPM(Transformation,locator,
                                                        M,mode,
                                                        reference_pts,reference_labels,
                                                        target_pts,target_labels,
                                                        out_ref,out_target,out_weights,
                                                        temperature,numref,numtarget,debug);


  return 1;
}
// ----------------------------------------------------------------------------------------------------------------------------------------------------------

int bisRPMCorrespondenceFinder::computeCorrespodnencesICP(bisAbstractTransformation* Transformation,bisPointLocator* locator,
                                                          float* reference_pts,float* out_ref,float* out_target,float* out_weights,int numref,int debug) {

  if (debug)
    std::cout << "___ Computing ICP correspondences for " << numref << " points" << std::endl;
  
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
  }
  return 1;
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------
int bisRPMCorrespondenceFinder::computeCorrespondencesRPM(bisAbstractTransformation* Transformation,bisPointLocator* locator,
                                                          Eigen::SparseMatrix<float,Eigen::RowMajor>& M,
                                                          int mode,
                                                          float* reference_pts,int* reference_labels,
                                                          float* target_pts,int* target_labels,
                                                          float* out_ref, float* out_target,float* out_weights,
                                                          float temperature,int numrefpts,int numtargetpts,int debug) {

  // Create extra columns and rows
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
  
  std::cout << "___ nums=" << numrefpts << "*" << numtargetpts << " M=" <<  M.outerSize() << "*" << M.innerSize() << std::endl;
    
  for (int row=0;row<numrefpts;row++)
    {
      float x[3],tx[3];
      for (int ia=0;ia<=2;ia++) {
	x[ia]=reference_pts[row*3+ia];
        out_ref[row*3+ia]=x[ia];
      }
      Transformation->transformPoint(x,tx);

      int db=0;
      if (row == 7 && debug && mode==1)
        db=1;
      
      int nump=locator->getPointsWithinRadius(tx,threshold,pointlist,db);
      for (int i=0;i<nump;i++) {
        int col=pointlist[i];
        if (reference_labels[row] == target_labels[col]) {
          float dist2=0.0;
          for (int ia=0;ia<=2;ia++) 
            dist2+=pow(x[ia]-target_pts[col*3+ia],2.0f);
          M.coeffRef(row,col)=exp(-dist2/T2);
          if (  row== 7 && debug)  {
            std::cout << "M(" << row << "," << col << ")=" << M.coeffRef(row,col) << " dist2=" << dist2 << " threshold=" << threshold << std::endl;
            std::cout << "\t\t " << x[0] << "," << x[1] << "," << x[2] << "-->";
            std::cout << "\t\t " << tx[0] << "," << tx[1] << "," << tx[2] << "-->";
            std::cout << "\t\t " << target_pts[col*3] << "," << target_pts[col*3+1] << "," << target_pts[col*3+2] << std::endl;
          }
        }
      }
    }
    
   
  int numiterations=4;
  if (mode==1)
    numiterations=1;

  for (int iteration=0;iteration<numiterations;iteration++) {

    if (debug) {
      std::cout << "____ iteration= " << iteration << "/" << numiterations << std::endl;
    }

    // Initialize Sums with Outlier Column/Rows
    for (int i=0;i<numrefpts;i++)
      sum_ref(i)=outlier_ref(i);
    for (int i=0;i<numtargetpts;i++)
      sum_targ(i)=outlier_targ(i);
    

    if (debug) {
      std::cout << "____ computing sums " << std::endl;
      std::cout << "____ initial = " << outlier_ref(7) << " " << outlier_targ(7) << "\t (" << sum_ref(7) << "," << sum_targ(7) << ")" << std::endl;
    }

    // Compute Sums
    for (int k=0; k<M.outerSize(); ++k) {
      for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it ; ++it)
        {
          float v=it.value();
          int row=it.row();   // row index
          int col=it.col();   // col index (here it is equal to k)
          if ( (row==7 || col==7) && debug)
            std::cout << "Row=" << row << " Col=" << col << " val=" << v << " k=" << k << std::endl;
          
          sum_ref(row)+=v;
          sum_targ(col)+=v;
        }
    }

    if (debug)  {
      std::cout << std::endl << "_______ Row=" << 7 << " sum=" << sum_ref(7) << " col7=" << sum_targ(7) << std::endl;
      std::cout << "____ normalizing rows or columns " << std::endl;
    }

    // Normalize either row or column
    for (int k=0; k<M.outerSize(); ++k) {
      Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k);
      for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it; ++it)
        {
          float v=it.value();
          int row=it.row();   // row index
          int col=it.col();   // col index (here it is equal to k)
          
          if ( (row==7 || col==7) && debug)
            std::cout << "row=" << row << " col=" << col << " v=" << v << std::endl;
          
          if (iteration % 2 == 0) 
            M.coeffRef(row,col)=v/sum_ref(row);
          else
            M.coeffRef(row,col)=v/sum_targ(col);
        }
    }

    // Normalize outlier column or row
    if (iteration % 2 == 0) {
      for (int i=0;i<numrefpts;i++) {
        if (debug && i==7)
          std::cout << "\t Before normalization outlier=" << outlier_ref(i) << " -->";
        outlier_ref(i)=outlier_ref(i)/sum_ref(i);
        if (debug && i==7)
          std::cout << " (after)" << outlier_ref(i) << " (" << sum_ref(i) << ")" << std::endl;
      }
    } else {
      for (int i=0;i<numtargetpts;i++)
        outlier_targ(i)=outlier_targ(i)/sum_targ(i);
    }
  }

  // Compute Correspondences
  if (debug)
    std::cout << "___ normalizing correspondences " << std::endl;
  
  for (int i=0;i<numrefpts;i++)
    out_weights[i]=0.000001;

  for (int k=0; k<M.outerSize(); ++k) {
    for (Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(M,k); it; ++it)
      {
        float v=it.value();
        int row=it.row();   // row index
        int col=it.col();   // col index (here it is equal to k)
        out_weights[row]+=v;
        for (int ia=0;ia<=2;ia++) 
          out_target[row*3+ia]+=v*target_pts[col*3+ia];

        if (row==7 && debug) {
          std::cout << "Row=" << row << " Col=" << col << " val=" << v << " out_weights[row]=" << out_weights[row] << " k=" << k << std::endl;
          std::cout << "\t Out_target=" << out_target[row*3] << "," << out_target[row*3+1] << "," << out_target[row*3+2] << std::endl;
        }
      }
  }
  
  for (int i=0;i<numrefpts;i++) {
    for (int ia=0;ia<=2;ia++)  {
      out_target[i*3+ia]/=out_weights[i];
    }
     if (i==7 && debug) {
       std::cout << "Row=" << i << " w=" << out_weights[i] << " pt=" << out_target[i*3] << "," << out_target[i*3+1] << "," << out_target[i*3+2] << std::endl;
     }

   
  }

  return 1;
}

                 
