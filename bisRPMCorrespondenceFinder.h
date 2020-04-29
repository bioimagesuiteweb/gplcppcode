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



#ifndef __bisRPMCorrespondenceFinder_h
#define __bisRPMCorrespondenceFinder_h

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "bisAbstractTransformation.h"
#include "bisSimpleDataStructures.h"

class bisPointLocator;

class bisRPMCorrespondenceFinder : public bisObject  {

public:

  bisRPMCorrespondenceFinder(std::string n);
  virtual ~bisRPMCorrespondenceFinder();

  // Description:
  // Specify the source and target data sets.
  int initialize(bisSimpleMatrix<float>* source,
                 bisSimpleMatrix<float>* target,
                 int maxnumlandmarks=1000,
                 int samplingweight=1, // 1 = equal sampling, if  > 1 get extra points from points whose label > 0
                 bisSimpleVector<int>* source_labels=0,
                 bisSimpleVector<int>* target_labels=0,
                 int debug=0);
  

  // Description:
  // Estimate Correspondences and Weights
  int estimateCorrespondence(bisAbstractTransformation* Transformation,
                             float temperature,
                             int mode,
                             bisSimpleMatrix<float>* OutputSourceLandmarks,
                             bisSimpleMatrix<float>* OutputTargetLandmarks,
                             bisSimpleMatrix<float>* OutputWeights,
                            int debug=0);

  // Sample Points
  static int samplePoints(bisSimpleMatrix<float>* sourcePoints,
                          bisSimpleVector<int>* sourceLabels,
                          int maxnumpoints,
                          int samplingweight,
                          bisSimpleMatrix<float>* outputPoints,
                          bisSimpleVector<int>* outputLabels,
                          int debug=0);


  /** Compute ICP Correspondence */
  static int computeCorrespodnencesICP(bisAbstractTransformation* Transformation,bisPointLocator* locator,
                                       float* reference_pts,float* out_ref,float* out_target,float* out_weights,int numref,int debug=0);

  /** Compute And Normalize Distance Matrix */
  static int computeCorrespondencesRPM(bisAbstractTransformation* Transformation,bisPointLocator* locator,
                                       Eigen::SparseMatrix<float,Eigen::RowMajor>& M,
                                       int mode,
                                       float* reference_pts,int* reference_labels,
                                       float* target_pts,int* target_labels,
                                       float* out_ref,float* out_target,float* out_weights,
                                       float temperature,int numrefpts,int numtargetpts,int debug=0);


 protected:
  std::shared_ptr<bisSimpleMatrix<float> > SampledReferencePoints;
  bisSimpleMatrix<float>* SampledTargetPoints;
  bisSimpleVector<int>* SampledReferenceLabels;
  bisSimpleVector<int>* SampledTargetLabels;
  int numpoints_ref;
  int numpoints_target;
  bisPointLocator* locator;

  // Cleanup
  void cleanup();

private:

  /** Copy constructor disabled to maintain shared/unique ptr safety */
  bisRPMCorrespondenceFinder(const bisRPMCorrespondenceFinder&);

  /** Assignment disabled to maintain shared/unique ptr safety */
  void operator=(const bisRPMCorrespondenceFinder&);  
                             
  
};


#endif
