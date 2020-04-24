/*  LICENSE
 
 _This file is Copyright 2018 by the Image Processing and Analysis Group (BioImage Suite Team). Dept. of Radiology & Biomedical Imaging, Yale School of Medicine._
 
 BioImage Suite Web is licensed under the Apache License, Version 2.0 (the "License");
 
 - you may not use this software except in compliance with the License.
 - You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)
 
 __Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.__
 
 ENDLICENSE */

// While the rest of the files in this directory are GPL v2, this particular file is released under Apache V2 (which can be merged with GPL code)


#ifndef _bis_GPLExportedFunctions_h
#define _bis_GPLExportedFunctions_h

#include "bisDefinitions.h"


#ifdef __cplusplus
extern "C" {
#endif


  /** @namespace bisGPLExportedFunctions
      Functions exported to JS and Python. 
      See \link bisExportedFunctions.h \endlink and secondarily
      \link bisTesting.h \endlink.
  */
  
  
  /** @file bisGPLExportedFunctions.h
      Functions exported to JS and Python
  */

  /** Returns 1*/
  // BIS: { 'uses_gpl', 'Int' } 
  BISEXPORT int uses_gpl();
  
  
  /** run Linear Image Registration using \link bisLinearImageRegistration  \endlink
   * @param reference serialized reference image as unsigned char array 
   * @param target    serialized target image as unsigned char array 
   * @param initial_xform serialized initial transformation as unsigned char array 
   * @param jsonstring the parameter string for the algorithm including return_vector which if true returns a length-28 vector
   * containing the 4x4 matrix and the 12 transformation parameters
   * @param debug if > 0 print debug messages
   * @returns a pointer to a serialized vector or matrix depending on the value of return_vector
   */
  // BIS: { 'runLinearRegistrationWASM', 'bisLinearTransformation', [ 'bisImage', 'bisImage', 'bisLinearTransformation_opt', 'ParamObj', 'debug' ], {"checkorientation" : "python matlab"} } 
  BISEXPORT  unsigned char*  runLinearRegistrationWASM(unsigned char* reference,
						       unsigned char* target,
						       unsigned char* initial_xform,
						       const char* jsonstring,
						       int debug);
  
  
  /** run Non Linear Image Registration using \link bisNonLinearImageRegistration  \endlink
   * @param reference serialized reference image as unsigned char array 
   * @param target    serialized target image as unsigned char array 
   * @param initial_xform serialized initial transformation as unsigned char array 
   * @param jsonstring the parameter string for the algorithm 
   * @param debug if > 0 print debug messages
   * @returns a pointer to a serialized combo transformation (bisComboTransformation)
   */
  // BIS: { 'runNonLinearRegistrationWASM', 'bisComboTransformation', [ 'bisImage', 'bisImage', 'bisLinearTransformation_opt', 'ParamObj', 'debug' ], {"checkorientation" : "python matlab"}  } 
  BISEXPORT unsigned char* runNonLinearRegistrationWASM(unsigned char* reference,
							unsigned char* target,
							unsigned char* initial_xform,
							const char* jsonstring,
							int debug);



  /** Approximate Displacement Field with Grid Transformation (pre initialized)
   * @param dispfield serialized target displacement field
   * @param initial_grid serialized grid transformation as unsigned char array 
   * @param jsonstring the parameter string for the algorithm 
   * @param debug if > 0 print debug messages
   * @returns a pointer to the updated grid (bisGridTransformation)
   */
  // BIS: { 'approximateDisplacementFieldWASM', 'bisGridTransformation', [ 'bisImage', 'bisGridTransformation', 'ParamObj', 'debug' ] } 
  BISEXPORT unsigned char* approximateDisplacementFieldWASM(unsigned char* dispfield,
							    unsigned char* initial_grid,
							    const char* jsonstring,
							    int debug);
  
  /** Approximate Displacement Field with Grid Transformation -- initialized using the sapcing parameter
   * @param dispfield serialized target displacement field
   * @param jsonstring the parameter string for the algorithm  -- key is spacing : --> this defines the spacing for the grid transformation
   * @param debug if > 0 print debug messages
   * @returns a pointer to the updated grid (bisGridTransformation)
   */
  // BIS: { 'approximateDisplacementFieldWASM2', 'bisGridTransformation', [ 'bisImage', 'ParamObj', 'debug' ] } 
  BISEXPORT unsigned char* approximateDisplacementFieldWASM2(unsigned char* dispfield,
							     const char* jsonstring,
							     int debug);
  

  /** Perform image segmentation either histogram based or plus mrf segmentation if smoothness > 0.0
   * @param input serialized input as unsigned char array 
   * @param jsonstring the parameter string for the algorithm { "numclasses" : 3, "maxsigmaratio":0.2, "robust" : true, "numbins": 256, "smoothhisto": true, "smoothness" : 0.0, "mrfconvergence" : 0.2, "mrfiterations" : 8, "noisesigma2" : 0.0 }
   * @param debug if > 0 print debug messages
   * @returns a pointer to a serialized segmented image 
   */
  // BIS: { 'segmentImageWASM', 'bisImage', [ 'bisImage', 'ParamObj', 'debug' ] } 
  BISEXPORT unsigned char* segmentImageWASM(unsigned char* input,const char* jsonstring,int debug);

  /** Perform objectmap regularization 
   * @param input serialized input as unsigned char array 
   * @param jsonstring the parameter string for the algorithm { "smoothness" : 2.0, "convergence" : 0.2, "terations" : 8, "internaliterations" : 4 }
   * @param debug if > 0 print debug messages
   * @returns a pointer to a (short) serialized segmented image 
   */
  // BIS: { 'regularizeObjectmapWASM', 'bisImage', [ 'bisImage', 'ParamObj', 'debug' ] } 
  BISEXPORT unsigned char* regularizeObjectmapWASM(unsigned char* input,const char* jsonstring,int debug);

  
  /** Tests Optimizer with numdof = 1 or 2 and all three modes 
   * @param numdof number of degrees of freedom for simple quadratic function (1 or 2)
   * @returns number of failed tests
   */
  // BIS: { 'test_optimizer', 'Int', [ 'Int'] } 
  BISEXPORT int test_optimizer(int numdof);


  /** Compute DTI Tensor
   * @param input_ptr the images as a serialized array
   * @param baseline_ptr the "Baseline" T2 Image as a serialized array
   * @param mask_ptr the Mask Image (optional, set this to 0) as a serialized array
   * @param directions_ptr the directions matrix
   * @param jsonstring { "bvalue": 1000, "numbaseline:" 1 }
   * @param debug if > 0 print debug messages
   * @returns a pointer to the tensor image */
  // BIS: { 'computeDTITensorFitWASM', 'bisImage', [ 'bisImage', 'bisImage',  'bisImage_opt' ,'Matrix', 'ParamObj', 'debug'] } 
  BISEXPORT unsigned char* computeDTITensorFitWASM(unsigned char* input_ptr,
						   unsigned char* baseline_ptr,
						   unsigned char* mask_ptr,
						   unsigned char* directions_ptr,
						   const char* jsonstring,
						   int debug);


  /** Compute DTI Tensor EigenSystem
   * @param input_ptr the image tensor as a serialized array
   * @param mask_ptr the Mask Image (optional, set this to 0) as a serialized array
   * @param debug if > 0 print debug messages
   * @returns a pointer to the eigensystem image */
  // BIS: { 'computeTensorEigenSystemWASM', 'bisImage', [ 'bisImage', 'bisImage_opt' , 'debug'] } 
  unsigned char* computeTensorEigenSystemWASM(unsigned char* input_ptr,
					      unsigned char* mask_ptr,
					      int debug);


  /** Compute DTI Tensor Invariants
   * @param input_ptr the image tensor eigensystem as a serialized array
   * @param mask_ptr the Mask Image (optional, set this to 0) as a serialized array
   * @param jsonstring { "mode": 0 } // mode 0=FA, 1=RA etc. -- see bisDTIAlgorithms::computeTensorInvariants
   * @param debug if > 0 print debug messages
   * @returns a pointer to the invarient image */
  // BIS: { 'computeDTITensorInvariantsWASM', 'bisImage', [ 'bisImage', 'bisImage_opt' , 'ParamObj', 'debug'] } 
  BISEXPORT unsigned char* computeDTITensorInvariantsWASM(unsigned char* input_ptr,
							  unsigned char* mask_ptr,
							  const char* jsonstring,
							  int debug);

  /** Compute DTI Orientation Map
   * @param input_ptr the image tensor eigensystem as a serialized array
   * @param mask_ptr the Mask Image (optional, set this to 0) as a serialized array
   * @param magnitude_ptr the Magnitude Image (e.g. FA map) (optional, set this to 0) as a serialized array
   * @param jsonstring { "scaling": 1.0 } Optional extra scaling
   * @param debug if > 0 print debug messages
   * @returns a pointer to the colormap image */
  // BIS: { 'computeDTIColorMapImageWASM', 'bisImage', [ 'bisImage', 'bisImage_opt' ,'bisImage_opt', 'ParamObj', 'debug'] } 
  BISEXPORT unsigned char* computeDTIColorMapImageWASM(unsigned char* input_ptr,
						       unsigned char* mask_ptr,
						       unsigned char* magnitude_ptr,
						       const char* jsonstring,
						       int debug);



    /** runWeighted Linear Image Registration using \link bisLinearImageRegistration  \endlink
   * @param reference serialized reference image as unsigned char array 
   * @param target    serialized target image as unsigned char array 
   * @param ref_weight  serialized reference weight image as unsigned char array 
   * @param targ_weight  serialized target weight image as unsigned char array 
   * @param initial_xform serialized initial transformation as unsigned char array 
   * @param jsonstring the parameter string for the algorithm including return_vector which if true returns a length-28 vector
   * containing the 4x4 matrix and the 12 transformation parameters
   * @param debug if > 0 print debug messages
   * @returns a pointer to a serialized vector or matrix depending on the value of return_vector
   */
  // BIS: { 'runWeightedLinearRegistrationWASM', 'bisLinearTransformation', [ 'bisImage', 'bisImage', 'bisImage', 'bisImage_opt', 'bisLinearTransformation_opt', 'ParamObj', 'debug' ], {"checkorientation" : "python matlab"} } 
  BISEXPORT  unsigned char*  runWeightedLinearRegistrationWASM(unsigned char* reference,
                                                               unsigned char* target,
                                                               unsigned char* reference_weight,
                                                               unsigned char* target_weight,
                                                               unsigned char* initial_xform,
                                                               const char* jsonstring,
                                                               int debug);
  
  
  /** runWeighted Non Linear Image Registration using \link bisNonLinearImageRegistration  \endlink
   * @param reference serialized reference image as unsigned char array 
   * @param target    serialized target image as unsigned char array 
   * @param ref_weight  serialized reference weight image as unsigned char array 
   * @param targ_weight  serialized target weight image as unsigned char array 
   * @param initial_xform serialized initial transformation as unsigned char array 
   * @param jsonstring the parameter string for the algorithm 
   * @param debug if > 0 print debug messages
   * @returns a pointer to a serialized combo transformation (bisComboTransformation)
   */
  // BIS: { 'runWeightedNonLinearRegistrationWASM', 'bisComboTransformation', [ 'bisImage', 'bisImage', 'bisImage', 'bisImage_opt', 'bisLinearTransformation_opt', 'ParamObj', 'debug' ], {"checkorientation" : "python matlab"}  } 
  BISEXPORT unsigned char* runWeightedNonLinearRegistrationWASM(unsigned char* reference,
                                                                unsigned char* target,
                                                                unsigned char* reference_weight,
                                                                unsigned char* target_weight,
                                                                unsigned char* initial_xform,
                                                                const char* jsonstring,
                                                                int debug);


  /** Tests Optimizer with numdof = 1 or 2 and all three modes 
   * @param reference serialized reference points as  unsigned char array 
   * @param target    serialized target points as unsigned char array 
   * @param spacing   grid_spacing for transformation
   * @returns a pointer to the updated grid (bisGridTransformation)
   */
  // BIS: { 'test_landmarkApproximationWASM', 'bisGridTransformation', [ 'Matrix', 'Matrix', 'ParamObj', debug] } 
  BISEXPORT unsigned char*  test_landmarkApproximationWASM(unsigned char* in_source,
                                                           unsigned char* in_target,
                                                           const char* jsonstring,
                                                           int debug);

  
#ifdef __cplusplus
}
#endif

#endif
