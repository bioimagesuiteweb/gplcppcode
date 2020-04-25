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

  
  /** Tests Optimizer with numdof = 1 or 2 and all three modes 
   * @param numdof number of degrees of freedom for simple quadratic function (1 or 2)
   * @returns number of failed tests
   */
  // BIS: { 'test_optimizer', 'Int', [ 'Int'] } 
  BISEXPORT int test_optimizer(int numdof);


  /** Tests Landmark Approximation BSpline Code
   * @param reference serialized reference points as  unsigned char array 
   * @param target    serialized target points as unsigned char array 
   * @param spacing   grid_spacing for transformation
   * @param debug  debug flag
   * @returns a pointer to the updated grid (bisGridTransformation)
   */
  // BIS: { 'test_landmarkApproximationWASM', 'Matrix', [ 'Matrix', 'Matrix', 'ParamObj', debug] } 
  BISEXPORT unsigned char*  test_landmarkApproximationWASM(unsigned char* in_source,
                                                           unsigned char* in_target,
                                                           const char* jsonstring,
                                                           int debug);

  /** Tests RPM Correspondence Finder Code
   * @param reference serialized reference points as  unsigned char array 
   * @param target    serialized target points as unsigned char array 
   * @param ParamObj  JSON string (mode  0=icp 1=mixture, 2=full rpm, temperature, numlandmarks)
   * @param debug  debug flag
   * @returns a pointer to a matrix of either ICP closest or match matrix
   */
  // BIS: { 'test_rpmCorrespondenceEstimatorWASM', 'Matrix', [ 'Matrix', 'Matrix', 'ParamObj', debug] } 
  BISEXPORT unsigned char*  test_rpmCorrespondenceEstimatorWASM(unsigned char* in_source,
                                                                unsigned char* in_target,
                                                                const char* jsonstring,
                                                                int debug);


  /** Tests RPM Correspondence Sampling Code
   * @param points serialized points as  unsigned char array 
   * @param labels serialized label points as unsigned char array 
   * @param ParamObj  JSON string (numlandmakrs=100, preferentialsampling=1)
   * @param debug  debug flag
   * @returns a pointer to the sampled point matrix [ N x 3]
   */
  // BIS: { 'test_rpmSamplingWASM', 'Matrix', [ 'Matrix', 'Vector', 'ParamObj', debug] } 
  BISEXPORT unsigned char*  test_rpmSamplingWASM(unsigned char* in_source,
                                                 unsigned char* in_target,
                                                 const char* jsonstring,
                                                 int debug);
  
#ifdef __cplusplus
}
#endif

#endif
