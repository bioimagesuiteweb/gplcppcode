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


#ifndef __bisNonLinearRPMRegistration_h
#define __bisNonLinearRPMRegistration_h

#include "bisComboTransformation.h"
#include "bisPointRegistrationUtils.h"
#include "bisRPMCorrespondenceFinder.h"

class bisNonLinearRPMRegistration : public bisRPMCorrespondenceFinder
{
public:

  bisNonLinearRPMRegistration(std::string n="rpm_linear");
  virtual ~bisNonLinearRPMRegistration();


  // Description:
  // get Output Transformation  
  std::shared_ptr<bisComboTransformation> getOutput();


  // Description:
  // Specify the source and target data sets and the initial linear transformation
  int initializeWithLinear(bisMatrixTransformation* in_initialTransformation,
                           bisSimpleMatrix<float>* source,
                           bisSimpleMatrix<float>* target,
                           int maxnumlandmarks=3000,
                           int samplingweight=1, // 1 = equal sampling, if  > 1 get extra points from points whose label > 0
                           bisSimpleVector<int>* source_labels=0,
                           bisSimpleVector<int>* target_labels=0,
                           int debug=0);

  /** run transformation
   * @param correspondenceMode (0=ICP,1=Mixture RPM,2=Full RPM)
   * @param in_cps_begin initial control point spacing
   * @param in_cps_end final control point spacing
   * @param in_smoothness_begin initial smoothness
   * @param in_smoothness_end final smoothness
   * @param initialTemperature initial temperature
   * @param finalTemperature final temperature
   * @param iterationPerTemperature number of iterations at each temperature
   * @param annealRate anneal rate for RPM
   * @param debug if true print extra messages
   * @returns number of iterations if  success (>=1), or 0 if failed
   */
  virtual int run(int correspondenceMode=2,
                  float in_cps_begin=40.0,
                  float in_cps_end=20.0,
                  float in_smoothness_begin=1.0,
                  float in_smoothness_end=0.1,
                  float initialTemperature=5.0,
                  float finalTemperature=1.0,
                  int iterationPerTemperature=5,
                  float annealRate=0.93,
                  int debug=0);
  
protected:
  
  /** Description */
  std::shared_ptr<bisComboTransformation> Output;
  std::shared_ptr<bisGridTransformation> Grid;

  /** Initialize Grid */
  int initializeGrid(float cps,float cps_begin,float cps_end,bisSimpleMatrix<float>* RefLandmarks);

  /** approximate Grid */
  int approximateGrid(bisSimpleMatrix<float>* sourcePts,
                      bisSimpleMatrix<float>* targetPts,
                      bisSimpleVector<float>* weights,
                      float Smoothness=0.001,
                      int accurate=0,int debug=0);

  
private:

  /** Copy constructor disabled to maintain shared/unique ptr safety */
  bisNonLinearRPMRegistration(const bisNonLinearRPMRegistration&);

  /** Assignment disabled to maintain shared/unique ptr safety */
  void operator=(const bisNonLinearRPMRegistration&);  


};

#endif
