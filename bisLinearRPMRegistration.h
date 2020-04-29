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


#ifndef __bisLinearRPMRegistration_h
#define __bisLinearRPMRegistration_h

#include "bisMatrixTransformation.h"
#include "bisPointRegistrationUtils.h"
#include "bisRPMCorrespondenceFinder.h"

class bisLinearRPMRegistration : public bisRPMCorrespondenceFinder
{
public:

  bisLinearRPMRegistration(std::string n="rpm_linear");
  virtual ~bisLinearRPMRegistration();


  // Description:
  // get Output Transformation  
  bisSimpleMatrix<float>* getOutputMatrix();

  /** run transformation
   * @param transformMode (0=rigid,1=similarity,2=affine)
   * @param correspondenceMode (0=ICP,1=Mixture RPM,2=Full RPM)
   * @param initialTemperature initial temperature
   * @param finalTemperature final temperature
   * @param annealRate anneal rate for RPM
   * @param useCentroids if true center points first
   * @param initialTransformaiton use this to initialize the mappings
   * @param debug if true print extra messages
   * @returns number of iterations if  success (>=1), or 0 if failed
   */
  virtual int run(int transformMode=2,
                  int correspondenceMode=2,
                  float initialTemperature=5.0,
                  float finalTemperature=1.0,
                  float annealRate=0.93,
                  int useCentroids=1,
                  bisMatrixTransformation* initialTransformation=NULL,
                  int debug=0);
  
protected:
  
  /** Description */
  bisMatrixTransformation* Output;

private:

  /** Copy constructor disabled to maintain shared/unique ptr safety */
  bisLinearRPMRegistration(const bisLinearRPMRegistration&);

  /** Assignment disabled to maintain shared/unique ptr safety */
  void operator=(const bisLinearRPMRegistration&);  


};

#endif
