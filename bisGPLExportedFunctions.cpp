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

#include "bisSimpleDataStructures.h"
#include "bisGPLExportedFunctions.h"
#include "bisJSONParameterList.h"
#include "bisMatrixTransformation.h"
#include "bisGridTransformation.h"
#include "bisComboTransformation.h"
#include "bisDataObjectFactory.h"
#include "bisDTIAlgorithms.h"
#include "bisImageSegmentationAlgorithms.h"
#include "bisLinearImageRegistration.h"
#include "bisNonLinearImageRegistration.h"
#include "bisApproximateDisplacementField.h"
#include "bisLinearRPMRegistration.h"
#include <memory>


int uses_gpl() {
  return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------
// Linear Image Registration
// --------------------------------------------------------------------------------------------------------------------------------------------------------

unsigned char* runLinearRegistrationWASM(unsigned char* reference,
					 unsigned char* target,
					 unsigned char* initial_ptr,
					 const char* jsonstring,
					 int debug)
{  
  if (debug)
    std::cout << "_____ Beginning runLinearRegistrationJSON" << std::endl;
  
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  int return_vector=params->getBooleanValue("return_vector",0);
  
  //  if(debug)
  //    params->print("from runLinearRegistrationJSON","_____");

  std::shared_ptr<bisSimpleImage<float> > reference_image(new bisSimpleImage<float>("reference_image_json"));
  if (!reference_image->linkIntoPointer(reference))
    return 0;
  
  std::shared_ptr<bisSimpleImage<float> > target_image(new bisSimpleImage<float>("target_image_json"));
  if (!target_image->linkIntoPointer(target))
    return 0;

  std::unique_ptr<bisMatrixTransformation> initial_transformation(new bisMatrixTransformation("parse_initial"));
  initial_transformation->identity();
  
  std::unique_ptr<bisSimpleMatrix<float> > initial_matrix(new bisSimpleMatrix<float>("initial_matrix_json"));
  if (initial_ptr!=0)
    {
      if (!initial_matrix->linkIntoPointer(initial_ptr))
	return 0;
      if (!initial_transformation->setSimpleMatrix(initial_matrix.get()))
	return 0;
    }

  std::unique_ptr<bisLinearImageRegistration> reg(new bisLinearImageRegistration("linear registration"));
  reg->setReferenceImage(reference_image);
  reg->setTargetImage(target_image);
  reg->setInitialTransformation(initial_transformation.get());
  reg->run(params.get());


  
  if (return_vector==0)
    {
      std::unique_ptr<bisSimpleMatrix<float> > output(reg->getOutputMatrix());
      bisUtil::mat44 m; output->exportMatrix(m);
      return output->releaseAndReturnRawArray();
    }

  std::unique_ptr<bisSimpleVector<float> > output(reg->getTransformationParameterVector());
  int length=output->getLength();

  Eigen::MatrixXf outmat=Eigen::MatrixXf::Zero(length,1);
  for (int i=0;i<length;i++)
    outmat(i,0)=output->getData()[i];
  return bisEigenUtil::serializeAndReturn(outmat,"param_vector");
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------
// Non Linear Image Registration
// --------------------------------------------------------------------------------------------------------------------------------------------------------

unsigned char* runNonLinearRegistrationWASM(unsigned char* reference,
					    unsigned char* target,
					    unsigned char* initial_ptr,
					    const char* jsonstring,
					    int debug)
{

  if (debug)
    std::cout << "_____ Beginning runNonLinearRegistrationJSON" << std::endl;
  
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  //  if(debug)
  //    params->print("from runNonLinearRegistrationJSON","_____");

  std::shared_ptr<bisSimpleImage<float> > reference_image(new bisSimpleImage<float>("reference_image_json"));
  if (!reference_image->linkIntoPointer(reference))
    return 0;
  
  std::shared_ptr<bisSimpleImage<float> > target_image(new bisSimpleImage<float>("target_image_json"));
  if (!target_image->linkIntoPointer(target))
    return 0;

  std::unique_ptr<bisNonLinearImageRegistration> reg(new bisNonLinearImageRegistration("nonlinear"));
  reg->setReferenceImage(reference_image);
  reg->setTargetImage(target_image);

  if (initial_ptr) {
    std::shared_ptr<bisAbstractTransformation> initial_transformation=bisDataObjectFactory::deserializeTransformation(initial_ptr,"initialxform");
    reg->setInitialTransformation(initial_transformation);
  }
  reg->run(params.get());

  std::shared_ptr<bisComboTransformation> output(reg->getOutputTransformation());
  unsigned char* pointer=output->serialize();
  
 
  return pointer;
}


// -----------------------------------------------------------------------------------------------------
// Image Segmentation
// -----------------------------------------------------------------------------------------------------
template <class BIS_TT> unsigned char* segmentImageTemplate(unsigned char* input,bisJSONParameterList* params,int debug,BIS_TT*)
{
  std::unique_ptr<bisSimpleImage<BIS_TT> > inp_image(new bisSimpleImage<BIS_TT>("inp_image"));
  if (!inp_image->linkIntoPointer(input))
    return 0;

  int frame=params->getIntValue("frame",0);
  int component=params->getIntValue("component",0);
  int numclasses=params->getIntValue("numclasses",3);
  float maxsigmaratio=params->getFloatValue("maxsigmaratio",0.2f);
  int robust=params->getBooleanValue("robust",1);
  int smhisto=params->getBooleanValue("smoothhisto",1);
  int num_bins=params->getIntValue("numbins",256);
  int max_iterations=params->getIntValue("maxiterations",30);
  float convergence=params->getFloatValue("convergence",0.05f);
  int use_variance=params->getBooleanValue("usevariance",1);
  float smoothness=params->getFloatValue("smoothness",0.0f);
  int mrfiterations=params->getIntValue("mrfiterations",8);
  int internaliterations=params->getIntValue("internaliterations",4);
  float noisesigma2=params->getFloatValue("noisesigma2",25.0f);
  float mrfconvergence=params->getFloatValue("mrfconvergence",0.2f);
  
  
  if (debug) {
    std::cout << "Image Segmentation Parameters: smoothness=" << smoothness << std::endl;
    std::cout << "Parsed parameters  frame=" << frame << " comp=" << component << std::endl;
    std::cout << "\t numclasses=" << numclasses << " maxsigmaratio=" << maxsigmaratio << " maxiter=" << max_iterations << " conv=" << convergence << " numbins=" << num_bins << std::endl;
    std::cout << "\t robust=" << robust << " smoothisto=" << smhisto << std::endl;
    if (smoothness>0.0)
      std::cout << "\t MRF -- mrfiterations=" << mrfiterations << "internaliter=" << internaliterations << " noisesigma2=" << noisesigma2 << " mrfconvergence=" << mrfconvergence << std::endl;
    std::cout << "-----------------------------------" << std::endl;
  }
  
  
  if (debug)
    std::cout << std::endl << "..... Begin Histogram Segmentation" << std::endl;
  std::unique_ptr<bisSimpleImage<short> > out_image(bisImageSegmentationAlgorithms::histogramSegmentation(inp_image.get(),
													  numclasses,maxsigmaratio,
													  max_iterations,convergence,use_variance,
													  num_bins,robust,smhisto,
													  frame,component));
  if (debug)
    std::cout << std::endl << "..... Histogram Segmentation done" << std::endl;

  if (smoothness>0.0)
    {
      if (debug)
	std::cout << std::endl << "..... Begin MRF Segmentation" << std::endl;
      
      bisImageSegmentationAlgorithms::doMRFSegmentation(inp_image.get(),
							out_image.get(),
							smoothness,
							noisesigma2,
							convergence,
							mrfiterations,internaliterations,
							frame,component);
      if (debug)
	std::cout << std::endl << "..... MRF Segmentation done" << std::endl;
    }
  
  return out_image->releaseAndReturnRawArray();
}

unsigned char* segmentImageWASM(unsigned char* input,
				     const char* jsonstring,int debug)
{
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  int ok=params->parseJSONString(jsonstring);
  if (!ok) 
    return 0;

  if(debug)
    params->print();

  int* header=(int*)input;
  int in_type=header[1];

  switch (in_type)
      {
	bisvtkTemplateMacro( return segmentImageTemplate(input,params.get(),debug, static_cast<BIS_TT*>(0)));
      }
  return 0;
}

// ------------------------------------------------------------------------------------
unsigned char* approximateDisplacementFieldWASM(unsigned char* dispfield_ptr,
						unsigned char* initial_grid_ptr,
						const char* jsonstring,
						int debug)
{
  if (debug)
    std::cout << "_____ Beginning approximateDisplacementFieldJSON" << std::endl;
  
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  //  if(debug)
  //    params->print("from runApproximateDisplacementField","_____");

  std::unique_ptr<bisSimpleImage<float> > disp_field(new bisSimpleImage<float>("disp_field_json"));
  if (!disp_field->linkIntoPointer(dispfield_ptr))
    return 0;

  
  std::unique_ptr<bisGridTransformation> initial_grid(new bisGridTransformation("initial_grid_json"));
  if (!initial_grid->deSerialize(initial_grid_ptr))
    return 0;


  std::unique_ptr<bisApproximateDisplacementField> reg(new bisApproximateDisplacementField("approx"));
  reg->run(disp_field.get(),initial_grid.get(),params.get());


  unsigned char* pointer=initial_grid->serialize();
  return pointer;
}

// ------------------------------------------------------------------------------------
unsigned char* approximateDisplacementFieldWASM2(unsigned char* dispfield_ptr,
						 const char* jsonstring,
						 int debug)
{
  if (debug)
    std::cout << "_____ Beginning approximateDisplacementFieldJSON" << std::endl;
  
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  if(debug)
    params->print("from runApproximateDisplacementField","_____");

  std::unique_ptr<bisSimpleImage<float> > disp_field(new bisSimpleImage<float>("disp_field_json"));
  if (!disp_field->linkIntoPointer(dispfield_ptr))
    return 0;


  float spacing=params->getFloatValue("spacing",10.0);
  int dim[3]; disp_field->getImageDimensions(dim);
  float spa[3]; disp_field->getImageSpacing(spa);
  int griddim[3];
  float gridspa[3],gridori[3];
  
  for (int ia=0;ia<=2;ia++) {
    griddim[ia] = int((dim[ia]*spa[ia])/spacing)+1;
    gridspa[ia] = spacing;
    gridori[ia] = -0.5*(griddim[ia]*gridspa[ia]-dim[ia]*spa[ia]);
  }
  
  if (debug)  {
    std::cout << "\t input spacing of grid=" << spacing << std::endl;
    std::cout << "\t initialized grid: dim=" << griddim[0] << "," << griddim[1] << "," << griddim[2] << ",";
    std::cout << "\t spa=" << gridspa[0] << "," << gridspa[1] << "," << gridspa[2] << ",";
    std::cout << "\t ori=" << gridori[0] << "," << gridori[1] << "," << gridori[2] << std::endl;
  }

  std::unique_ptr<bisGridTransformation> output_grid(new bisGridTransformation("output_grid"));
  output_grid->initializeGrid(griddim,gridspa,gridori,1);

  std::unique_ptr<bisApproximateDisplacementField> reg(new bisApproximateDisplacementField("approx"));
  reg->run(disp_field.get(),output_grid.get(),params.get());

  unsigned char* pointer=output_grid->serialize();
  return pointer;
}
// -----------------------------------------------------------------------------------------------------
// Regularize Objectmap
// -----------------------------------------------------------------------------------------------------
  unsigned char* regularizeObjectmapWASM(unsigned char* input,const char* jsonstring,int debug)
{
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  int ok=params->parseJSONString(jsonstring);
  if (!ok) 
    return 0;

  if(debug)
    params->print();

  std::unique_ptr<bisSimpleImage<short> > in_image(new bisSimpleImage<short>("in_image"));
  // This 1 at the end means "copy data" as this is an in/out function
  if (!in_image->linkIntoPointer(input,1))
    return 0;

  float smoothness=params->getFloatValue("smoothness",2.0f);
  int maxiter=params->getIntValue("iterations",8);
  int internal_iter=params->getIntValue("internaliterations",4);
  float mrfconvergence=params->getFloatValue("convergence",0.2f);
  
  
  if (debug) {
    std::cout << "Objectmap Regularization Parameters: smoothness=" << smoothness << std::endl;
    std::cout << "\t iterations=" << maxiter << "internaliter=" << internal_iter <<  " convergence=" << mrfconvergence << std::endl;
    std::cout << "-----------------------------------" << std::endl;
  }
  
  
  if (debug)
    std::cout << std::endl << "..... Begin Objectmap Regularization" << std::endl;
  std::unique_ptr<bisSimpleImage<short> > out_image(bisImageSegmentationAlgorithms::doObjectMapRegularization(in_image.get(),
                                                                                                              smoothness,
                                                                                                              mrfconvergence,
                                                                                                              maxiter,internal_iter));
  if (debug)
    std::cout << std::endl << "..... Objectmap Regularization done " << std::endl;

  
  return out_image->releaseAndReturnRawArray();
}

// ---------------------------------------- DTI Stuff --------------------------------
unsigned char* computeDTITensorFitWASM(unsigned char* input_ptr,
				       unsigned char* baseline_ptr,
				       unsigned char* mask_ptr,
				       unsigned char* directions_ptr,
				       const char* jsonstring,
				       int debug)
{

  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  if(debug)
    params->print("computeCorrelationMatrixJSON","_____");

  
  std::unique_ptr<bisSimpleImage<short> > in_image(new bisSimpleImage<short>("input_dti_data"));
  if (!in_image->linkIntoPointer(input_ptr))
    return 0;

  std::unique_ptr<bisSimpleImage<short> > baseline_image(new bisSimpleImage<short>("baseline_dti_data"));
  if (!baseline_image->linkIntoPointer(baseline_ptr))
    return 0;
  
  

  Eigen::MatrixXf directions;
  std::unique_ptr<bisSimpleMatrix<float> > s_matrix(new bisSimpleMatrix<float>("directions"));
  if (!bisEigenUtil::deserializeAndMapToEigenMatrix(s_matrix.get(),directions_ptr,directions,debug))
    return 0;

  std::unique_ptr<bisSimpleImage<float> > out_image(new bisSimpleImage<float>("output_dti_data"));


  float bvalue=params->getFloatValue("bvalue",1000.0f);
  if (debug)
    std::cout << "Beginning Fit " << bvalue << std::endl;

  int ok=0;
  if (mask_ptr==0)
    {
      if (debug)
	std::cout << "Not using mask " << std::endl;
      ok=bisDTIAlgorithms::computeTensorFit(in_image.get(),
					    baseline_image.get(),
					    0,
					    directions,
					    bvalue,
					    out_image.get());
    }
  else
    {
      std::unique_ptr<bisSimpleImage<unsigned char> > mask_image(new bisSimpleImage<unsigned char>("mask_dti_data"));
      if (!mask_image->linkIntoPointer(mask_ptr))
	return 0;
      if (debug)
	std::cout << "Using mask " << std::endl;
      ok=bisDTIAlgorithms::computeTensorFit(in_image.get(),
					    baseline_image.get(),
					    mask_image.get(),
					    directions,
					    bvalue,
					    out_image.get());
    }


  if (debug)
    std::cout << "Fitting Done " << ok << std::endl;

  return out_image->releaseAndReturnRawArray();
}


/** Computes Eigenvalues and Eigenvector as a single image of 4 components x 3 frames
 * component 0 = eigenvalues
 * components 1-3 eigenvectors
 * frames are x,y,z
 * @param tensor the input dti tensor (from computeTensorFit)
 * @param mask the input mask image (can be NULL,0)
 * @param eigenSystem the output images as defined above
 * @returns 1 if success, 0 if failed */
unsigned char* computeTensorEigenSystemWASM(unsigned char* input_ptr,
					    unsigned char* mask_ptr,
					    int debug) 
{
  std::unique_ptr<bisSimpleImage<float> > in_image(new bisSimpleImage<float>("input_eigensystem_data"));
  if (!in_image->linkIntoPointer(input_ptr))
    return 0;


  std::unique_ptr<bisSimpleImage<unsigned char> > mask_image(new bisSimpleImage<unsigned char>("mask_dti_data"));
  bisSimpleImage<unsigned char>* mask=0;
  if (mask_ptr!=0)
    {
      if (!mask_image->linkIntoPointer(mask_ptr))
	return 0;
      mask=mask_image.get();
    }
  
  if (debug)
    {
      std::cout << "Beginning Compute Tensor Eigen System ";
      if (mask!=0)
	std::cout <<  "using mask";
      std::cout << std::endl;
    }
  
  std::unique_ptr<bisSimpleImage<float> > output(new bisSimpleImage<float>("out_dti_eigensystem"));
  
  int ok=bisDTIAlgorithms::computeTensorEigenSystem(in_image.get(),mask,output.get());

  if (debug)
    std::cout << "Done Computing ok=" << ok << std::endl;
  
  return output->releaseAndReturnRawArray();
  
}
/** Compute DTI Tensor Invariants
 * @param input_ptr the image tensor eigensystem as a serialized array
 * @param mask_ptr the Mask Image (optional, set this to 0) as a serialized array
 * @param jsonstring { "mode": 0 } // mode 0=FA, 1=RA etc. -- see bisDTIAlgorithms::computeTensorInvariants
 * @param debug if > 0 print debug messages
 * @returns a pointer to the invarient image */
BISEXPORT unsigned char* computeDTITensorInvariantsWASM(unsigned char* input_ptr,
							unsigned char* mask_ptr,
							const char* jsonstring,
							int debug)
{
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  int ok=params->parseJSONString(jsonstring);
  if (!ok) 
    return 0;

  if(debug)
    params->print();

  std::unique_ptr<bisSimpleImage<float> > in_image(new bisSimpleImage<float>("input_eigensystem_data"));
  if (!in_image->linkIntoPointer(input_ptr))
    return 0;


  std::unique_ptr<bisSimpleImage<unsigned char> > mask_image(new bisSimpleImage<unsigned char>("mask_dti_data"));
  if (mask_ptr!=0)
    {
      if (!mask_image->linkIntoPointer(mask_ptr))
	return 0;
    }
  else
    {
      std::unique_ptr<bisSimpleImage<unsigned char> >m(bisImageAlgorithms::createMaskImage<float>(in_image.get()));
      mask_image=std::move(m);
    }

  int mode=params->getIntValue("mode",0);
  
  if (debug)
    std::cout << "Beginning Compute Tensor Invariants mode=" << mode  << std::endl;
  
  std::unique_ptr<bisSimpleImage<float> > output(new bisSimpleImage<float>("out_dti_eigensystem"));
  
  ok=bisDTIAlgorithms::computeTensorInvariants(in_image.get(),mask_image.get(),mode,output.get());

  if (debug)
    std::cout << "Done Computing ok=" << ok << std::endl;
  
  return output->releaseAndReturnRawArray();


}

/** Compute DTI Orientation Map
 * @param input_ptr the image tensor eigensystem as a serialized array
 * @param mask_ptr the Mask Image (optional, set this to 0) as a serialized array
 * @param magnitude_ptr the Magnitude Image (e.g. FA map) (optional, set this to 0) as a serialized array
 * @param jsonstring { "scaling": 1.0 } Optional extra scaling
 * @param debug if > 0 print debug messages
 * @returns a pointer to the colormap image */
BISEXPORT unsigned char* computeDTIColorMapImageWASM(unsigned char* input_ptr,
						     unsigned char* mask_ptr,
						     unsigned char* magnitude_ptr,
						     const char* jsonstring,
						     int debug)
{
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  int ok=params->parseJSONString(jsonstring);
  if (!ok) 
    return 0;

  if(debug)
    params->print();

  std::unique_ptr<bisSimpleImage<float> > in_image(new bisSimpleImage<float>("input_eigensystem_data"));
  if (!in_image->linkIntoPointer(input_ptr))
    return 0;


  std::unique_ptr<bisSimpleImage<unsigned char> > mask_image(new bisSimpleImage<unsigned char>("mask_dti_data"));
  if (mask_ptr!=0)
    {
      if (!mask_image->linkIntoPointer(mask_ptr))
	return 0;
    }
  else
    {
      std::unique_ptr<bisSimpleImage<unsigned char> >m(bisImageAlgorithms::createMaskImage<float>(in_image.get()));
      mask_image=std::move(m);
    }
  

  std::unique_ptr<bisSimpleImage<float> > magn_image(new bisSimpleImage<float>("magnitude_dti_data"));
  bisSimpleImage<float>* magn=0;
  if (magnitude_ptr!=0)
    {
      if (!magn_image->linkIntoPointer(magnitude_ptr))
	return 0;
      magn=magn_image.get();
    }


  float scaling=params->getFloatValue("scaling",1.0);
  
  if (debug)
    std::cout << "Beginning Compute Tensor Colormap scaling=" << scaling  << std::endl;
  
  std::unique_ptr<bisSimpleImage<unsigned char> > output(new bisSimpleImage<unsigned char>("out_dti_colormap"));
  
  ok=bisDTIAlgorithms::computeTensorColormap(in_image.get(),mask_image.get(),magn,scaling,output.get());
  
  if (debug)
    std::cout << "Done Computing ok=" << ok << std::endl;
  
  return output->releaseAndReturnRawArray();
}




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
unsigned char*  runWeightedLinearRegistrationWASM(unsigned char* reference_ptr,
                                                  unsigned char* target_ptr,
                                                  unsigned char* reference_weight_ptr,
                                                  unsigned char* target_weight_ptr,
                                                  unsigned char* initial_ptr,
                                                  const char* jsonstring,
                                                  int debug)
{
  if (debug)
    std::cout << "_____ Beginning runWeightedLinearRegistrationJSON" << std::endl;
  
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  int return_vector=params->getBooleanValue("return_vector",0);
  
  //  if(debug)
  //    params->print("from runWeightedLinearRegistrationJSON","_____");

  std::shared_ptr<bisSimpleImage<float> > reference_image(new bisSimpleImage<float>("reference_image_json"));
  if (!reference_image->linkIntoPointer(reference_ptr))
    return 0;
  
  std::shared_ptr<bisSimpleImage<float> > target_image(new bisSimpleImage<float>("target_image_json"));
  if (!target_image->linkIntoPointer(target_ptr))
    return 0;

  std::shared_ptr<bisSimpleImage<short> > reference_weight_image(new bisSimpleImage<short>("reference_weight_json"));
  if (!reference_weight_image->linkIntoPointer(reference_weight_ptr))
    return 0;

  
  std::unique_ptr<bisMatrixTransformation> initial_transformation(new bisMatrixTransformation("parse_initial"));
  initial_transformation->identity();
  
  std::unique_ptr<bisSimpleMatrix<float> > initial_matrix(new bisSimpleMatrix<float>("initial_matrix_json"));
  if (initial_ptr!=0)
    {
      if (!initial_matrix->linkIntoPointer(initial_ptr))
	return 0;
      if (!initial_transformation->setSimpleMatrix(initial_matrix.get()))
	return 0;
    }

  std::unique_ptr<bisLinearImageRegistration> reg(new bisLinearImageRegistration("linear registration"));
  reg->setReferenceImage(reference_image);
  reg->setTargetImage(target_image);
  reg->setReferenceWeightImage(reference_weight_image);
  if (target_weight_ptr!=0) {
    std::shared_ptr<bisSimpleImage<short> > target_weight_image(new bisSimpleImage<short>("target_weight_json"));
    if (!target_weight_image->linkIntoPointer(reference_weight_ptr))
      return 0;
    reg->setTargetWeightImage(target_weight_image);
  }
  reg->setInitialTransformation(initial_transformation.get());
  reg->run(params.get());
  
  if (return_vector==0)
    {
      std::unique_ptr<bisSimpleMatrix<float> > output(reg->getOutputMatrix());
      bisUtil::mat44 m; output->exportMatrix(m);
      return output->releaseAndReturnRawArray();
    }

  std::unique_ptr<bisSimpleVector<float> > output(reg->getTransformationParameterVector());
  int length=output->getLength();

  Eigen::MatrixXf outmat=Eigen::MatrixXf::Zero(length,1);
  for (int i=0;i<length;i++)
    outmat(i,0)=output->getData()[i];
  return bisEigenUtil::serializeAndReturn(outmat,"param_vector");


}
  
  
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
unsigned char* runWeightedNonLinearRegistrationWASM(unsigned char* reference,
                                                    unsigned char* target,
                                                    unsigned char* reference_weight_ptr,
                                                    unsigned char* target_weight_ptr,
                                                    unsigned char* initial_ptr,
                                                    const char* jsonstring,
                                                    int debug)
{

  if (debug)
    std::cout << "_____ Beginning runWeightedNonLinearRegistrationJSON" << std::endl;
  
  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  //  if(debug)
  //    params->print("from runWeightedNonLinearRegistrationJSON","_____");

  std::shared_ptr<bisSimpleImage<float> > reference_image(new bisSimpleImage<float>("reference_image_json"));
  if (!reference_image->linkIntoPointer(reference))
    return 0;
  
  std::shared_ptr<bisSimpleImage<float> > target_image(new bisSimpleImage<float>("target_image_json"));
  if (!target_image->linkIntoPointer(target))
    return 0;

  std::shared_ptr<bisSimpleImage<short> > reference_weight_image(new bisSimpleImage<short>("reference_weight_json"));
  if (!reference_weight_image->linkIntoPointer(reference_weight_ptr))
    return 0;



  std::unique_ptr<bisNonLinearImageRegistration> reg(new bisNonLinearImageRegistration("nonlinear"));
  reg->setReferenceImage(reference_image);
  reg->setTargetImage(target_image);
  reg->setReferenceWeightImage(reference_weight_image);
  if (target_weight_ptr!=0) {
    std::shared_ptr<bisSimpleImage<short> > target_weight_image(new bisSimpleImage<short>("target_weight_json"));
    if (!target_weight_image->linkIntoPointer(reference_weight_ptr))
      return 0;
    reg->setTargetWeightImage(target_weight_image);
  }

  if (initial_ptr) {
    std::shared_ptr<bisAbstractTransformation> initial_transformation=bisDataObjectFactory::deserializeTransformation(initial_ptr,"initialxform");
    reg->setInitialTransformation(initial_transformation);
  }
  reg->run(params.get());

  std::shared_ptr<bisComboTransformation> output(reg->getOutputTransformation());
  unsigned char* pointer=output->serialize();
  
 
  return pointer;
}

/** run Linear RPM Registration using \link bisLinearImageRegistration  \endlink
 * @param reference serialized reference points unsigned char array 
 * @param target    serialized target points as unsigned char array 
 * @param initial_xform serialized initial transformation as unsigned char array 
 * @param reference_labels serialized reference labels unsigned char array 
 * @param target_labels serialized target labels unsigned char array 
 * @param jsonstring the parameter string for the algorithm  { numLandmarks: 1000, initialTemperature: 10.0, finalTemperature: 1.0,annealRate : 0.93, prefSampling : 1, 
 *                                                             trnasformMode :2 , correspondenceMode :2 , useCentroids : 1 }
 * @param debug if > 0 print debug messages
 * @returns a pointer to a serialized 4x4 matrix
 */
unsigned char*  runLinearRPMRegistrationWASM(unsigned char* in_reference,
                                             unsigned char* in_target,
                                             unsigned char* initial_ptr,
                                             unsigned char* in_reference_labels,
                                             unsigned char* in_target_labels,
                                             const char* jsonstring,
                                             int debug)
{

  std::unique_ptr<bisJSONParameterList> params(new bisJSONParameterList());
  if (!params->parseJSONString(jsonstring))
    return 0;

  int numlandmarks=params->getIntValue("numLandmarks",1000);
  int transformMode=params->getIntValue("transformMode",2);
  int correspondenceMode=params->getIntValue("correspondenceMode",2);
  int useCentroids=params->getIntValue("useCentroids",1);
  float initialTemperature=params->getFloatValue("initialTemperature",10.0);
  float finalTemperature=params->getFloatValue("finalTemperature",10.0);
  int iterPerTemp=params->getIntValue("iterPerTemp",5);
  float annealRate=params->getFloatValue("annealRate",10.0);
  int prefSampling=params->getIntValue("prefSampling",10.0);

  if (debug) {
    std::cout << " Linear RPM Parameters " << std::endl;
    std::cout << "          numlandmarks=" << numlandmarks << " transformMode=" << transformMode << " corrMode=" << correspondenceMode << " useCent=" << useCentroids << std::endl;
    std::cout << "          initialT=" << initialTemperature << " finalT=" << finalTemperature << " annealRate=" << annealRate << " iterPerT=" << iterPerTemp << " prefSampling=" << prefSampling << std::endl;
  }

  std::unique_ptr<bisSimpleMatrix<float> > ref_points(new bisSimpleMatrix<float>("ref_points_json"));
  if (!ref_points->linkIntoPointer(in_reference))
    return 0;
  if (debug)
    std::cout << "___ Ref Points = " << ref_points->getNumRows() << "*" << ref_points->getNumCols() << std::endl;

  std::unique_ptr<bisSimpleMatrix<float> > target_points(new bisSimpleMatrix<float>("target_points_json"));
  if (!target_points->linkIntoPointer(in_target))
    return 0;
  if (debug)
    std::cout << "___ Target Points = " << target_points->getNumRows() << "*" << target_points->getNumCols() << std::endl;

  std::unique_ptr<bisMatrixTransformation> initial_transformation(new bisMatrixTransformation("parse_initial"));
  initial_transformation->identity();
  std::unique_ptr<bisSimpleMatrix<float> > initial_matrix(new bisSimpleMatrix<float>("initial_matrix_json"));
  if (initial_ptr!=0)
    {
      if (!initial_matrix->linkIntoPointer(initial_ptr))
	return 0;
      if (!initial_transformation->setSimpleMatrix(initial_matrix.get()))
	return 0;
    }
  else
    {
      initial_transformation=0;
    }

  std::unique_ptr<bisLinearRPMRegistration> RPM(new bisLinearRPMRegistration());

  if (in_reference_labels!=0 && in_target_labels!=0) {
    // Use labels
    std::unique_ptr<bisSimpleVector<int> > ref_labels(new bisSimpleVector<int>("ref_labels_json"));
    std::unique_ptr<bisSimpleVector<int> > targ_labels(new bisSimpleVector<int>("targ_labels_json"));
    if (!ref_labels->linkIntoPointer(in_reference_labels) || 
        !targ_labels->linkIntoPointer(in_target_labels) ) {
      std::cerr << "__ bad labels" << std::endl;
      return 0;
    }
    RPM->initialize(ref_points.get(),target_points.get(),numlandmarks,prefSampling,ref_labels.get(),targ_labels.get(),debug);
  } else {
    RPM->initialize(ref_points.get(),target_points.get(),numlandmarks,1,0,0,debug);
  }

  RPM->run(transformMode,
           correspondenceMode,
           initialTemperature,
           finalTemperature,
           iterPerTemp,
           annealRate,
           useCentroids,
           initial_transformation.get(),
           debug);

  std::unique_ptr<bisSimpleMatrix<float> > output(RPM->getOutputMatrix());
  bisUtil::mat44 m; output->exportMatrix(m);
  return output->releaseAndReturnRawArray();
}

