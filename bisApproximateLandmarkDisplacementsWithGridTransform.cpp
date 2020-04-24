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

#include "bisApproximateLandmarkDisplacementsWithGridTransform.h"
#include "bisImageAlgorithms.h"
#include "bisUtil.h"
#include <sstream>
#include <iomanip>
#include <time.h>

static int count=0;

bisApproximateLandmarkDisplacementsWithGridTransform::bisApproximateLandmarkDisplacementsWithGridTransform(std::string s) : bisOptimizableAlgorithm(s)
{
  this->class_name="bisApproximateLandmarkDisplacementsWithGridTransform";
  this->enable_feedback=0;
}

bisApproximateLandmarkDisplacementsWithGridTransform::~bisApproximateLandmarkDisplacementsWithGridTransform()
{
  this->lastSmoothness=-1.0;
  this->lastSimilarity=-1.0;
  this->debug_flag=0;
}


void bisApproximateLandmarkDisplacementsWithGridTransform::generateFeedback(std::string input)
{

  std::cout << input << "  (" << std::fixed << std::setw(5) << this->lastSimilarity << "," << std::setw(5) << this->lastSmoothness << ")" << std::endl;
}


void bisApproximateLandmarkDisplacementsWithGridTransform::generateFeedback2(std::string input)
{
  std::cout << input << std::endl;
}

// Optimizer Stuff
float bisApproximateLandmarkDisplacementsWithGridTransform::computeValue(std::vector<float>& position)
{
  this->currentGridTransformation->setParameterVector(position);

  int numpoints=this->sourcePoints->getNumRows();
  float* inp_pts=this->sourcePoints->getData();
  float* out_pts=this->targetPoints->getData();
  float* wgt_pts=this->sourceWeights->getData();


  float sum=0.0;
  for (int i=0;i<numpoints;i++)
    {
      float x[3],tx[3],x2[3];
      for (int ia=0;ia<=2;ia++) {
        x[ia]=inp_pts[i*3+ia];
        tx[ia]=out_pts[i*3+ia];
      }
      this->currentGridTransformation->transformPoint(x,x2);
      if (count < 1 && i %1400 == 0 ) 
        std::cout << "x=" << x[0] << "," << x[1] << "," << x[2] << "-->" << x2[0] << "," << x2[1] << "," << x[2] << std::endl;
      

      double d=pow(tx[0]-x2[0],2.0f)+pow(tx[1]-x2[1],2.0f)+pow(tx[2]-x2[2],2.0f);
      float w=wgt_pts[i];

      if (count <1 && i %1400 == 0 ) 
        std::cout << "tx=" << tx[0] << "," << tx[1] << "," << tx[2] << " and d=" << d << " w=" << w << std::endl;

      
      sum+=w*sqrt(d);
    }

  this->lastSimilarity=sum/float(numpoints);
  float v=this->lastSimilarity;
  if (this->lambda>0.0)
    {
      this->lastSmoothness=this->currentGridTransformation->getTotalBendingEnergy();
      v+=this->lambda*this->lastSmoothness;
    }
  if (count < 11) {
    std::cout << "sum=" << sum << " last=" << this->lastSimilarity << " v=" << v << std::endl;
  }
  count=count+1;
  return v;
}

float bisApproximateLandmarkDisplacementsWithGridTransform::computeValueFunctionPiece(int cp)
{

  int numpoints=this->sourcePoints->getNumRows();
  int npc=this->gridPointList[cp].size();
  int numc=this->currentGridTransformation->getNumberOfControlPoints();
  float* inp_pts=this->sourcePoints->getData();
  float* out_pts=this->targetPoints->getData();
  float* wgt_pts=this->sourceWeights->getData();
  float sum=0.0;

  for (int l=0;l<npc;l++)
    {
      int pt=this->gridPointList[cp][l];
      float wb=this->gridPointWeight[cp][l];

      float x[3],tx[3],x2[3];

      for (int ia=0;ia<=2;ia++) {
        x[ia]=inp_pts[pt*3+ia];
        tx[ia]=out_pts[pt*3+ia];
      }
      
      this->currentGridTransformation->transformPoint(x,x2);
      double d=pow(tx[0]-x2[0],2.0f)+pow(tx[1]-x2[1],2.0f)+pow(tx[2]-x2[2],2.0f);
      float w=wgt_pts[pt];
      sum+=wb*w*sqrt(d);
    }

  float v=sum/float(numpoints);
  if (this->lambda>0.0)
    {
      float scale=0.01f*(numc);
      float sm=this->currentGridTransformation->getBendingEnergyAtControlPoint(cp,scale);
      v+=this->lambda*sm;
    }
  
  return v;
}

float bisApproximateLandmarkDisplacementsWithGridTransform::computeGradient(std::vector<float>& params,std::vector<float>& grad)
{
  unsigned int numc=this->currentGridTransformation->getNumberOfControlPoints();
  if (params.size()!=numc*3 || grad.size()!=params.size()) {
    std::cerr << "Bad dimensions for computing grdient optimization in grid transform";
    return 0;
  }

  this->currentGridTransformation->setParameterVector(params);
  float* dispfield=this->currentGridTransformation->getData();
  float GradientNorm = 0.000001f;

  for (unsigned int cp_index=0;cp_index<numc;cp_index++) { 
    for (int coord=0;coord<=2;coord++)
      {
        int index=cp_index+coord*numc;
        dispfield[index]=params[index]-this->stepsize;
        float a=this->computeValueFunctionPiece(cp_index);
        dispfield[index]=params[index]+this->stepsize;
        float b=this->computeValueFunctionPiece(cp_index);
        dispfield[index]=params[index];
        float g=-0.5f*(b-a)/this->stepsize;
        grad[index]=g;
        GradientNorm+=g*g;
      }
  }
  
  GradientNorm = float( sqrt(GradientNorm));
  for (unsigned int i=0;i<grad.size(); i++)
    grad[i]=grad[i]/GradientNorm;
  return GradientNorm;
}

  
int bisApproximateLandmarkDisplacementsWithGridTransform::checkInputParameters(bisJSONParameterList* plist)
{
  std::unique_ptr<bisJSONParameterList> tmp(new bisJSONParameterList(this->name+":plist"));
  this->internalParameters=std::move(tmp);

  
  this->internalParameters->setFloatValue("lambda",bisUtil::frange(plist->getFloatValue("lambda",0.0f),0.0f,1.0f));
  this->internalParameters->setFloatValue("stepsize",bisUtil::frange(plist->getFloatValue("stepsize",1.0f),0.05f,4.0f));
  this->internalParameters->setIntValue("steps",bisUtil::irange(plist->getIntValue("steps",2),1,100));
  this->internalParameters->setIntValue("iterations",bisUtil::irange(plist->getIntValue("iterations",15),1,100));
  this->internalParameters->setFloatValue("tolerance",bisUtil::frange(plist->getFloatValue("tolerance",0.001f),0.0f,0.5f));
  this->lambda=this->internalParameters->getFloatValue("lambda",0.0f);
  this->stepsize=  this->internalParameters->getFloatValue("stepsize",1.0f);
  if (this->enable_feedback)
    this->internalParameters->print("Approximate Landmark Displacement Field");
  return 1;
}

double B(int i, double t)
{
  switch (i) 
    {
    case 0:
      return (1-t)*(1-t)*(1-t)/6.0;
      break;
      
    case 1:
      return (3*t*t*t - 6*t*t + 4)/6.0;
    case 2:
      return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
      
    case 3:
      return (t*t*t)/6.0;
    }
  return 0;
}


void bisApproximateLandmarkDisplacementsWithGridTransform::initializePointLists() {

  const float thr=0.01;
 
  int numc=this->currentGridTransformation->getNumberOfControlPoints();
  float ori[3],spa[3];
  int dim[3];
  this->currentGridTransformation->getGridOrigin(ori);
  this->currentGridTransformation->getGridSpacing(spa);
  this->currentGridTransformation->getGridDimensions(dim);

  this->gridPointList.clear();
  this->gridPointWeight.clear();
  for (int i=0;i<numc;i++) {
    std::vector<int> a;
    this->gridPointList.push_back(a);
    std::vector<float> b;
    this->gridPointWeight.push_back(b);

  }

  int gridslicedims=dim[0]*dim[1];
  int numpoints=this->sourcePoints->getNumRows();
  float* pts=this->sourcePoints->getData();
  for (int node=0;node<numpoints;node++)
    {
      double p1[3],s[3];
      int  lat[3];
      for (int ia=0;ia<=2;ia++) {
        p1[ia]=pts[node*3+ia];
        float x=(p1[ia]-ori[ia])/spa[ia];
        lat[ia]=int(floor(x+0.0001));
        s[ia]=x-lat[ia];
      }

      if (node==10 && this->debug_flag>1) {
        std::cout << "__ Point 100 " << p1[0] << "," << p1[1] << "," << p1[2] << std::endl;
        std::cout << "__   lat = " << lat[0] << "," << lat[1] << "," << lat[2] << std::endl;
        std::cout << "__     s = " << s[0] << "," << s[1] << "," << s[2] << std::endl;
      }
      
      for (int k = 0; k < 4; k++)
        {
          int K = bisUtil::irange(k + lat[2] - 1,0,dim[2]-1);
          for (int j = 0; j < 4; j++)
            {
              int J = bisUtil::irange(j + lat[1] - 1,0,dim[1]-1);
              for (int i = 0; i < 4; i++)
                {
                  int I = bisUtil::irange(i + lat[0] - 1,0,dim[0]-1);
                  int   cpoint=I+J*dim[0]+K*gridslicedims;
                  double wgt= B(i, s[0]) * B(j, s[1]) * B(k, s[2]);
                  
                  if (wgt>=thr)
                    {
                      unsigned int found=0 ,index=0;
                      
                      while (found ==0 && index <this->gridPointList[cpoint].size()) {
                        if (this->gridPointList[cpoint][index]==node) {
                          found=1;
                          this->gridPointWeight[cpoint][index]+=wgt;
                          if (node==10 && this->debug_flag>1) {
                            std::cout << "____ " << i << "," << j << "," << k << " ---> " << I << "," << J << "," << K << "  --> " << cpoint << " " << wgt << std::endl;
                            std::cout << "Incrementing " << cpoint << std::endl;
                          }
                        } else {
                          index=index+1;
                        }
                      }
                      
                      if (!found) {
                        if (node==10 && this->debug_flag>1) {
                          std::cout << "____ " << i << "," << j << "," << k << " ---> " << I << "," << J << "," << K << "  --> " << cpoint << " " << wgt << std::endl;
                          std::cout << "Adding to " << cpoint << std::endl;
                        }
                        
                        this->gridPointList[cpoint].push_back(node);
                        this->gridPointWeight[cpoint].push_back(wgt);
                      }
                    }
                }
            }
        }
    }
  
  /*int lst[3]={ 12,12,12 };
  for (int i=0;i<1;i++) {
    int node=lst[i];
    int sz=this->gridPointList[node].size();
    std::cout << "For cp = " << node << std::endl << "\t";
    for (int i=0;i<sz;i++) {
      std::cout << "(" << this->gridPointList[node][i] << "," << this->gridPointWeight[node][i] << ") ";
    }
    std::cout << std::endl;
    }*/
}


// Set Parameters and Run
float bisApproximateLandmarkDisplacementsWithGridTransform::run(bisSimpleMatrix<float>* in_sourceLandmarks,
                                                                bisSimpleMatrix<float>* in_targetLandmarks,
                                                                bisSimpleMatrix<float>* in_sourceWeights,
                                                                bisGridTransformation* transformation,
                                                                bisJSONParameterList* plist,
                                                                int dbg)
{

  this->debug_flag=dbg;
  if (this->debug_flag<2)
    count=10;
  
  this->currentGridTransformation=transformation;
  this->sourcePoints=in_sourceLandmarks;
  this->targetPoints=in_targetLandmarks;
  this->sourceWeights=in_sourceWeights;

  int rows1=this->sourcePoints->getNumRows();
  int rows2=this->targetPoints->getNumRows();
  int rows3=this->sourceWeights->getNumRows();

  int cols1=this->sourcePoints->getNumCols();
  int cols2=this->targetPoints->getNumCols();
  int cols3=this->sourceWeights->getNumCols();

  if (rows1 < 4 || rows1 != rows2 || rows1 !=rows3 || cols1 !=3  || cols2 !=3 || cols3!=1 ) {
    std::cerr << "Bad Landmark sets =" << rows1 << "*" << cols1 << " and " << rows2 << "*" << cols2 << " and " << rows3 << "*" << cols3 << std::endl;
    return -1.0;
  }
  this->initializePointLists();

  
  if (this->enable_feedback)
    this->generateFeedback2("++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ +");
  this->checkInputParameters(plist);

  std::stringstream strss;
  strss.precision(5);

  int numsteps=   this->internalParameters->getIntValue("steps");
  int iterations=this->internalParameters->getIntValue("iterations");
  float tolerance=this->internalParameters->getFloatValue("tolerance",0.001f);

  strss.clear();
  std::stringstream strss2;
  strss2 << "++   Beginning to appproximate landmark displacement field . numsteps= " << numsteps << "  tolerance=" << tolerance << " lambda=" << this->lambda;
  this->generateFeedback2(strss2.str());
  
  int numdof=this->currentGridTransformation->getNumberOfDOF();
  this->generateFeedback2("++  ");
  std::stringstream strss3;
  strss3 << "++   Approx numdof=" << numdof << " step=" << this->stepsize;
  this->generateFeedback2(strss3.str());
  this->generateFeedback2("++  ");


  std::unique_ptr<bisOptimizer> optimizer(new bisOptimizer(this));
  std::vector<float> position(numdof);
  // Get current state ...
  this->currentGridTransformation->getParameterVector(position);
  float last=0.0;
  for (int step=numsteps;step>=1;step=step-1)
    {
      std::cout << "~~~~ optimizing. step = " << step << ", iterations = " << iterations << " cur=" << this->stepsize;
      strss.clear();
      this->generateFeedback2(strss.str());
      last=optimizer->computeConjugateGradient(position,iterations,tolerance);
      this->stepsize*=0.5;
    }
  this->generateFeedback2("++  ");
  return last;
}




