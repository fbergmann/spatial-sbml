/**
 * \file    local.i
 * \brief   Java-specific SWIG support code for wrapping spatialSBML API
 * \author  Frank T. Bergmann
 *
 */

%include "arrays_java.i";
%include "typemaps.i";
%apply double* OUTPUT {double* array}
 
%extend SpatialSimulator
{

  void __getX(double* array)
  {
	int length; 
	array = self->getX(length);  
  }
  
  void __getZ(double* array)
  {
	int length; 
	array = self->getZ(length);
  }
  
  void __getY(double* array)
  {
	int length; 
	array = self->getY(length);
  }
  
  void __getVariable(const std::string &speciesId, double* array)
  {
	int length; 
	array = self->getVariable(speciesId, length);
  }
  void __getGeometry(const std::string &compartmentId, double* array)
  {
	int length; 
	array = self->getGeometry(compartmentId, length);
  }

}

%typemap(javacode) SpatialSimulator
%{

   public double[] getVariable(String speciesId) {
	int length = getVariableLength();
	double[] result = new double[length];
	__getVariable(speciesId, result);
	return result;
   }
   
   public double[] getGeometry(String compartmentId) {
	int length = getGeometryLength();
	double[] result = new double[length];
	__getGeometry(compartmentId, result);
	return result;
   }

   public double[] getX() {
     int length = getVariableLength();
	 double[] result = new double[length];
	 __getX(result);
	 return result;
   }
   
   public double[] getY() {
     int length = getVariableLength();
	 double[] result = new double[length];
	 __getY(result);
	 return result;
   }
   
   public double[] getZ() {
     int length = getVariableLength();
	 double[] result = new double[length];
	 __getZ(result);
	 return result;
   }

%}