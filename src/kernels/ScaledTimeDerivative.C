//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ScaledTimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"

#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", ScaledTimeDerivative);

InputParameters
ScaledTimeDerivative::validParams()
{
  InputParameters params = TimeDerivative::validParams();
  params.addParam<Real>("factor",1,"The scale factor for the time derivative");
  return params;
}

ScaledTimeDerivative::ScaledTimeDerivative(const InputParameters & parameters)
  : TimeDerivative(parameters), _factor(getParam<Real>("factor"))
{
}

Real
ScaledTimeDerivative::computeQpResidual()
{
  return _factor * _test[_i][_qp] * _u_dot[_qp];
}

Real
ScaledTimeDerivative::computeQpJacobian()
{
  return _factor * _test[_i][_qp] * _phi[_j][_qp] * _du_dot_du[_qp];
}
