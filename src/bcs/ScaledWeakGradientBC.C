//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ScaledWeakGradientBC.h"

registerMooseObject("MooseApp", ScaledWeakGradientBC);

InputParameters
ScaledWeakGradientBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addParam< Real >( "factor", 1.0, "The scaling factor for the BC" );
  params.addClassDescription(
      "Computes a boundary residual contribution consistent with the Diffusion Kernel. "
      "Does not impose a boundary condition; instead computes the boundary "
      "contribution corresponding to the current value of grad(u) and accumulates "
      "it in the residual vector.");
  return params;
}

ScaledWeakGradientBC::ScaledWeakGradientBC(const InputParameters & parameters) : IntegratedBC(parameters),
    _factor( getParam< Real >( "factor" ) ) {}

Real
ScaledWeakGradientBC::computeQpResidual()
{
  return _factor * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
}

Real
ScaledWeakGradientBC::computeQpJacobian()
{
  return _factor * (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp];
}
