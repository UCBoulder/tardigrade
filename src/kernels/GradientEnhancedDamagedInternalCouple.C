#include <GradientEnhancedDamagedInternalCouple.h>

registerMooseObject("tardigradeApp", GradientEnhancedDamagedInternalCouple);

InputParameters
GradientEnhancedDamagedInternalCouple::validParams()
{
  InputParameters params = InternalCouple::validParams();
  params.addRequiredCoupledVar("nonlocal_damage", "The nonlocal damage field");
  return params;
}

GradientEnhancedDamagedInternalCouple::GradientEnhancedDamagedInternalCouple(
    const InputParameters & parameters)
  : InternalCouple(parameters),
    _nonlocal_damage_var(coupled("nonlocal_damage")),
    _mat_omega(getMaterialPropertyByName<Real>("mat_omega")),
    _domega_dk(getMaterialPropertyByName<Real>(+"domega_dk"))
{
}

Real
GradientEnhancedDamagedInternalCouple::computeQpOffDiagJacobian(unsigned int jvar)
{

  if (jvar == _nonlocal_damage_var)
    return -InternalCouple::computeQpResidual() * _domega_dk[_qp] * _phi[_j][_qp] /
           (1. - _mat_omega[_qp]);

  return InternalCouple::computeQpOffDiagJacobian(jvar);
}
