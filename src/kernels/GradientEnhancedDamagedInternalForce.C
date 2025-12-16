#include <GradientEnhancedDamagedInternalForce.h>

// We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", GradientEnhancedDamagedInternalForce);

InputParameters
GradientEnhancedDamagedInternalForce::validParams()
{
  InputParameters params = InternalForce::validParams();
  params.addRequiredCoupledVar("nonlocal_damage", "The nonlocal damage field");
  return params;
}

GradientEnhancedDamagedInternalForce::GradientEnhancedDamagedInternalForce(
    const InputParameters & parameters)
  : InternalForce(parameters),
    _nonlocal_damage_var(coupled("nonlocal_damage")),
    _mat_omega(getMaterialPropertyByName<Real>("mat_omega")),
    _domega_dk(getMaterialPropertyByName<Real>("domega_dk"))
{
}

// TODO: This Jacobian needs to be fixed
Real
GradientEnhancedDamagedInternalForce::computeQpOffDiagJacobian(unsigned int jvar)
{

	/*	
	std::cout << "Info" << std::endl;
	std::cout << _domega_dk[_qp] << std::endl;
	std::cout << "--------" << std::endl;
  */

  if (jvar == _nonlocal_damage_var)
    return -InternalForce::computeQpResidual() * _domega_dk[_qp] * _phi[_j][_qp] /
           (1. - _mat_omega[_qp]);

  return InternalForce::computeQpOffDiagJacobian(jvar);
}
