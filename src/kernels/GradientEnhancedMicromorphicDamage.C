#include "GradientEnhancedMicromorphicDamage.h"

registerMooseObject("tardigradeApp", GradientEnhancedMicromorphicDamage);

InputParameters
GradientEnhancedMicromorphicDamage::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("2. order Helmholtz like PDE for describing nonlocal damage");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addRequiredCoupledVar("micro_displacement_gradient",
                               "The components of the micro displacement gradient");
  params.addRequiredCoupledVar("nonlocal_damage", "The nonlocal damage field");
  params.addRequiredParam<MaterialPropertyName>("nonlocal_radius",
                                                "Non local radius used in the damage model.");
  return params;
}

GradientEnhancedMicromorphicDamage::GradientEnhancedMicromorphicDamage(
    const InputParameters & parameters)
  : ALEKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _k_local(getMaterialPropertyByName<Real>(_base_name + "k_local")),
    _nonlocal_radius(getMaterialProperty<Real>(_base_name + "nonlocal_radius")),
    //_dk_local_dF( getMaterialPropertyByName< Arr33 >( _base_name + "dk_local_dF" ) ),
    _dk_local_dF(getMaterialPropertyByName<Tensor33R>(_base_name + "dk_local_dF")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _nmicro_disp_gradient(coupledComponents("micro_displacement_gradient")),
    _micro_disp_gradient_var(_nmicro_disp_gradient),
    _nonlocal_damage_var(coupled("nonlocal_damage"))
{
  if (_ndisp != 3 || _nmicro_disp_gradient != 9)
    mooseError("Gradient-enhanced micromorphic kernels are implemented only for 3D!");

  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
  for (unsigned int i = 0; i < _nmicro_disp_gradient; ++i)
    _micro_disp_gradient_var[i] = coupled("micro_displacement_gradient", i);
}

Real
GradientEnhancedMicromorphicDamage::computeQpResidual()
{
  return std::pow(_nonlocal_radius[_qp], 2) * _grad_test[_i][_qp] * _grad_u[_qp] +
         _test[_i][_qp] * (_u[_qp] - _k_local[_qp]);
}

Real
GradientEnhancedMicromorphicDamage::computeQpJacobian()
{
  return computeQpJacobianNonlocalDamage();
}

Real
GradientEnhancedMicromorphicDamage::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (unsigned int j = 0; j < _ndisp; ++j)
    if (jvar == _disp_var[j])
      return computeQpJacobianDisplacement(j);

  for (unsigned int j = 0; j < _nmicro_disp_gradient; ++j)
    if (jvar == _micro_disp_gradient_var[j])
      return computeQpJacobianMicroDisplacementGradient(j);

  return 0.0;
}

Real
GradientEnhancedMicromorphicDamage::computeQpJacobianDisplacement(unsigned int comp_j)
{
  Real df_du_j = 0;

  for (int K = 0; K < 3; K++)
    df_du_j += -1 * (_dk_local_dF[_qp])(comp_j, K) * _grad_phi[_j][_qp](K);

  return _test[_i][_qp] * df_du_j;
}

Real
GradientEnhancedMicromorphicDamage::computeQpJacobianMicroDisplacementGradient(unsigned int comp_j)
{
  return 0;
}

Real
GradientEnhancedMicromorphicDamage::computeQpJacobianNonlocalDamage()
{
  return std::pow(_nonlocal_radius[_qp], 2) * _grad_test[_i][_qp] * _grad_phi[_j][_qp] +
         _test[_i][_qp] * _phi[_j][_qp];
}
