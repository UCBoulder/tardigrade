#pragma once

#include "ALEKernel.h"
#include "FastorHelper.h"

using Arr3 = std::vector<Real>;
using Arr33 = std::vector<Arr3>;
using Arr333 = std::vector<Arr33>;
using Arr3333 = std::vector<Arr333>;

/**
 * Computes the classical 2.o Helmholtz like equation for nonlocal damage
 */
class GradientEnhancedMicromorphicDamage : public ALEKernel
{
public:
  static InputParameters validParams();

  GradientEnhancedMicromorphicDamage(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  Real computeQpJacobianDisplacement(unsigned int comp_j);
  Real computeQpJacobianMicroDisplacementGradient(unsigned int comp_j);
  Real computeQpJacobianNonlocalDamage();

  /// Base name of the material system that this kernel applies to
  const std::string _base_name;

  const MaterialProperty<Real> & _k_local;
  const MaterialProperty<Real> & _nonlocal_radius;
  /// Derivatives of the w.r.t. deformation gradient, micro rotations, material gradient of the micro rotations and the nonlocal damage driving field
  // const MaterialProperty< Arr33 > & _dk_local_dF;
  const MaterialProperty<Tensor33R> & _dk_local_dF;
  /* const MaterialProperty< Arr3 > & _dk_local_dPhi */
  /* const MaterialProperty< Arr33 > & _dk_local_dgrad_phi; */
  /* const MaterialProperty< Real > & _dk_local_dk; */

  /// Coupled displacement variables
  unsigned int _ndisp;
  /// Displacement variables IDs
  std::vector<unsigned int> _disp_var;

  /// Coupled micro rotation variables
  unsigned int _nmicro_disp_gradient;
  /// Micro rotation variables IDs
  std::vector<unsigned int> _micro_disp_gradient_var;

  /// The MOOSE variable number of the nonlocal damage variable
  unsigned int _nonlocal_damage_var;
};
