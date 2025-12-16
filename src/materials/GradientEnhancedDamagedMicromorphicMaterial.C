#include "GradientEnhancedDamagedMicromorphicMaterial.h"
#include "GMCDPFiniteStrainDamage.h"

// MOOSE includes
#include "Function.h"

#include <Fastor/Fastor.h>
using namespace Fastor;

registerMooseObject("tardigradeApp", GradientEnhancedDamagedMicromorphicMaterial);

using vecReal = std::vector<double>;
using vecVecReal = std::vector<vecReal>;

template <int nRows, int nCols>
Eigen::Matrix<Real, nRows, nCols, Eigen::RowMajor>
copyVecVecToMatrix(const vecVecReal & vecvec)
{
  Eigen::Matrix<Real, nRows, nCols, Eigen::RowMajor> theMatrix;

  for (int i = 0; i < nRows; i++)
    theMatrix.row(i) = Eigen::Map<const Eigen::Matrix<Real, nCols, 1>>(vecvec[i].data());
  return theMatrix;
};

template <int nRows, int nCols>
auto copyMatrixToVecVec = [](const Eigen::Matrix<Real, nRows, nCols, Eigen::RowMajor> & theMatrix) {
  vecVecReal theVecVec(nRows);
  for (int i = 0; i < nRows; i++)
  {
    theVecVec[i].resize(nCols);
    Eigen::Map<Eigen::Matrix<Real, nCols, 1>>(theVecVec[i].data()) = theMatrix.row(i);
  }
  return theVecVec;
};

template <int nRows, int nCols>
vecVecReal
copyMatrixToVecVec2(const Eigen::Matrix<Real, nRows, nCols, Eigen::RowMajor> & theMatrix)
{
  vecVecReal theVecVec(nRows);
  for (int i = 0; i < nRows; i++)
  {
    theVecVec[i].resize(nCols);
    Eigen::Map<Eigen::Matrix<Real, nCols, 1>>(theVecVec[i].data()) = theMatrix.row(i);
  }
  return theVecVec;
};

template <int nRows, int nCols>
vecVecReal
copyMatrixToVecVec3(const Eigen::Matrix<Real, nRows, nCols, Eigen::RowMajor> & theMatrix)
{
  vecVecReal theVecVec(nRows);
  auto ptr = theMatrix.data();
  for (auto & vec : theVecVec)
  {
    vec = std::vector<Real>(ptr, ptr + theMatrix.cols());
    ptr += theMatrix.cols();
  }
  return theVecVec;
};

InputParameters
GradientEnhancedDamagedMicromorphicMaterial::validParams()
{
  InputParameters params = Material::validParams();

  // Vectors of material properties
  params.addRequiredParam<std::vector<Real>>(
      "material_fparameters",
      "The vector of floating point material parameters required for the stiffness matrices");

  params.addRequiredParam<std::string>("model_name", "The material model name");

  params.addParam<int>(
      "number_ADD_DOF",
      0,
      "The number of additional degrees of freedom beyond u and phi in the problem");
  params.addParam<int>("number_ADD_TERMS",
                       0,
                       "The number of additional balance equations being solved beyond the balance "
                       "of linear momentum and first moment of momentum");
  params.addParam<int>(
      "number_ADD_JACOBIANS",
      0,
      "The number of additional jacobians being provided beyond those of the stress measures");

  // params.addParam<bool>(
  //    "MMS", false, "Flag for whether to compute the method of manufactured solutions");

  params.addRequiredParam<std::vector<Real>>(
      "gradient_enhanced_damage_fparameters",
      "The vector of floating point material parameters for the gradient-enhanced damage");

  params.addRequiredCoupledVar("displacements", "The 3 displacement components");
  params.addRequiredCoupledVar("micro_displacement_gradient",
                               "The 9 components of the micro displacement gradient");
  params.addRequiredCoupledVar("nonlocal_damage", "The nonlocal damage variable");

  // The state variable array
  params.addParam<int>("number_SDVS", 0, "The number of solution-dependent state variables");

  // Functions for method of manufactured solutions
  params.addParam<FunctionName>("u1_fxn", "0", "The function for displacement in the 1 direction.");

  params.addParam<FunctionName>("u2_fxn", "0", "The function for displacement in the 2 direction.");

  params.addParam<FunctionName>("u3_fxn", "0", "The function for displacement in the 3 direction.");
  params.addParam<FunctionName>(
      "phi_11_fxn", "0", "The function for the 11 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_22_fxn", "0", "The function for the 22 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_33_fxn", "0", "The function for the 33 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_23_fxn", "0", "The function for the 23 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_13_fxn", "0", "The function for the 13 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_12_fxn", "0", "The function for the 12 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_32_fxn", "0", "The function for the 32 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_31_fxn", "0", "The function for the 31 component of the phi tensor.");
  params.addParam<FunctionName>(
      "phi_21_fxn", "0", "The function for the 21 component of the phi tensor.");

  params.addRequiredParam<std::vector<std::string>>(
      "user_material_prop_names", "Names of material properties obtained from file");
  params.addRequiredParam<std::vector<unsigned int>>(
      "user_material_prop_indices", "Index of material properties obtained from file");

  return params;
}

void
GradientEnhancedDamagedMicromorphicMaterial::initQpStatefulProperties()
{
  // Initialize the internal state variable array
  _SDVS[_qp] = std::vector<Real>(_n_SDVS, 0);
  _ge_damage_statevars[_qp] = std::vector<Real>(2, 0);

  return;
}

GradientEnhancedDamagedMicromorphicMaterial::GradientEnhancedDamagedMicromorphicMaterial(
    const InputParameters & parameters)
  : Material(parameters),

    // Declare that this material is going to provide Eigen matrices containing the PK2 stress and
    // jacobians that Kernels can use.
    _fparams(getParam<std::vector<Real>>("material_fparameters")),
    _n_ADD_DOF(getParam<int>("number_ADD_DOF")),
    _n_ADD_TERMS(getParam<int>("number_ADD_TERMS")),
    _n_ADD_JACOBIANS(getParam<int>("number_ADD_JACOBIANS")),
    _model_name(getParam<std::string>("model_name")),
    // _MMS(getParam<bool>("MMS")),
    _n_SDVS(getParam<int>("number_SDVS")),

    _ge_damage_parameters(getParam<std::vector<Real>>("gradient_enhanced_damage_fparameters")),
    _ge_damage_statevars(declareProperty<std::vector<double>>("ge_damage_statevars")),
    _ge_damage_statevars_old(getMaterialPropertyOld<std::vector<double>>("ge_damage_statevars")),

    _grad_disp(coupledGradients("displacements")),
    _grad_disp_old(coupledGradientsOld("displacements")),

    _micro_disp_gradient(coupledValues("micro_displacement_gradient")),
    _micro_disp_gradient_old(coupledValuesOld("micro_displacement_gradient")),

    _grad_micro_disp_gradient(coupledGradients("micro_displacement_gradient")),
    _grad_micro_disp_gradient_old(coupledGradientsOld("micro_displacement_gradient")),

    _k(coupledValue("nonlocal_damage")),
    _k_local(declareProperty<Real>("k_local")),
    _mat_omega(declareProperty<Real>("mat_omega")),
    _dk_local_dF(declareProperty<Tensor33R>("dk_local_dF")),
    _domega_dk(declareProperty<Real>("domega_dk")),

    _deformation_gradient(declareProperty<std::vector<double>>("MM_deformation_gradient")),
    _micro_deformation(declareProperty<std::vector<double>>("micro_deformation")),
    _gradient_micro_displacement(
        declareProperty<std::vector<std::vector<double>>>("gradient_micro_displacement")),
    _cauchy(declareProperty<std::vector<double>>("cauchy")),
    _s(declareProperty<std::vector<double>>("s")),
    _m(declareProperty<std::vector<double>>("m")),
    _PK2(declareProperty<std::vector<double>>("PK2")),
    _SIGMA(declareProperty<std::vector<double>>("SIGMA")),
    _M(declareProperty<std::vector<double>>("M")),
    _DPK2Dgrad_u(declareProperty<std::vector<std::vector<double>>>("DPK2Dgrad_u")),
    _DPK2Dphi(declareProperty<std::vector<std::vector<double>>>("DPK2Dphi")),
    _DPK2Dgrad_phi(declareProperty<std::vector<std::vector<double>>>("DPK2Dgrad_phi")),
    _DSIGMADgrad_u(declareProperty<std::vector<std::vector<double>>>("DSIGMADgrad_u")),
    _DSIGMADphi(declareProperty<std::vector<std::vector<double>>>("DSIGMADphi")),
    _DSIGMADgrad_phi(declareProperty<std::vector<std::vector<double>>>("DSIGMADgrad_phi")),
    _DMDgrad_u(declareProperty<std::vector<std::vector<double>>>("DMDgrad_u")),
    _DMDphi(declareProperty<std::vector<std::vector<double>>>("DMDphi")),
    _DMDgrad_phi(declareProperty<std::vector<std::vector<double>>>("DMDgrad_phi")),
    _ADD_TERMS(declareProperty<std::vector<std::vector<double>>>("ADD_TERMS")),
    _ADD_JACOBIANS(declareProperty<std::vector<std::vector<std::vector<double>>>>("ADD_JACOBIANS")),
    _SDVS(declareProperty<std::vector<double>>("SDVS")),
    _old_SDVS(getMaterialPropertyOld<std::vector<double>>("SDVS")),
    _u1_fxn(getFunction("u1_fxn")),
    _u2_fxn(getFunction("u2_fxn")),
    _u3_fxn(getFunction("u3_fxn")),
    _phi_11_fxn(getFunction("phi_11_fxn")),
    _phi_22_fxn(getFunction("phi_22_fxn")),
    _phi_33_fxn(getFunction("phi_33_fxn")),
    _phi_23_fxn(getFunction("phi_23_fxn")),
    _phi_13_fxn(getFunction("phi_13_fxn")),
    _phi_12_fxn(getFunction("phi_12_fxn")),
    _phi_32_fxn(getFunction("phi_32_fxn")),
    _phi_31_fxn(getFunction("phi_31_fxn")),
    _phi_21_fxn(getFunction("phi_21_fxn")),
    _PK2_MMS(declareProperty<std::vector<double>>("PK2_MMS")),
    _SIGMA_MMS(declareProperty<std::vector<double>>("SIGMA_MMS")),
    _M_MMS(declareProperty<std::vector<double>>("M_MMS")),
    _ADD_TERMS_MMS(declareProperty<std::vector<std::vector<double>>>("ADD_TERMS_MMS")),

    _user_material_prop_names(getParam<std::vector<std::string>>("user_material_prop_names")),
    _user_material_prop_indices(getParam<std::vector<unsigned int>>("user_material_prop_indices"))
{
  /*!
  ==============================
  |    GradientEnhancedDamagedMicromorphicMaterial    |
  ==============================
  The constructor for GradientEnhancedDamagedMicromorphicMaterial.
  */
  // Get the number of user material property names/indices
  unsigned int prop_names = _user_material_prop_names.size();
  unsigned int prop_indices = _user_material_prop_indices.size();

  // Resize the user_materials vector
  _user_materials.resize(prop_names);

  // Raise an error if the numbers of names and indices are not the same
  if (prop_names != prop_indices)
    mooseError("Number of prop names and prop indices do not match");

  // Get the corresponding material properties
  for (unsigned int i = 0; i < prop_names; i++)
    _user_materials[i] = &getMaterialProperty<Real>(_user_material_prop_names[i]);
}

void
GradientEnhancedDamagedMicromorphicMaterial::computeQpProperties()
{
  /*!
  =============================
  |    computeQpProperties    |
  =============================
  Evaluate the constitutive model at the quadrature point.
  */

  // Define required DOF values for the
  // balance of linear momentum and first moment of momentum
  double __grad_u[3][3], __old_grad_u[3][3];
  double __phi[9], __old_phi[9];
  double __grad_phi[9][3], __old_grad_phi[9][3];

  RealVectorValue tmp_grad;

  // Copy over the gradient of u
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      __grad_u[i][j] = (*_grad_disp[i])[_qp](j);

  // Copy over the old gradient of u
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      __old_grad_u[i][j] = (*_grad_disp_old[i])[_qp](j);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      __phi[i + j * 3] = (*_micro_disp_gradient[i + j * 3])[_qp];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      __old_phi[i + j * 3] = (*_micro_disp_gradient_old[i + j * 3])[_qp];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        __grad_phi[i + j * 3][k] = (*_grad_micro_disp_gradient[i + j * 3])[_qp](k);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        __old_grad_phi[i + j * 3][k] = (*_grad_micro_disp_gradient_old[i + j * 3])[_qp](k);

  // Copy the deformation gradient
  _deformation_gradient[_qp].resize(9);
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      _deformation_gradient[_qp][3 * i + j] = __grad_u[i][j];
    }
  }
  _deformation_gradient[_qp][0] += 1;
  _deformation_gradient[_qp][4] += 1;
  _deformation_gradient[_qp][8] += 1;

  // Copy the micro-deformation
  _micro_deformation[_qp].resize(9);
  for (int i = 0; i < 9; i++)
  {
    _micro_deformation[_qp][i] = __phi[i];
  }
  _micro_deformation[_qp][0] += 1;
  _micro_deformation[_qp][4] += 1;
  _micro_deformation[_qp][8] += 1;

  // Copy the gradient of the micro-displacement
  _gradient_micro_displacement[_qp].resize(9);
  for (int i = 0; i < 9; i++)
  {
    _gradient_micro_displacement[_qp][i].resize(3);
    for (int j = 0; j < 3; j++)
    {
      _gradient_micro_displacement[_qp][i][j] = __grad_phi[i][j];
    }
  }

  // Evaluate the stresses and their jacobians
  auto & factory = micromorphic_material_library::MaterialFactory::Instance();
  auto material = factory.GetMaterial(_model_name);

  // Extract the time
  std::vector<double> time;
  time.resize(2);

  time[0] = _t;
  time[1] = _dt;

  /* //TODO: Remove these hardcoded values */
  if (_n_ADD_DOF > 0)
  {
    mooseError("GradientEnhancedDamagedMicromorphicMaterial does not support additional degrees of "
               "freedom");
  }

  std::vector<double> ADD_DOF, old_ADD_DOF;
  ADD_DOF.resize(_n_ADD_DOF);
  old_ADD_DOF.resize(_n_ADD_DOF);

  std::vector<std::vector<double>> ADD_grad_DOF, old_ADD_grad_DOF;
  ADD_grad_DOF.resize(_n_ADD_DOF);
  old_ADD_grad_DOF.resize(_n_ADD_DOF);

  // END of hardcoded values

  // Set the sizes of the additional term vectors
  _ADD_TERMS[_qp].resize(_n_ADD_TERMS);
  _ADD_JACOBIANS[_qp].resize(_n_ADD_JACOBIANS);

  // Evaluate the model
  std::string output_message;

#ifdef DEBUG_MODE
  std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double>>>> debug;
#endif

  // Set the state variables to the previously converged values
  _SDVS[_qp] = _old_SDVS[_qp];

  // These three lines are just for sizing
  // effPK2, effSigma, and effM are re-written in tardigrade
  vecReal effPK2 = _PK2[_qp];
  vecReal effSigma = _SIGMA[_qp];
  vecReal effM = _M[_qp];

  vecVecReal dEffPK2_dGrad_u = _DPK2Dgrad_u[_qp], dEffPK2_dPhi = _DPK2Dphi[_qp],
             dEffPK2_dGrad_phi = _DPK2Dgrad_phi[_qp], dEffSigma_dGrad_u = _DSIGMADgrad_u[_qp],
             dEffSigma_dPhi = _DSIGMADphi[_qp], dEffSigma_dGrad_phi = _DSIGMADgrad_phi[_qp],
             dEffM_dGrad_u = _DMDgrad_u[_qp], dEffM_dPhi = _DMDphi[_qp],
             dEffM_dGrad_phi = _DMDgrad_phi[_qp];

  auto q_fparams = _fparams;
  // Rewrite the corresponding indices with the material values at each Gauss point
  for (unsigned int i = 0; i < _user_material_prop_indices.size(); i++)
    q_fparams[_user_material_prop_indices[i]] = (*_user_materials[i])[_qp];

  int errorCode =
      material->evaluate_model(time,
                               q_fparams,
                               __grad_u,
                               __phi,
                               __grad_phi,
                               __old_grad_u,
                               __old_phi,
                               __old_grad_phi,
                               _SDVS[_qp],
                               ADD_DOF,
                               ADD_grad_DOF,
                               old_ADD_DOF,
                               old_ADD_grad_DOF,
                               effPK2,
                               effSigma,
                               effM,
                               dEffPK2_dGrad_u,
                               dEffPK2_dPhi,
                               dEffPK2_dGrad_phi,
                               dEffSigma_dGrad_u,
                               dEffSigma_dPhi,
                               dEffSigma_dGrad_phi,
                               dEffM_dGrad_u,
                               dEffM_dPhi,
                               dEffM_dGrad_phi,
                               /* _PK2[_qp],           _SIGMA[_qp],        _M[_qp], */
                               /* _DPK2Dgrad_u[_qp],   _DPK2Dphi[_qp],     _DPK2Dgrad_phi[_qp], */
                               /* _DSIGMADgrad_u[_qp], _DSIGMADphi[_qp],   _DSIGMADgrad_phi[_qp], */
                               /* _DMDgrad_u[_qp],     _DMDphi[_qp],       _DMDgrad_phi[_qp], */
                               _ADD_TERMS[_qp],
                               _ADD_JACOBIANS[_qp],
                               output_message
#ifdef DEBUG_MODE
                               ,
                               debug
#endif
      );

  if (errorCode == 1)
  {
    std::string error_message =
        "Convergence not achieved in material model. Requesting timestep cutback.\n";
    error_message += output_message;
    mooseException(error_message.c_str());
  }

  if (errorCode == 2)
  {
    std::string error_message = "FATAL ERROR IN MICROMORPHIC MATERIAL MODEL\n";
    error_message += output_message;
    mooseError(error_message.c_str());
  }

  const double * ptr_Fp_np = &_SDVS[_qp][0];
  const double * ptr_Fp_n = &_old_SDVS[_qp][0];

  const auto & _dFp_dF = _ADD_JACOBIANS[_qp][0];

  // we don't use Eigen's default ColMajor layout but rather RowMajor for convenient interfacing
  // with STL
  using RMMatrix3d = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
  using RMMatrix9d = Eigen::Matrix<double, 9, 9, Eigen::RowMajor>;
  using RMMatrix9_27d = Eigen::Matrix<double, 9, 27, Eigen::RowMajor>;
  using RMMatrix27d = Eigen::Matrix<double, 27, 27, Eigen::RowMajor>;
  using RMMatrix27_9d = Eigen::Matrix<double, 27, 9, Eigen::RowMajor>;
  using RMMapMatrix3d = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>;
  using RMMapConstMatrix3d = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>;
  using Vector9d = Eigen::Matrix<double, 9, 1>;
  using Vector27d = Eigen::Matrix<double, 27, 1>;
  using MapVector9d = Eigen::Map<Vector9d>;
  using MapVector27d = Eigen::Map<Vector27d>;
  using MapConstVector9d = Eigen::Map<const Eigen::Matrix<double, 9, 1>>;
  using RMMapMatrix9d = Eigen::Map<Eigen::Matrix<double, 9, 9, Eigen::RowMajor>>;
  using RMMapConstMatrix9d = Eigen::Map<const Eigen::Matrix<double, 9, 9, Eigen::RowMajor>>;

  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  //                                   +
  //             DAMAGE PART           +
  //                                   +
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++

  const auto & As = _ge_damage_parameters[0];
  const auto & softeningModulus = _ge_damage_parameters[1];
  const auto & weightingParameter = _ge_damage_parameters[2];
  const auto & maxDamage = _ge_damage_parameters[3];

  const auto & alphaDLocalOld = _ge_damage_statevars_old[_qp][0];
  auto & alphaDLocal = _ge_damage_statevars[_qp][0];

  const auto & omegaOld = _ge_damage_statevars_old[_qp][1];
  auto & omega = _ge_damage_statevars[_qp][1];

  double delta_alphaDLocal;
  Tensor33R Fastor_dAlphaLocal_dF;

  if (_SDVS[_qp][50] < As)
  {
    delta_alphaDLocal = 0.;
    Fastor_dAlphaLocal_dF = 0.;
  }
  else
  {
    const auto dFp_np_dF = copyVecVecToMatrix<9, 9>(_dFp_dF);

    const RMMatrix3d Fp_np = RMMapConstMatrix3d(ptr_Fp_np) + RMMatrix3d::Identity();
    const RMMatrix3d Fp_n = RMMapConstMatrix3d(ptr_Fp_n) + RMMatrix3d::Identity();

    const Tensor33R Fastor_Fp(Fp_np.data());
    const Tensor33R Fastor_Fp_n(Fp_n.data());

    const RMMatrix3d invFp_n = Fp_n.inverse();
    const RMMatrix3d deltaFp = Fp_np * invFp_n;

    const Tensor33R Fastor_deltaFp(deltaFp.data());

    Tensor33R Fastor_I;
    Fastor_I.eye();

    /*
    const Tensor33R Fastor_Ep =
        0.5 * (Fastor::einsum<Index<0, 1>, Index<0, 2>>(Fastor_Fp, Fastor_Fp) - Fastor_I);

    Tensor33R Fastor_dev_Ep = Fastor_Ep - 1. / 3. * Fastor::trace(Fastor_Ep) * Fastor_I;
    double norm_dev_Ep = Fastor::norm(Fastor_dev_Ep);

    Tensor3333R Fastor_dEp_dFp;
    for (int i1 = 0; i1 < 3; i1++)
      for (int j1 = 0; j1 < 3; j1++)
        for (int k1 = 0; k1 < 3; k1++)
          for (int l1 = 0; l1 < 3; l1++)
            Fastor_dEp_dFp(i1, j1, k1, l1) =
                0.5 * (Fastor_Fp(k1, j1) * Fastor_I(i1, l1) + Fastor_Fp(k1, i1) * Fastor_I(j1, l1));

    Tensor33R Fastor_Fp_n_inv = Fastor::inverse(Fastor_Fp_n);

    Tensor3333R Fastor_dDeltaFp_dF =
        Fastor::einsum<Index<0, 1, 2, 3>, Index<1, 4>, OIndex<0, 4, 2, 3>>(Fastor_dFp_dF,
                                                                           Fastor_Fp_n_inv);
    */

    const Tensor33R Fastor_FpT_Fp = Fastor::einsum<Index<0, 1>, Index<0, 2>>(Fastor_Fp, Fastor_Fp);
    const Tensor33R Fastor_FpT_Fp_n =
        Fastor::einsum<Index<0, 1>, Index<0, 2>>(Fastor_Fp_n, Fastor_Fp_n);

    const Tensor33R Fastor_deltaEp = 0.5 * (Fastor_FpT_Fp - Fastor_FpT_Fp_n);

    double norm_deltaEp = Fastor::norm(Fastor_deltaEp);

    if (MooseUtils::absoluteFuzzyEqual(norm_deltaEp, 0))
      norm_deltaEp = libMesh::TOLERANCE * libMesh::TOLERANCE;

    Tensor33R Fastor_dAlphaLocal_ddeltaEp = Fastor_deltaEp / norm_deltaEp;

    Tensor3333R Fastor_ddeltaEp_dFp =
        0.5 * (Fastor::einsum<Index<2, 0>, Index<1, 3>, OIndex<0, 1, 2, 3>>(Fastor_Fp, Fastor_I) +
               Fastor::einsum<Index<2, 1>, Index<0, 3>, OIndex<0, 1, 2, 3>>(Fastor_Fp, Fastor_I));

    Tensor3333R Fastor_dFp_dF;
    for (int i1 = 0; i1 < 3; i1++)
      for (int j1 = 0; j1 < 3; j1++)
        for (int k1 = 0; k1 < 3; k1++)
          for (int l1 = 0; l1 < 3; l1++)
            Fastor_dFp_dF(i1, j1, k1, l1) = dFp_np_dF(3 * i1 + j1, 3 * k1 + l1);

    delta_alphaDLocal = norm_deltaEp;

    Fastor_dAlphaLocal_dF =
        Fastor::einsum<Index<0, 1>, Index<0, 1, 2, 3>, Index<2, 3, 4, 5>, OIndex<4, 5>>(
            Fastor_dAlphaLocal_ddeltaEp, Fastor_ddeltaEp_dFp, Fastor_dFp_dF);
  }

  alphaDLocal = alphaDLocalOld + delta_alphaDLocal;

  _k_local[_qp] = alphaDLocal;
  _dk_local_dF[_qp] = Fastor_dAlphaLocal_dF;

  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // Feedback Part
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++

  // omega = damageResult.omega;
  // const Vector9d dOmega_dF = damageResult.dOmega_dAlphaLocal * dAlphaLocal_dF;
  //_domega_dk[_qp] = damageResult.dOmega_dAlphaNonLocal;
  // Vector9d dOmega_dDeformationGradient = damageResult.dOmega_dAlphaLocal * dAlphaLocal_dDeltaFp;

  const Real alphaDNonLocal = _k[_qp];

  omega = 1 - std::exp(-alphaDNonLocal / softeningModulus);

  if (omega > maxDamage)
    omega = maxDamage;

  // storing the value of omega to be used in the kernel
  _mat_omega[_qp] = omega;
  _domega_dk[_qp] = 1 / softeningModulus * std::exp(-alphaDNonLocal / softeningModulus);

  Vector9d dOmega_dDeformationGradient = Vector9d::Zero();

  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++
  // +++++++++++++++++++++++++++++++++++

  MapVector9d _effPK2(&effPK2[0]);
  MapVector9d _effSigma(&effSigma[0]);
  MapVector27d _effM(&effM[0]);

  const auto _dEffPK2_dGrad_u = copyVecVecToMatrix<9, 9>(dEffPK2_dGrad_u);
  const auto _dEffPK2_dPhi = copyVecVecToMatrix<9, 9>(dEffPK2_dPhi);
  const auto _dEffPK2_dGrad_phi = copyVecVecToMatrix<9, 27>(dEffPK2_dGrad_phi);

  const auto _dEffSigma_dGrad_u = copyVecVecToMatrix<9, 9>(dEffSigma_dGrad_u);
  const auto _dEffSigma_dPhi = copyVecVecToMatrix<9, 9>(dEffSigma_dPhi);
  const auto _dEffSigma_dGrad_phi = copyVecVecToMatrix<9, 27>(dEffSigma_dGrad_phi);

  const auto _dEffM_dGrad_u = copyVecVecToMatrix<27, 9>(dEffM_dGrad_u);
  const auto _dEffM_dPhi = copyVecVecToMatrix<27, 9>(dEffM_dPhi);
  const auto _dEffM_dGrad_phi = copyVecVecToMatrix<27, 27>(dEffM_dGrad_phi);

  const Vector9d PK2 = (1 - omega) * _effPK2;
  const Vector9d Sigma = (1 - omega) * _effSigma;
  const Vector27d M = (1 - omega) * _effM;

  const RMMatrix9d dPK2_dGrad_u =
      (1 - omega) * _dEffPK2_dGrad_u - _effPK2 * dOmega_dDeformationGradient.transpose();
  const RMMatrix9d dPK2_dPhi = (1 - omega) * _dEffPK2_dPhi;
  const RMMatrix9_27d dPK2_dGrad_phi = (1 - omega) * _dEffPK2_dGrad_phi;

  const RMMatrix9d dSigma_dGrad_u =
      (1 - omega) * _dEffSigma_dGrad_u - _effSigma * dOmega_dDeformationGradient.transpose();
  const RMMatrix9d dSigma_dPhi = (1 - omega) * _dEffSigma_dPhi;
  const RMMatrix9_27d dSigma_dGrad_phi = (1 - omega) * _dEffSigma_dGrad_phi;

  const RMMatrix27_9d dM_dGrad_u =
      (1 - omega) * _dEffM_dGrad_u - _effM * dOmega_dDeformationGradient.transpose();
  const RMMatrix27_9d dM_dPhi = (1 - omega) * _dEffM_dPhi;
  const RMMatrix27d dM_dGrad_phi = (1 - omega) * _dEffM_dGrad_phi;

  _PK2[_qp].resize(9);
  _SIGMA[_qp].resize(9);
  _M[_qp].resize(27);
  MapVector9d(_PK2[_qp].data()) = PK2.transpose();
  MapVector9d(_SIGMA[_qp].data()) = Sigma.transpose();
  MapVector27d(_M[_qp].data()) = M.transpose();

  auto inter_test = copyMatrixToVecVec<9, 9>(dPK2_dGrad_u);
  _DPK2Dgrad_u[_qp] = inter_test;

  _DPK2Dphi[_qp] = copyMatrixToVecVec<9, 9>(dPK2_dPhi);
  _DPK2Dgrad_phi[_qp] = copyMatrixToVecVec<9, 27>(dPK2_dGrad_phi);

  _DSIGMADgrad_u[_qp] = copyMatrixToVecVec<9, 9>(dSigma_dGrad_u);
  _DSIGMADphi[_qp] = copyMatrixToVecVec<9, 9>(dSigma_dPhi);
  _DSIGMADgrad_phi[_qp] = copyMatrixToVecVec<9, 27>(dSigma_dGrad_phi);

  _DMDgrad_u[_qp] = copyMatrixToVecVec<27, 9>(dM_dGrad_u);
  _DMDphi[_qp] = copyMatrixToVecVec<27, 9>(dM_dPhi);
  _DMDgrad_phi[_qp] = copyMatrixToVecVec<27, 27>(dM_dGrad_phi);
}
