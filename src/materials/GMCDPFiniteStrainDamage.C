#include "GMCDPFiniteStrainDamage.h"
#include <iostream>

namespace Marmot::Materials
{
using namespace Eigen;

using Vector6d = Matrix<double, 6, 1>;
using Matrix36 = Matrix<double, 3, 6>;
using Matrix66 = Matrix<double, 6, 6>;

Vector6d
stressToVoigt(const Matrix3d & stressTensor)
{
  Vector6d stress;
  // clang-format off
            stress << stressTensor( 0, 0 ), 
                      stressTensor( 1, 1 ), 
                      stressTensor( 2, 2 ), 
                      stressTensor( 0, 1 ),
                      stressTensor( 0, 2 ), 
                      stressTensor( 1, 2 );
  // clang-format on
  return stress;
}

Matrix3d
voigtToStrain(const Vector6d & voigt)
{
  Matrix3d strain;
  strain << voigt[0], voigt[3] / 2, voigt[4] / 2, voigt[3] / 2, voigt[1], voigt[5] / 2,
      voigt[4] / 2, voigt[5] / 2, voigt[2];
  return strain;
}

std::pair<Eigen::Vector3d, Eigen::Matrix<double, 3, 6>>
principalValuesAndDerivatives(const Eigen::Matrix<double, 6, 1> & S)
{
  // This is a fast implementation of the classical algorithm for determining
  // the principal components of a symmetric 3x3 Matrix in Voigt notation (off
  // diagonals expected with factor 1) as well as its respective derivatives

  constexpr double Pi = 3.141592653589793238463;

  using namespace Eigen;

  const static auto dS0_dS = (Vector6d() << 1, 0, 0, 0, 0, 0).finished();
  const static auto dS1_dS = (Vector6d() << 0, 1, 0, 0, 0, 0).finished();
  const static auto dS2_dS = (Vector6d() << 0, 0, 1, 0, 0, 0).finished();

  const static auto I = (Vector6d() << 1, 1, 1, 0, 0, 0).finished();

  const double p1 = S(3) * S(3) + S(4) * S(4) + S(5) * S(5);

  if (p1 <= 1e-16)
  { // matrix is already diagonal
    // clang-format off
                    const static auto dE_dS = ( Matrix36() << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ).finished();
                    return { S.head( 3 ), dE_dS };
    // clang-format on
  }

  const Vector6d dP1_dS = (Vector6d() << 0, 0, 0, 2 * S(3), 2 * S(4), 2 * S(5)).finished();

  const double q = S.head(3).sum() / 3;
  const Vector6d dQ_dS = I / 3;

  const double p2 = std::pow(S(0) - q, 2) + std::pow(S(1) - q, 2) + std::pow(S(2) - q, 2) + 2 * p1;
  const Vector6d dP2_dS = 2 * (S(0) - q) * (dS0_dS - dQ_dS) + 2 * (S(1) - q) * (dS1_dS - dQ_dS) +
                          2 * (S(2) - q) * (dS2_dS - dQ_dS) + 2 * dP1_dS;

  const double p = std::sqrt(p2 / 6);
  const Vector6d dP_dS = 1. / 12 * 1. / p * dP2_dS;

  const Vector6d B = 1. / p * (S - q * I);
  const Matrix66 dB_dS =
      -B * 1. / p * dP_dS.transpose() + 1. / p * (Matrix66::Identity() - I * dQ_dS.transpose());

  const double detB = B(0) * B(1) * B(2) + B(3) * B(4) * B(5) * 2 - B(2) * B(3) * B(3) -
                      B(1) * B(4) * B(4) - B(0) * B(5) * B(5);

  Vector6d dDetB_dB;
  dDetB_dB(0) = B(1) * B(2) - B(5) * B(5);
  dDetB_dB(1) = B(0) * B(2) - B(4) * B(4);
  dDetB_dB(2) = B(0) * B(1) - B(3) * B(3);
  dDetB_dB(3) = B(4) * B(5) * 2 - B(2) * B(3) * 2;
  dDetB_dB(4) = B(3) * B(5) * 2 - B(1) * B(4) * 2;
  dDetB_dB(5) = B(3) * B(4) * 2 - B(0) * B(5) * 2;

  const double r = detB * 1. / 2;
  const Vector6d dR_dS = 1. / 2 * dDetB_dB.transpose() * dB_dS;

  double phi;
  double dPhi_dR;
  if (r <= -1)
  {
    phi = Pi / 3;
    dPhi_dR = 0.0;
  }
  else if (r >= 1)
  {
    phi = 1.0;
    dPhi_dR = 0.0;
  }
  else
  {
    phi = std::acos(r) / 3;
    dPhi_dR = 1. / 3 * -1. / (std::sqrt(1 - r * r));
  }
  const Vector6d dPhi_dS = dPhi_dR * dR_dS;

  Vector3d e;
  Matrix36 dE_dS;

  e(0) = q + 2 * p * std::cos(phi);
  dE_dS.row(0) = dQ_dS + 2 * (dP_dS * std::cos(phi) - p * std::sin(phi) * dPhi_dS);

  const double phiShifted = phi + (2. / 3 * Pi);
  e(2) = q + 2 * p * std::cos(phiShifted);
  dE_dS.row(2) = dQ_dS + 2 * (dP_dS * std::cos(phiShifted) - p * std::sin(phiShifted) * dPhi_dS);

  e(1) = 3 * q - e(0) - e(2);
  dE_dS.row(1) = 3 * dQ_dS - dE_dS.row(0).transpose() - dE_dS.row(2).transpose();

  return {e, dE_dS};
}

GMCDPFiniteStrainDamage::GMCDPFiniteStrainDamage(const ModelParameters & modelParameters_)
  : modelParameters(modelParameters_)
{
}

GMCDPFiniteStrainDamage::DamageResponse
GMCDPFiniteStrainDamage::computeDamage(double alphaLocal,
                                       double alphaNonLocal,
                                       const Eigen::Matrix3d & deltaFp) const
{
  const auto [deltaAlphaLocal, dAlphaLocal_dDeltaFp] = computeDeltaAlphaLocal(deltaFp);

  const double alphaLocalNew = alphaLocal + deltaAlphaLocal;

  const auto [alphaWeighted_, dAlphaWeighted_dAlphaLocal_, dAlphaWeighted_dAlphaNonLocal_] =
      alphaWeighted(alphaLocalNew, alphaNonLocal);

  const auto [omega_, dOmega_dAlphaWeighted_] = this->omega(alphaWeighted_);

  return DamageResponse{.omega = omega_,
                        .alphaLocal = alphaLocalNew,
                        .dOmega_dAlphaLocal = dOmega_dAlphaWeighted_ * dAlphaWeighted_dAlphaLocal_,
                        .dOmega_dAlphaNonLocal =
                            dOmega_dAlphaWeighted_ * dAlphaWeighted_dAlphaNonLocal_,
                        .dAlphaLocal_dDeltaFp = dAlphaLocal_dDeltaFp};
}

GMCDPFiniteStrainDamage::ReducedPureNonLocalDamageResponse
GMCDPFiniteStrainDamage::computeDamage(double alphaLocal, double alphaNonLocal) const
{
  const auto [alphaWeighted_, dAlphaWeighted_dAlphaLocal_, dAlphaWeighted_dAlphaNonLocal_] =
      alphaWeighted(alphaLocal, alphaNonLocal);

  const auto [omega_, dOmega_dAlphaWeighted_] = this->omega(alphaWeighted_);

  return ReducedPureNonLocalDamageResponse{.omega = omega_,
                                           .dOmega_dAlphaNonLocal = dOmega_dAlphaWeighted_ *
                                                                    dAlphaWeighted_dAlphaNonLocal_};
}

std::pair<double, double>
GMCDPFiniteStrainDamage::omega(double alphaWeighted) const
{
  if (alphaWeighted <= 1e-15)
    return {0.0, 0.0};

  const double omega = 1 - std::exp(-alphaWeighted / modelParameters.softeningModulus);

  if (omega >= modelParameters.maxDamage)
    return {modelParameters.maxDamage, 0.0};

  else
    return {omega, 1 / modelParameters.softeningModulus * (-omega + 1)};
}

std::tuple<double, double, double>
GMCDPFiniteStrainDamage::alphaWeighted(double alphaLocal, double alphaNonLocal) const
{
  const double dAlphaWeighted_dAlphaLocal = 1 - modelParameters.weightingParameter;
  const double dAlphaWeighted_dAlphaNonLocal = modelParameters.weightingParameter;

  return {dAlphaWeighted_dAlphaNonLocal * alphaNonLocal + dAlphaWeighted_dAlphaLocal * alphaLocal,
          dAlphaWeighted_dAlphaLocal,
          dAlphaWeighted_dAlphaNonLocal};
}

std::pair<double, Eigen::Matrix3d>
GMCDPFiniteStrainDamage::computeDeltaAlphaLocal(const Eigen::Matrix3d & deltaFp) const
{

  const Matrix3d deltaCp = deltaFp.transpose() * deltaFp;

  // TODO: stress to voigt since principalsOfVoigt takes shear terms with factor 1 ... will be fixed
  // anyways!
  const Vector6d deltaCpVoigt = stressToVoigt(Eigen::Matrix3d(deltaCp.data()));

  const auto [deltaCpPrincipal, dSortedStrainPrincipal_dStrain] =
      principalValuesAndDerivatives(deltaCpVoigt);

  const Eigen::Vector3d dEPPrinc = deltaCpPrincipal.array().sqrt().log();

  const double dEPVol = dEPPrinc.sum();

  const auto [Rs_, dRs_dEPPrinc_] = this->Rs(dEPPrinc);
  const auto [xs_, dXs_dRs_] = this->xs(Rs_);

  // deltaAlphaLocal:
  const double deltaAlphaLocal = dEPVol / xs_;

  // deltaAlphaLocal_dDeltaFp:
  const Eigen::Vector3d dAlphaLocal_dEPPrinc =
      Eigen::Vector3d::Constant(1. / xs_) - dEPVol / (xs_ * xs_) * dXs_dRs_ * dRs_dEPPrinc_;

  const Eigen::Vector3d dAlphaLocal_dDeltaCpPrincipal =
      1. / 2 * dAlphaLocal_dEPPrinc.array() / deltaCpPrincipal.array();

  const Vector6d dAlphaLocal_dDeltaCpVoigt =
      dAlphaLocal_dDeltaCpPrincipal.transpose() * dSortedStrainPrincipal_dStrain;

  // Note voigtToStrain (instead of voigtToStress) since we are computing
  // Eigenvalues of the 'symmetric' part of Cp (which is always symmetric)
  // for the derivatives w.r.t. Cp the factor 0.5 accounts for the
  // correction!
  const Matrix3d dAlphaLocal_dDeltaCp(voigtToStrain(dAlphaLocal_dDeltaCpVoigt).data());

  const Matrix3d dAlphaLocal_dDeltaFp =
      deltaFp * (dAlphaLocal_dDeltaCp + dAlphaLocal_dDeltaCp.transpose());

  return {deltaAlphaLocal, dAlphaLocal_dDeltaFp};
}

std::pair<double, double>
GMCDPFiniteStrainDamage::xs(double Rs_) const
{
  if (Rs_ < 1)
    return {1 + modelParameters.As * Rs_ * Rs_, 2 * modelParameters.As * Rs_};
  else
    return {1 + modelParameters.As * (4 * sqrt(Rs_) - 3),
            2 * modelParameters.As * 1. / (sqrt(Rs_))};
}

std::pair<double, Eigen::Vector3d>
GMCDPFiniteStrainDamage::Rs(const Eigen::Vector3d & dEPP) const
{

  const double dEPVol = std::max(dEPP.sum(), 1e-16);

  Vector3d dEPVolNeg_dEPP = Vector3d::Zero();
  double dEPVolNeg = 0.0;

  for (int i = 0; i < 3; i++)
  {
    const double negDE = -dEPP(i);
    dEPVolNeg += negDE > 0 ? negDE : 0.0;
    dEPVolNeg_dEPP(i) = negDE > 0 ? -1 : 0.0;
  }

  const double Rs_ = dEPVolNeg / dEPVol;

  const Vector3d dRs_dEPPrinc_ =
      1. / dEPVol * dEPVolNeg_dEPP - Vector3d::Constant(dEPVolNeg / (dEPVol * dEPVol));

  return {Rs_, dRs_dEPPrinc_};
}
} // namespace Marmot::Materials
