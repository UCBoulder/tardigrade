#pragma once

#include <Eigen/Core>

namespace Marmot::Materials
{

class GMCDPFiniteStrainDamage
{

private:
  std::tuple<double, double, double> alphaWeighted(double alphaLocal, double alphaNonLocal) const;
  std::pair<double, double> omega(double alphaWeighted) const;

  std::pair<double, double> xs(double Rs) const;
  std::pair<double, Eigen::Vector3d> Rs(const Eigen::Vector3d & dHenckyPlasticPrincipal) const;
  std::pair<double, Eigen::Matrix3d> computeDeltaAlphaLocal(const Eigen::Matrix3d & deltaFp) const;

public:
  struct DamageResponse
  {
    double omega;
    double alphaLocal;

    double dOmega_dAlphaLocal;
    double dOmega_dAlphaNonLocal;

    Eigen::Matrix3d dAlphaLocal_dDeltaFp;
  };

  struct ReducedPureNonLocalDamageResponse
  {
    double omega;
    double dOmega_dAlphaNonLocal;
  };

  struct ModelParameters
  {
    double As;
    double softeningModulus;
    double weightingParameter;
    double maxDamage;
  };

  const ModelParameters modelParameters;

  GMCDPFiniteStrainDamage(const ModelParameters & modelParameters);

  DamageResponse
  computeDamage(double alphaLocal, double alphaNonLocal, const Eigen::Matrix3d & deltaFp) const;

  ReducedPureNonLocalDamageResponse computeDamage(double alphaLocal, double alphaNonLocal) const;
};
} // namespace Marmot::Materials
