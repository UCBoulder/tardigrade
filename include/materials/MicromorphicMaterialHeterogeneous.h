/*!
==================================================================
|              MicromorphicMaterialHeterogeneous.h               |
==================================================================
*/

#ifndef MICROMORPHICMATERIALHETEROGENEOUS_H
#define MICROMORPHICMATERIALHETEROGENEOUS_H

#include "Material.h"
#include "MooseException.h"

#include "MicromorphicMaterial.h"

//Forward declarations
class MicromorphicMaterialHeterogeneous;

class MicromorphicMaterialHeterogeneous : public MicromorphicMaterial
{
    public:
        MicromorphicMaterialHeterogeneous(const InputParameters & parameters);

        static InputParameters validParams();

    protected:
        virtual void initQpStatefulProperties() override;
        virtual void computeQpProperties() override;
        virtual void defineMaterialParameters( std::vector< double > &parameters ) override;

    private:
        const std::vector<std::string> & _user_material_prop_names; //!< string mapping _user_material_prop_indices to GenericFunctionMaterial prop_names
        const std::vector<unsigned int> & _user_material_prop_indices; //!< string of indices specifying which values in material_fparameters will be replaced by GenericFunctionMaterial
        std::vector<const MaterialProperty<Real> *> _user_materials; //! array for storing unique material properties
};

#endif