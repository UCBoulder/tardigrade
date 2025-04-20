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

    private:
        // For varying material parameters
        const std::vector<std::string> & _user_material_prop_names;
        const std::vector<unsigned int> & _user_material_prop_indices;
        std::vector<const MaterialProperty<Real> *> _user_materials;
};

#endif