/*!
====================================================================
|              MicromorphicMaterialHeterogeneous.cpp               |
====================================================================
*/


#include "MicromorphicMaterial.h"
#include "MicromorphicMaterialHeterogeneous.h"

//MOOSE includes
#include "Function.h"
#include "micromorphic_material_library.h"

registerMooseObject("tardigradeApp", MicromorphicMaterialHeterogeneous);

InputParameters
MicromorphicMaterialHeterogeneous::validParams()
{
    InputParameters params = MicromorphicMaterial::validParams();

    // Add new parameters for unique material property names and indices
    params.addRequiredParam<std::vector<std::string>>(
        "user_material_prop_names", "Names of material properties obtained from file");
    params.addRequiredParam<std::vector<unsigned int>>(
        "user_material_prop_indices", "Index of material properties obtained from file");

    return params;
}

void MicromorphicMaterialHeterogeneous::initQpStatefulProperties(){
    // call base method
    MicromorphicMaterial::initQpStatefulProperties(); 
}

MicromorphicMaterialHeterogeneous::MicromorphicMaterialHeterogeneous(const InputParameters & parameters)
    : MicromorphicMaterial(parameters),

    // Add new parameters for unique material property names and indices
    _user_material_prop_names(getParam<std::vector<std::string>>("user_material_prop_names")),
    _user_material_prop_indices(getParam<std::vector<unsigned int>>("user_material_prop_indices"))

{
    /*!
    ===========================================
    |    MicromorphicMaterialHeterogeneous    |
    ===========================================

    The constructor for MicromorphicMaterialHeterogeneous.
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

void MicromorphicMaterialHeterogeneous::computeQpProperties(){
    // call base method
    MicromorphicMaterial::computeQpProperties(); 
}