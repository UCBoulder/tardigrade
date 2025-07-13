#include<CylindricalSurfaceDirichletBC.h>
#include<cmath>

//Register the object
registerMooseObject( "tardigradeApp", CylindricalSurfaceDirichletBC );

InputParameters
CylindricalSurfaceDirichletBC::validParams()
{
    InputParameters params = DirichletBCBase::validParams();
    params.addRequiredParam< Point >("center", "A point on the cylindrical axis prior to motion");
    params.addRequiredParam< Real >("radius", "Radius of the cylindrical surface");
    params.addRequiredParam< Real >("velocity", "Magnitude of velocity of the plane along its normal");
    params.addRequiredParam< RealVectorValue >("axis", "Axis direction of the cylinder");
    params.addRequiredParam< RealVectorValue >("normal", "Direction of motion");
    params.addParam< bool >("use_sector", false, "Limit BC to a sector between angle_min and angle_max");
    params.addParam< Real >("angle_min", 0.0, "Minimum angle in radians for active sector");
    params.addParam< Real >("angle_max", 6.28318530718, "Maximum angle in radians for active sector");
    params.addParam< bool >("invert_displacement", false, "Reverse the sign of the displacement correction");

    return params;
}

CylindricalSurfaceDirichletBC::CylindricalSurfaceDirichletBC(const InputParameters & parameters)
    : // Call the constructor for the base class
        DirichletBCBase( parameters ),
        _center( getParam< Point >( "center" ) ),
        _radius( getParam< Real >( "radius" ) ),
        _velocity( getParam< Real >( "velocity" ) ),
        _axis( getParam< RealVectorValue >( "axis" ) ),
        _normal( getParam< RealVectorValue >( "normal" ) ),
        _use_sector( getParam< bool >( "use_sector" ) ),
        _angle_min( getParam< Real >( "angle_min" ) ),
        _angle_max( getParam< Real >( "angle_max" ) ),
        _invert_displacement( getParam< bool >( "invert_displacement" ) )
    {
        // Normalize direction of motion vector
        const Real len_normal = _normal.norm();
        if (len_normal < 1e-8)
            mooseError("Normal vector magnitude too small to normalize");
        _normal /= len_normal;

        // Normalize axis vector
        const Real len_axis = _axis.norm();
        if (len_axis < 1e-8)
            mooseError("Axis vector magnitude too small to normalize");
        _axis /= len_axis;

        // Determine the displacement direction based on the defined displacement variable
        // TODO: Figure out if we can combine displacement for multiple components
        const std::string & var_name = _var.name();
        if (var_name.find("x") != std::string::npos)
            _disp_dir = RealVectorValue(1.0, 0.0, 0.0);
        else if (var_name.find("y") != std::string::npos)
            _disp_dir = RealVectorValue(0.0, 1.0, 0.0);
        else if (var_name.find("z") != std::string::npos)
            _disp_dir = RealVectorValue(0.0, 0.0, 1.0);
        else
            mooseError("Cannot determine direction from variable name: ", var_name);

}

Real CylindricalSurfaceDirichletBC::computeQpValue(){

    Real d = computeSignedDistanceToSurface(*_current_node);

    RealVectorValue displacement_vector = -d * _normal;

    if (_invert_displacement)
        displacement_vector *= -1.0;

    return displacement_vector * _disp_dir;
}

bool CylindricalSurfaceDirichletBC::shouldApply() const{

    if (!inAngularRange(*_current_node))
        return false;

    Real d = computeSignedDistanceToSurface(*_current_node);

    if (_invert_displacement)
        d *= -1.0;

    return d <= 1e-8;  // Apply only if node is outside or on the moving surface
}

bool CylindricalSurfaceDirichletBC::inAngularRange(const Point & pt) const{

    if (!_use_sector)
        return true;

    // Perform calculation using original center location
    RealVectorValue r_vec = pt - _center;
    Real axial_comp = r_vec * _axis;
    RealVectorValue radial_vec = r_vec - axial_comp * _axis;

    Real theta = std::atan2(radial_vec(1), radial_vec(0));
    if (theta < 0.0)
        theta += 6.28318530718;

    return theta >= _angle_min && theta <= _angle_max;
}

Real CylindricalSurfaceDirichletBC::computeSignedDistanceToSurface(const Point & pt) const{

    // Current location of center point
    const Point new_center = _center + _normal * (_velocity * _t);

    // Vector from new center to point
    RealVectorValue r_vec = pt - new_center;
    Real axial_comp = r_vec * _axis;
    RealVectorValue radial_vec = r_vec - axial_comp * _axis;

    Real r_node = radial_vec.norm();
    return r_node - _radius;
}
