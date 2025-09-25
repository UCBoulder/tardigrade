#include<CylindricalSurfaceDirichletBC.h>
#include<cmath>

//Register the object
registerMooseObject( "tardigradeApp", CylindricalSurfaceDirichletBC );

InputParameters
CylindricalSurfaceDirichletBC::validParams()
{
    InputParameters params = DirichletBCBase::validParams();
    params.addRequiredParam< Point >( "center", "A point on the cylindrical axis prior to motion" );
    params.addRequiredParam< Real >( "radius", "Radius of the cylindrical surface" );
    params.addRequiredParam< Real >( "velocity", "Magnitude of velocity of the plane along its normal" );
    params.addRequiredParam< RealVectorValue >( "axis", "Axis direction of the cylinder" );
    params.addRequiredParam< RealVectorValue >( "normal", "Direction of motion of surface" );
    params.addParam< bool >( "use_sector", false, "Limit BC to a sector between angle_min and angle_max" );
    params.addParam< Real >( "angle_min", 0.0, "Minimum angle in radians for active sector" );
    params.addParam< Real >( "angle_max", 6.28318530718, "Maximum angle in radians for active sector" );
    params.addParam< bool >( "invert_displacement", false, "Reverse the sign of the displacement correction" );
    params.addCoupledVar("displacements", "The displacements");
    params.addParam< Real >( "beta", 1e-2, "A damping coefficient to aid in stability" );

    return params;
}

CylindricalSurfaceDirichletBC::CylindricalSurfaceDirichletBC( const InputParameters & parameters )
    : // Call the constructor for the base class
        DirichletBCBase( parameters ),
        _u_dot( _var.uDot( ) ),
        _du_dot_du( _var.duDotDu( ) ),
        _center( getParam< Point >( "center" ) ),
        _radius( getParam< Real >( "radius" ) ),
        _velocity( getParam< Real >( "velocity" ) ),
        _axis( getParam< RealVectorValue >( "axis" ) ),
        _normal( getParam< RealVectorValue >( "normal" ) ),
        _use_sector( getParam< bool >( "use_sector" ) ),
        _angle_min( getParam< Real >( "angle_min" ) ),
        _angle_max( getParam< Real >( "angle_max" ) ),
        _invert_displacement( getParam< bool >( "invert_displacement" ) ),
        _beta( getParam< Real >( "beta" ) )
    {
        // Normalize direction of motion vector
        const Real len_normal = _normal.norm();
        if ( len_normal < 1e-8 ){
            mooseError( "Normal vector magnitude too small to normalize" );
        }
        _normal /= len_normal;

        // Normalize axis vector
        const Real len_axis = _axis.norm();
        if ( len_axis < 1e-8 ){
            mooseError( "Axis vector magnitude too small to normalize" );
        }
        _axis /= len_axis;

        // Determine the displacement direction based on the defined displacement variable
        // TODO: Figure out if we can combine displacement for multiple components
        const std::string & var_name = _var.name();
        if ( var_name.find( "x" ) != std::string::npos ){
            _disp_dir = RealVectorValue( 1.0, 0.0, 0.0 );
        }
        else if ( var_name.find( "y" ) != std::string::npos ){
            _disp_dir = RealVectorValue( 0.0, 1.0, 0.0 );
        }
        else if ( var_name.find( "z" ) != std::string::npos ){
            _disp_dir = RealVectorValue( 0.0, 0.0, 1.0 );
        }
        else{
            mooseError( "Cannot determine direction from variable name: ", var_name );
        }

        // Set the displacements
        auto num_disp = coupledComponents( "displacements" );
        for ( decltype( num_disp ) i = 0; i < num_disp; ++i ){
            _displacements.push_back( coupled( "displacements", i ) );
        }

}

Real CylindricalSurfaceDirichletBC::computeQpValue( ){

    mooseError( "This approach isn't allowed!" );

}

Real CylindricalSurfaceDirichletBC::computeQpResidual( ){

    if ( isOverclosed( ) ){

        Real d = computeSignedDistanceToSurface( *_current_node );

        const Point new_center = _center + _normal * ( _velocity * _t );
        RealVectorValue r_vec = *_current_node - new_center;
        Real axial_comp = r_vec * _axis;
        RealVectorValue radial_vec = r_vec - axial_comp * _axis;

        Real s = 1.0;
        if ( _invert_displacement ){
            s *= -1.0;
        }

//        std::cerr << _current_node->id( ) << ", " << *_current_node << ", " << s * displacement_vector << ", " << _disp_dir << ", " << displacement_vector * _disp_dir << "\n";
//        if ( _current_node->id( ) == 131 ){
//            std::cerr << "\n" << _var.number( ) << ", " << _u[ _qp ] << ", " << s * d * radial_vec / ( radial_vec.norm( ) + 1e-9 ) * _disp_dir << ", " << + _beta * _u_dot[ _qp ] << " : " << _u[ _qp ] + s * d * radial_vec / ( radial_vec.norm( ) + 1e-9 ) * _disp_dir + _beta * _u_dot[ _qp ] << "\n";
//        }

        return -( _u[ _qp ] + _beta * _u_dot[ _qp ] ) + s * d * radial_vec / ( radial_vec.norm( ) + 1e-9 ) * _disp_dir;

    }

    return 0;

}

Real CylindricalSurfaceDirichletBC::computeQpJacobian( ){

    return computeQpOffDiagJacobian( _var.number( ) );

}

Real CylindricalSurfaceDirichletBC::computeQpOffDiagJacobian( const unsigned int jvar_num ){

    if ( !isOverclosed( ) ){

        return 1e-9;

    }

    auto val = std::find( _displacements.begin(), _displacements.end(), jvar_num );

    Real J = 0;

    if ( jvar_num == _var.number( ) ){

        J -= 1 + _beta * _du_dot_du[ _qp ];

//        if ( _current_node->id( ) == 131 ){
//            std::cerr << "first jacobian for var.number: " << _var.number( ) << ", " << jvar_num << "\n";
//        }
    }

    if ( val != _displacements.end( ) ){

        const Point new_center = _center + _normal * ( _velocity * _t );
        RealVectorValue r_vec = *_current_node - new_center;
        Real axial_comp = r_vec * _axis;
        RealVectorValue radial_vec = r_vec - axial_comp * _axis;

        Real s = 1.0;
        if ( _invert_displacement ){
            s *= -1.0;
        }

        J += s * ( _radius / ( radial_vec.norm( ) * radial_vec.norm( ) * radial_vec.norm( ) + 1e-9 ) ) * radial_vec( val - _displacements.begin( ) ) * ( radial_vec * _disp_dir )
           + s * ( 1 - _radius / ( radial_vec.norm( ) + 1e-9 ) ) * ( _disp_dir - _axis * _disp_dir * _axis )( val - _displacements.begin( ) );

//        if ( _current_node->id( ) == 1 ){
//
//            std::cerr << "d: " << _current_node->id( ) << ", " << computeSignedDistanceToSurface( *_current_node ) << "\n";
//
//        }
//
//        if ( _current_node->id( ) == 131 ){
//            std::cerr << "var.number: " << _var.number( ) << ", " << jvar_num << ", " << s << ", " << *_current_node << ", " << _center << ", " << radial_vec << ", " << _normal << ", " << _disp_dir << ", " << J << "\n";
//        }

    }

    return J;

}

bool CylindricalSurfaceDirichletBC::shouldApply( ) const{ return isOverclosed( ); }

bool CylindricalSurfaceDirichletBC::isOverclosed() const{

    if ( !inAngularRange( *_current_node ) ){
        return false;
    }

    Real d = computeSignedDistanceToSurface( *_current_node );

    if ( _invert_displacement ){
        d *= -1.0;
    }

    return d <= 0;  // Apply only if node is outside or on the moving surface
}

bool CylindricalSurfaceDirichletBC::inAngularRange( const Point & pt ) const{

    if ( !_use_sector ){
        return true;
    }

    // Perform calculation using original center location
    RealVectorValue r_vec = pt - _center;
    Real axial_comp = r_vec * _axis;
    RealVectorValue radial_vec = r_vec - axial_comp * _axis;

    Real theta = std::atan2( radial_vec( 1 ), radial_vec( 0 ) );
    if ( theta < 0.0 ){
        theta += 6.28318530718;
    }

    return theta >= _angle_min && theta <= _angle_max;
}

Real CylindricalSurfaceDirichletBC::computeSignedDistanceToSurface( const Point & pt ) const{

    // Current location of center point
    const Point new_center = _center + _normal * ( _velocity * _t );

    // Vector from new center to point
    RealVectorValue r_vec = pt - new_center;
    Real axial_comp = r_vec * _axis;
    RealVectorValue radial_vec = r_vec - axial_comp * _axis;

    Real r_node = radial_vec.norm();
    return r_node - _radius;
}
