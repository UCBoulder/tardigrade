#include<MovingPlaneDirichletBC.h>
#include<cmath>

//Register the object
registerMooseObject( "tardigradeApp", MovingPlaneDirichletBC );

InputParameters
MovingPlaneDirichletBC::validParams()
{
    InputParameters params = DirichletBCBase::validParams();
    params.addRequiredParam< Point >("point", "A point on the plane prior to motion");
    params.addRequiredParam< Real >("velocity", "Magnitude of velocity of the plane along its normal");
    params.addRequiredParam< RealVectorValue >("normal", "Normal vector of the plane");
    return params;
}

MovingPlaneDirichletBC::MovingPlaneDirichletBC(const InputParameters & parameters)
    : // Call the constructor for the base class
        DirichletBCBase( parameters ),
        _p0( getParam< Point >( "point" ) ),
        _velocity( getParam< Real >( "velocity" ) ),
        _normal( getParam< RealVectorValue >( "normal" ) )
    {
        const Real len = _normal.norm();
        if (len < 1e-8)
            mooseError("Normal vector magnitude too small to normalize");
        _normal /= len;
}

Real MovingPlaneDirichletBC::computeQpValue( ){

    return signedDistanceToPlane( *_current_node );
}

bool MovingPlaneDirichletBC::shouldApply( ) const{

    const Real d = signedDistanceToPlane( *_current_node );

    // Only apply if the node is on the "normal side" of the plane (ahead of it)
    return d <= 1e-8;
}

Real MovingPlaneDirichletBC::signedDistanceToPlane(const Point & pt) const{

    const Point plane_position = _p0 + _normal * (_velocity * _t);
    return _normal * ( pt - plane_position );
}
