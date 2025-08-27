#ifndef CYLINDRICALSURFACEDIRICHLETBC_H
#define CYLINDRICALSURFACEDIRICHLETBC_H

#include<DirichletBCBase.h>

class CylindricalSurfaceDirichletBC : public DirichletBCBase
{

    public:

        CylindricalSurfaceDirichletBC( const InputParameters &parameters );

        static InputParameters validParams();


    protected:

        virtual bool shouldApply() const override;

        virtual Real computeQpValue() override;

    private:

        Point _center; //!< A point on the cylindrical axis prior to motion
        Real _radius; //!< Radius of the cylindrical surface
        Real _velocity; //!< Magnitude of velocity of the plane along its normal
        RealVectorValue _axis; //!< Axis direction of the cylinder
        RealVectorValue _normal; //!< Direction of motion of surface
        RealVectorValue _disp_dir; //!< Direction associated with the Dirichlet BC degree of freedom

        bool _use_sector; //!< Limit BC to a sector between angle_min and angle_max
        Real _angle_min; //!< Minimum angle in radians for active sector
        Real _angle_max; //!< Maximum angle in radians for active sector
        bool _invert_displacement; //!< Reverse the sign of the displacement correction

        bool inAngularRange( const Point & pt ) const;
        Real computeSignedDistanceToSurface( const Point & pt ) const;
};

#endif