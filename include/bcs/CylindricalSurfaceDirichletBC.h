#ifndef CYLINDRICALSURFACEDIRICHLETBC_H
#define CYLINDRICALSURFACEDIRICHLETBC_H

#include<DirichletBCBase.h>

class CylindricalSurfaceDirichletBC : public DirichletBCBase
{

    public:

        CylindricalSurfaceDirichletBC(const InputParameters &parameters);

        static InputParameters validParams();


    protected:

        virtual bool shouldApply() const override;

        virtual Real computeQpValue() override;

    private:

        Point _center;
        Real _radius;
        Real _velocity;
        RealVectorValue _axis;
        RealVectorValue _normal;
        RealVectorValue _disp_dir;

        bool _use_sector;
        Real _angle_min;
        Real _angle_max;
        bool _invert_displacement;

        bool inAngularRange(const Point & pt) const;
        Real computeSignedDistanceToSurface(const Point & pt) const;
};

#endif