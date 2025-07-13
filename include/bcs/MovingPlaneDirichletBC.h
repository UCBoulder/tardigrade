#ifndef MOVINGPLANEDIRICHLETBC_H
#define MOVINGPLANEDIRICHLETBC_H

#include<DirichletBCBase.h>

class MovingPlaneDirichletBC : public DirichletBCBase
{

    public:

        MovingPlaneDirichletBC(const InputParameters &parameters);

        static InputParameters validParams();


    protected:

        virtual bool shouldApply() const override;

        virtual Real computeQpValue() override;

    private:

        Point _p0;
        Real _velocity;
        RealVectorValue _normal;

        Real signedDistanceToPlane(const Point & pt) const;
};

#endif