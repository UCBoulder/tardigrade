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

        Point _p0; //!< A point on the plane prior to motion
        Real _velocity; //!< Magnitude of velocity of the plane along its normal
        RealVectorValue _normal; //!< Normal vector of the plane

        Real signedDistanceToPlane(const Point & pt) const;
};

#endif