#ifndef CSVTimeHistoryDirichletBC_H
#define CSVTimeHistoryDirichletBC_H

#include<DirichletBCBase.h>
#include<vector>
#include<string>

class CSVTimeHistoryDirichletBC : public DirichletBCBase
{

    public:

        CSVTimeHistoryDirichletBC(const InputParameters & parameters);

        static InputParameters validParams();


    protected:

        virtual void initialSetup() override;

        virtual Real computeQpValue() override;

    private:

        struct Series
        {
            std::vector< Real > t;
            std::vector< Real > v;
            Real interpolate( Real time ) const;
        };

        // parameters
        const std::string _csv_path;
        const bool _enforce_zero_beyond;

        // constructed at setup
        std::unordered_map< dof_id_type, Series > _series_by_node;
        std::unordered_set< dof_id_type > _boundary_node_ids;

};

#endif