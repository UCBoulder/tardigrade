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

        struct Series //!< Struct for storing nodal Dirichlet BC time series unpacked from a CSV file
        {
            std::vector< Real > t; //!< The time values
            std::vector< Real > v; //!< The nodal Dirichlet BC values
            Real interpolate( Real time ) const; //!< Function for interpolating values between time values
        };

        // parameters
        const std::string _csv_path; //!< Path to CSV with time in first column, node data in subsequent columns, and header containing 'time, node_id_1, node_id_2, ..., node_id_n'
        const bool _enforce_zero_beyond; //!< If true, prescribe 0.0 outside the CSV time range, otherwise hold end values

        // constructed at setup
        std::unordered_map< dof_id_type, Series > _series_by_node; //!< Map between node IDs and their Dirichlet BC time series
        std::unordered_set< dof_id_type > _boundary_node_ids; //!< The set of node IDs contained in the requested nodeset

};

#endif