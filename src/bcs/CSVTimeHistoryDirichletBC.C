#include "CSVTimeHistoryDirichletBC.h"

#include "Moose.h"
#include "libmesh/boundary_info.h"
#include "DelimitedFileReader.h"

#include<algorithm>

// Register the object
registerMooseObject( "tardigradeApp", CSVTimeHistoryDirichletBC );

InputParameters
CSVTimeHistoryDirichletBC::validParams()
{
    InputParameters params = DirichletBCBase::validParams();
    params.addRequiredParam< std::string >(
        "csv_file", "Path to CSV with time in first column, node data in subsequent columns, "
        "and header containing 'time, node_id_1, node_id_2, ..., node_id_n'." );
    params.addParam< bool >(
        "enforce_zero_beyond_range", true, "If true, prescribe 0.0 outside the CSV time range, "
        "otherwise hold end values." );
    params.addClassDescription(
        "Time varying Dirichlet BC specified for each node in a nodeset. Read data from a CSV "
        "using DelimitedFileReader into a time-varying series." );

    return params;
}

CSVTimeHistoryDirichletBC::CSVTimeHistoryDirichletBC( const InputParameters & parameters )
    : // Call the constructor for the base class
        DirichletBCBase( parameters ),
        _csv_path( getParam< std::string >( "csv_file" ) ),
        _enforce_zero_beyond( getParam< bool >( "enforce_zero_beyond_range" ))
    {

}

void CSVTimeHistoryDirichletBC::initialSetup(){
    DirichletBCBase::initialSetup();

    _boundary_node_ids.clear();
    const auto & mesh = _mesh.getMesh();

    // 1) Collect unique node IDs in the requested boundary set(s)
    const libMesh::BoundaryInfo & boundary_info = mesh.get_boundary_info();
    std::vector< dof_id_type > all_node_ids;
    std::vector< BoundaryID >  all_bids;
    boundary_info.build_node_list( all_node_ids, all_bids );

    // turn boundary IDs into a set
    std::unordered_set< BoundaryID > target_boundary_IDs( boundaryIDs().begin(), boundaryIDs().end() );

    // collect nodes that belong to the boundary
    std::set< dof_id_type > tmp_nodes;
    for ( std::size_t i = 0; i < all_node_ids.size(); ++i ){
        if ( target_boundary_IDs.count( all_bids[ i ] ) ){
            tmp_nodes.insert( all_node_ids[ i ] );
        }
    }

    _boundary_node_ids.insert( tmp_nodes.begin(), tmp_nodes.end() );

    // 2) Read CSV with MOOSE utility
    MooseUtils::DelimitedFileReader reader( _csv_path );
    reader.setHeaderFlag( MooseUtils::DelimitedFileReader::HeaderFlag::ON );
    reader.setFormatFlag( MooseUtils::DelimitedFileReader::FormatFlag::COLUMNS );
    reader.read();

    const auto & node_names = reader.getNames();
    const auto & node_data  = reader.getData();

    // Check headers and columns
    if ( node_data.size() != node_names.size() ){
        mooseError( "CSVTimeHistoryDirichletBC: Header/column count mismatch (", node_names.size(), " node_names vs ", node_data.size(), " columns)." );
    }

    // Extract the time data
    std::vector< Real > times = node_data[ 0 ];
    if ( times.empty() ){
        mooseError( "CSVTimeHistoryDirichletBC: Time data missing from '", _csv_path, "'." );
    }

    // Remaining columns are one per node (same length as times)
    const std::size_t n_node_cols = node_data.size() - 1;

    // 3) Parse node IDs from header and make sure they're in the boundary node set
    std::vector< dof_id_type > header_node_ids;
    header_node_ids.reserve( n_node_cols );
    for ( std::size_t j = 1; j < node_names.size(); ++j ){
        // Check that node id is provided
        if ( node_names[ j ].empty() ){
            mooseError( "CSVTimeHistoryDirichletBC: Blank header at column ", j, "." );
        }

        // Convert node id string to integer
        dof_id_type nid = 0;
        try { nid = static_cast< dof_id_type >( std::stoll( node_names[ j ] ) ); }
        catch ( ... ) { mooseError( "CSVTimeHistoryDirichletBC: Cannot parse node id from header '", node_names[ j ], "'." ); }

        if ( !_boundary_node_ids.count( nid ) ){
            std::string range_str = "NA..NA";
            if ( !_boundary_node_ids.empty() ){
                auto mm = std::minmax_element(_boundary_node_ids.begin(), _boundary_node_ids.end() );
                range_str = std::to_string( *mm.first ) + " to " + std::to_string( *mm.second );
            }
            mooseError( "CSVTimeHistoryDirichletBC: Node id ", nid, " (header) is not in the target set. Captured IDs range from ", range_str, "." );
        }
        header_node_ids.push_back( nid );
    }

    // 4) Copy node data to columns
    std::vector< std::vector< Real > > cols( n_node_cols );
    for ( std::size_t j = 0; j < n_node_cols; ++j ){
        cols[ j ] = node_data[ j + 1 ];
        // Make sure there are the same number of rows
        if ( cols[ j ].size() != times.size() ){
            mooseError( "CSVTimeHistoryDirichletBC: Column ", j + 1, " length (", cols[ j ].size(), ") != time column length (", times.size(), ")." );
        }
    }

    // 5) Build time series map for each node
    _series_by_node.clear();
    for ( std::size_t j = 0; j < n_node_cols; ++j ){
        Series s;
        s.t = times;
        s.v = std::move( cols[ j ] );
        _series_by_node.emplace( header_node_ids[ j ], std::move( s ) );
    }
}

Real CSVTimeHistoryDirichletBC::computeQpValue(){

    const dof_id_type nid = _current_node->id();

    auto node_i = _series_by_node.find( nid );
    if ( node_i == _series_by_node.end() ){
        mooseError( "CSVTimeHistoryDirichletBC: No CSV series for node id ", nid, " on the applied boundary." );
    }

    Real prescribed = node_i->second.interpolate( _t );

    if ( _enforce_zero_beyond ){
        const auto & tvec = node_i->second.t;
        if ( _t < tvec.front() || _t > tvec.back() ){
            prescribed = 0.0;
        }
    }

    return prescribed;
}

Real CSVTimeHistoryDirichletBC::Series::interpolate( Real time ) const {

    const auto n = t.size();
    if ( n == 0 ) return 0.0; // if no values provided in csv
    if ( n == 1 ) return v[ 0 ]; // if only one value provided, return constant for all times

    if ( time <= t.front() ) return v.front(); // if t is before first time, clamp to first value
    if ( time >= t.back() )  return v.back(); // if t is after final time, clamp to final value

    auto it = std::upper_bound( t.begin(), t.end(), time );
    size_t hi = std::distance( t.begin(), it );
    size_t lo = hi - 1;

    const Real t0 = t[ lo ], t1 = t[ hi ];
    const Real v0 = v[ lo ], v1 = v[ hi ];
    const Real delta = ( time - t0 ) / ( t1 - t0 );
    return v0 + delta * ( v1 - v0 );
}
