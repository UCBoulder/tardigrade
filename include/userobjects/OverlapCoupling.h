/*=============================================================================
 *                               OverlapCoupling                              *
 *============================================================================*
 * Communicate with the overlap coupling libraries outside of MOOSE to        *
 * perform the overlap coupling.                                              *
 *============================================================================*/

#ifndef OVERLAPCOUPLINGUSEROBJECT_H
#define OVERLAPCOUPLINGUSEROBJECT_H
#include "GeneralUserObject.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"

class OverlapCoupling;

class OverlapCoupling : public GeneralUserObject 
{

    public:
        OverlapCoupling( const InputParameters & parameters );
        
        static InputParameters validParams();

        virtual void initialize( ) override;

        virtual void execute( ) override;

        virtual void threadJoin( const UserObject & y ) override;

        virtual void finalize( ) override;

        void setAttribute( const std::unordered_map< unsigned int, unsigned int > &attribute, const std::string &attributeName );

        void setAttribute( const std::vector< double > &attribute, const std::string &attributeName );

        const std::unordered_map< unsigned int, unsigned int > *getMicroGlobalLocalNodeMap( ) const;
        const std::vector< double > *getMicroDisplacementDOF( ) const;
        const std::vector< double > *getMicroAugmentedLagrangianForce( ) const;

        const std::unordered_map< unsigned int, unsigned int > *getMacroGlobalLocalNodeMap( ) const;
        const std::vector< double > *getMacroDisplacementDOF( ) const;
        const std::vector< double > *getMacroAugmentedLagrangianForce( ) const;

    protected:

        const bool &_isMacro;
        const std::string &_overlapCouplingFilename;

        unsigned int _dim;

        std::unordered_map< unsigned int, unsigned int > _microGlobalLocalNodeMap;
        std::unordered_map< unsigned int, unsigned int > _macroGlobalLocalNodeMap;

        std::vector< double > _updatedMicroDisplacementDOF;
        std::vector< double > _updatedMacroDisplacementDOF;

        std::vector< double > _microAugmentedLagrangianForce;
        std::vector< double > _macroAugmentedLagrangianForce;

};

#endif
