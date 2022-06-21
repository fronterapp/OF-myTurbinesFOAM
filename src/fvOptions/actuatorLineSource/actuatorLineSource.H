/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::actuatorLineSource

Description
    Actuator line class, which is a collection of actuator line elements.

SourceFiles
    actuatorLineSource.C

\*---------------------------------------------------------------------------*/

#ifndef actuatorLineSource_H
#define actuatorLineSource_H

#include "List.H"
#include "dictionary.H"
#include "vector.H"
#include "actuatorLineElement.H"
#include "cellSetOption.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class actuatorLineSource Declaration
\*---------------------------------------------------------------------------*/

class actuatorLineSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Total length of geometry in meters
        scalar totalLength_;

        //- Dictionary for the profiles' sectional coefficient data
        dictionary profileData_;

        //- Element geometry
        // point, spanDirection, chordLength, chordRefDir, pitch
        List<List<List<scalar> > > elementGeometry_;

        //- Element profiles
        List<word> elementProfiles_;

        //- Number of blade elements
        label nElements_;

        //- Free stream velocity
        vector freeStreamVelocity_;

        //- Free stream direction
        vector freeStreamDirection_;

        //- Total force vector from all elements
        vector force_;

        //- Force field from all elements
        volVectorField forceField_;

        //- List of actuator line elements
        PtrList<actuatorLineElement> elements_;

        //- Switch for writing performance
        bool writePerf_;

        //- Output file stream
        OFstream* outputFile_;

        //- Switch for harmonic pitching
        bool harmonicPitchingActive_;

        //- Reduced frequency of harmonic pitching = omega*c/(2*UInfty)
        scalar reducedFreq_;

        //- Amplitude of harmonic pitching in degrees
        scalar pitchAmplitude_;

        //- Time value to track whether to move
        scalar lastMotionTime_;

        //- Mean chord length of all elements
        scalar chordLength_;

        //- Aspect ratio
        scalar aspectRatio_;

        //- Switch for correcting end effects
        bool endEffectsActive_;

        //- Switch for multiphase simulations
        bool multiPhase_;

        //- Name of the phase which properties are to be used
        word phaseName_;

        //- Switch for harmonic floater motion
        bool harmonicFloaterActive_;

        //- Amplitude of floater translation motion (x y z, inertial frame) in meters 
        vector harmonicTraAmp_;

        //- Frequency of floater translation motion in Hz
        vector harmonicTraFreq_;

        //- Phase shift of floater translation motion in deg
        vector harmonicTraPhase_;

        //- Amplitude of floater rotation motion (roll pitch yaw, rotating frame) in deg
        vector harmonicRotAmp_;

        //- Frequency of floater rotation motion in Hz
        vector harmonicRotFreq_;

        //- Phase shift of floater translation motion in deg
        vector harmonicRotPhase_;

        //- Floater's center of rotation
        vector rotCenter_;

        //- Floater's center of rotation at previous time
        vector prevRotCenter_;

        //- Floater's orientation at t=0
        vector rot0_;

        //- Floater's orientation as a zy'x'' rotation Matrix
        tensor orientation_;

        //- Floater's orientation at previous time as a zy'x'' rotation Matrix
        tensor prevOrientation_;

    // Protected Member Functions

        //- Create actuator line elements
        void createElements();

        //- Read dictionary
        bool read(const dictionary& dict);

        //- Create the performance output file
        virtual void createOutputFile();

        //- Write performance to CSV
        void writePerf();

        //- Calculate end effects from lifting line theory
        void calcEndEffects();

        //- Execute harmonic pitching for a single time step
        void harmonicPitching();

        //- Read rigid body dictionary
        void readRigidBody(const fvMesh& mesh);

public:

    //- Runtime type information
    TypeName("actuatorLineSource");

    // Selectors

    //- Return a reference to the selected fvOption model
    static autoPtr<actuatorLineSource> New
    (
        const word& name,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Constructor
    actuatorLineSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~actuatorLineSource();


    // Member functions

        // Access

            //- Return const reference to the total force vector
            const vector& force();

            //- Return const reference to the force field
            const volVectorField& forceField();

            //- Return reference to element pointer list
            PtrList<actuatorLineElement>& elements();

            //- Return const reference to the current orientation
            const tensor& orientation();

            //- Return const reference to the previous orientation
            const tensor& prevOrientation();

            //- Return floater translation at current timestep
            const vector floaterTranslation();

            //- Return const reference to the current rotation center
            const vector& rotCenter();

            //- Return const reference to the previous rotation center
            const vector& prevRotCenter();

        // Edit

            //- Translate the actuator line
            void translate(vector translation);

            //- Rotate the actuator line about a specified axis
            void rotate(vector rotationPoint, vector axis, scalar radians);

            //- Rotate the actuator line using a rotation matrix
            void rotate(const vector &rotationPoint, const tensor &rotMatrix);

            //- Pitch the blade about its chord mount
            void pitch(scalar radians);

            //- Pitch the blade about some fraction of the chord
            void pitch(scalar radians, scalar chordFraction);

            //- Move the actuator line according to the floater motion. Account for translation velocity.
            void floaterMove(const vector &translation, const vector &velocity, const vector &omega);

            //- Return the intrinsic rotation matrix (zy'x'') given the roll-pitch-yaw angles
            tensor rotAngles2Matrix(const vector &rotAngles);

            //- Tranform a given axis-angle rotation into rotation matrix
            tensor rotAxis2Matrix(const vector &axis, const scalar &angle);

            //- Tranform a given rotation matrix into an axis-angle rotation
            void rotMatrix2Axis(const tensor &rotMatrix, vector &axis, scalar &angle);

            //- Scale the velocity of all elements
            void scaleVelocity(scalar scale);

            //- Set speed of all elements based on rotation
            void setSpeed(vector point, vector axis, scalar omega);

            //- Set floater velocity of all elements
            void setFloaterVelocity(const vector &velocity);

            //- Add velocity to all elements based on floater translation
            void addFloaterVelocity(const vector &velocity);

            //- Add velocity to all elements based on floater rotation
            void addFloaterOmega(const vector &point, const vector &axis, const scalar &omega);

            //- Set rotational velocity of all elements for flow curvature
            //  correction
            void setOmega(scalar omega);

             //- Execute harmonic floating motion for a single time step
            void harmonicFloaterMotion();

            //- Rotate floater according to initial orientation
            void floaterInitialise();
            
        // Evaluation

            //- Compute the moment about a given point
            vector moment(vector point);


        // IO

            //- Print dictionary values
            virtual void printCoeffs() const;


        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Source term to turbulence scalars
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldI
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
