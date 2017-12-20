/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Application
    tles_new

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
	using Temporal Large Eddy Simulation (TLES) turbulence modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"  	// Initialize p U V nu filterWidthRatio modelSwitch 
    #include "createFields1.H"  // 1: Initialize tau for model1
    #include "createFields2.H"  // 2: Initialize W, chi, filterWidthRatioReg
	#include "createFields3.H"  // 3: Initialize U_ for divergence clearing  
 
    int model = modelSwitch.value();
    //Info<< "modelSwitch: " <<model<< endl;
    
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
 
    const dimensionedScalar timeStep ("timeStep", dimensionSet(0,0,1,0,0,0,0), double(runTime.deltaTValue())); // 1&2: deltaT
    const dimensionedScalar filterWidth    (filterWidthRatio*timeStep); // 1&2: filter width
    const dimensionedScalar filterWidthReg (filterWidthRatioReg*filterWidth); // 2: filter width for regularization
  
	while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNoV.H"       

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );
        
        switch (model){
		case 0 : // DNS
			break;

        case 1 : //TEDM without regularization
			{
				fvVectorMatrix model1
				(
					  fvm::ddt(U)
					+ fvm::div(phi, U)
					- fvm::laplacian(nu, U)
					== 
					- fvc::div(tau)
				); // 1: UEqn
			
				UEqn = model1;
            }
			break;
            
        case 2 : // TEDM with regularization
		case 3 : // TEDM with divergence cleaning
			fvVectorMatrix model2
			(
				  fvm::ddt(U)
				+ fvm::div(phi, U)
				- fvm::laplacian(nu, U)
				+ fvm::Sp(chi,U)
				== 
				  fvc::Sp(chi,W)
//				- fvc::Sp(chi,V)
				- fvc::div(tau)
			); // 2: UEqn
		 
			UEqn = model2;
            
			break;      
        }
        
		Info << "Created UEqn Model" << endl;
 
        solve(UEqn == -fvc::grad(p));

		Info << "Solved UEqn Model" << endl;

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();

            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p); 

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA) 
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

		Info << "Finished PISO Loop" << endl;

        switch (model){
		case 0 : //DNS
			V = U;
			break;

        case 1 : 
            #include "adjust1.H" // 1: Calculate tau and V
			break;

		case 3 :
			#include "adjust3.H" // 3: Clean divergence from U
			// intentional fall-through
        case 2 :
            #include "adjust1.H" // 1: Calculate tau and V
            #include "adjust2.H" // 2: Calculate W
			break;
        }

		Info << "Completed Adjustments" << endl;

        phiV = linearInterpolate(V) & mesh.Sf(); // 1&2: Update phiV

        dimensionedVector V_vol_ave 	// (Calculating the kinetic energy for evaluating the transient time)
        (
            "V_vol_ave",
            dimensionSet(0,1,-1,0,0,0,0),
            vector(gAverage(V))
		);       
        
		KE =0.5*magSqr(V-V_vol_ave);
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
