    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	// kinematic viscosity
    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    Info<< "Reading turbulenceProperties\n" << endl;

    IOdictionary turbulenceProperties
    (
        IOobject
        (  
	    "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Switch: for model choice
    dimensionedScalar modelSwitch
    (
        turbulenceProperties.lookup("modelSwitch")
    );

    // Filter Width Ratio
    dimensionedScalar filterWidthRatio
    (
        turbulenceProperties.lookup("filterWidthRatio")
    );

    // Pressure field
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Velocity field
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Deconvoluted velocity field
    Info<< "Reading field V\n" << endl;
    volVectorField V
    (
        IOobject
        (
            "V",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );

    // Instantaneous Kinetic Energy
    volScalarField KE
    ( 
        IOobject
        (
	    "KE",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE   
		),
		mesh,
		dimensionedScalar("KE", dimensionSet(0,2,-2,0,0,0,0),scalar(0.0))
    ); 

#   include "createPhiV.H" 
#   include "createPhi.H" 

    label pRefCell = 0;
    scalar pRefValue = 0.0;

    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
