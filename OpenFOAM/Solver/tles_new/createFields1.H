 
// ========================================================================
// Nadja

// Info<< "Reading field tau\n" << endl;
    volTensorField tau
    (
        IOobject
        (
            "tau",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("tau", dimensionSet(0,2,-2,0,0,0,0),tensor::zero)
    );

// ========================================================================
