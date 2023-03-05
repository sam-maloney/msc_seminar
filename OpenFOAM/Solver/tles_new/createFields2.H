// ========================================================================
// TADM

// damping parameter
    dimensionedScalar chi
    (
        turbulenceProperties.lookup("chi")
    );

// Filter Width Ratio for regularization, must be >1
    dimensionedScalar filterWidthRatioReg
    (
        turbulenceProperties.lookup("filterWidthRatioReg")
    );

// Velocity field for regularization term
    volVectorField W
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		U
    );

