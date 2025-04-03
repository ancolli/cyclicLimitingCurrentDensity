singlePhaseTransportModel laminarTransport(U, phi);

// Parameters Common to both solvers
// max concentration [K]
dimensionedScalar cmax("cmax", dimless, laminarTransport); //dimMoles/dimVolume

// max temperature [K]
dimensionedScalar Tmax("Tmax", dimTemperature, laminarTransport);

// Laminar Sc number
dimensionedScalar Sc("Sc", dimless, laminarTransport);

// Turbulent Sc number
dimensionedScalar Sct("Sct", dimless, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// electrical conductivity
dimensionedScalar kf("kf", dimensionSet (-1,-3,3,0,0,2,0), laminarTransport);

// Domain length [m]
dimensionedScalar Lc("Lc", dimLength, laminarTransport);

// Ubar needed for flow direction [m/s]
dimensionedVector vel("flowVel", dimVelocity, laminarTransport);

// Initialize the flowDir vector as it doesn't change
const vector flowDir = (vel.value()/mag(vel.value()));

// desired delay before calculating lambda
dimensionedScalar delay("delay", dimless, laminarTransport);

// interval between calculations of lambda
//dimensionedScalar interval("interval", dimless, laminarTransport);

// Initiate Nu
//dimensionedScalar Sh("Sh", dimless, 0.0);

// Turbulent lambda number [1/m]
dimensionedScalar lambda("lambda", dimless/dimLength, 0.0);

// switch between different methods
const word& method = word(laminarTransport.lookup("method"));

// Thermal expansion coefficient [1/K]
//dimensionedScalar gamma("gamma", dimTemperature/dimLength, laminarTransport);

// Turbulent lambda number [1/m]
dimensionedScalar lambda0("lambda0", dimless/dimLength, laminarTransport);

// Reference temperature [K]
//dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);
