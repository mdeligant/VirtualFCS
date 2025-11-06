within VirtualFCS.Electrical;

model SimpleBLDC
  import SI = Modelica.Units.SI;
  // --- Paramètres ---
  parameter SI.Resistance R = 3.536e-4 "Résistance phase équivalente";
  parameter SI.Inductance L = 1e-6 "Inductance équivalente (à ajuster)";
  parameter SI.Inertia    J = 1e-6 "Inertie du rotor (à ajuster)";
  parameter SI.RotationalDampingConstant B = 3.54e-6 "Frottement visqueux";
  parameter SI.Current I0 = 18.4 "Courant à vide";
  parameter Real Ke(unit="V.s/rad") = 0.001644071 "Constante f.e.m.";
  parameter Real Kt(unit="N.m/A")   = 0.00132626 "Constante de couple";

  // --- Ports ---


  // --- Etats ---
  SI.Current i(start=0);
  SI.AngularVelocity w(start=0);
  SI.Angle phi(start=0);

  // --- Auxiliaires ---
  SI.Voltage v;
  SI.Torque tau_em;
  Modelica.Electrical.Analog.Interfaces.Pin p annotation(
    Placement(transformation(origin = {50, 90}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {50, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin n annotation(
    Placement(transformation(origin = {-50, 88}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-50, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft annotation(
    Placement(transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}})));
  //Modelica.Mechanics.Rotational.Interfaces.Flange_b support annotation(    Placement(transformation(origin = {90, -88}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {90, -88}, extent = {{-10, -10}, {10, 10}})));
equation
  // Tension aux bornes
  v = p.v - n.v;

  // Équation électrique (moyennée)
  v = R*i + L*der(i) + Ke*w;

  // Couple électromagnétique (corrigé du courant à vide)
  tau_em = Kt*(i - I0);

  // Dynamique mécanique
  //shaft.tau = -(B*w) + tau_em;   // couple positif vers l'extérieur
  J*der(w) -(B*w) - tau_em = shaft.tau; // équilibre des couples
  der(phi) = w;

  // Courant qui sort de la borne positive
  p.i + n.i = 0;
  p.i = i;
  // Vitesse au port mécanique
  shaft.phi = phi;

annotation(
    Diagram(graphics),
    Icon(graphics = {Rectangle(fillColor = {0, 128, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-40, 60}, {80, -60}}), Rectangle(fillColor = {128, 128, 128}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-40, 60}, {-60, -60}}), Rectangle(fillColor = {95, 95, 95}, fillPattern = FillPattern.HorizontalCylinder, extent = {{80, 10}, {100, -10}}), Rectangle(lineColor = {95, 95, 95}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, extent = {{-40, 70}, {40, 50}}), Polygon(fillPattern = FillPattern.Solid, points = {{-50, -90}, {-40, -90}, {-10, -20}, {40, -20}, {70, -90}, {80, -90}, {80, -100}, {-50, -100}, {-50, -90}})}),
  experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end SimpleBLDC;