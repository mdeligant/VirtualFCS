within VirtualFCS.Fluid;

model CentrifugalCompressor

  Modelica.Units.SI.PerUnit PR "Pressure ratio";
  Modelica.Units.SI.PerUnit eta "Isentropic Efficiency";
  Modelica.Fluid.Interfaces.FluidPort_a inlet(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {-88, -72}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-90, -68}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Interfaces.FluidPort_b outlet(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {70, 48}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {70, 52}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Units.SI.Angle phi "position angulaire de l'arbre";
  Modelica.Units.SI.Torque tau "couple net sur l'arbre";
  Modelica.Units.SI.AngularVelocity omega "vitesse angulaire de l'arbre";
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component" annotation(
    choicesAllMatching = true);
  Medium.BaseProperties in_props;
  //(p(start = Pout_start), T(start = Tout_start), Xi(start = Xstart[1:Medium.nXi]));
  Medium.BaseProperties out_props;
  //(p(start = Pout_start), T(start = Tout_start), Xi(start = Xstart[1:Medium.nXi]));
  Modelica.Units.SI.MassFlowRate qm;
  Modelica.Units.SI.SpecificEnthalpy hin, hout, h_out_is;
  Modelica.Units.SI.Power Pcomp;
  Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation(
    Placement(transformation(origin = {-82, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b annotation(
    Placement(transformation(origin = {82, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}})));//
  parameter String file = Modelica.Utilities.Files.loadResource("modelica://VirtualFCS.Resources.CompressorData/TableCompresseurDimLess.txt") "Compressor data file" annotation(
    Dialog(group = "Table data definition", enable = true, loadSelector(filter = "Text files (*.txt);;MATLAB MAT-files (*.mat)", caption = "Open file in which table is present")));
  Modelica.Blocks.Tables.CombiTable2Ds table_PR(fileName = file, tableName = "tablePR", tableOnFile = true) annotation(
    Placement(transformation(origin = {-50, -30}, extent = {{-6, -8}, {6, 8}})));
  Modelica.Blocks.Tables.CombiTable2Ds table_eta(fileName = file, tableName = "tableEtaC", tableOnFile = true, extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint) annotation(
    Placement(transformation(origin = {-51, 11}, extent = {{-7, -9}, {7, 9}})));
  Modelica.Blocks.Tables.CombiTable2Ds table_qm(fileName = file, tableName = "tableQm", tableOnFile = true) annotation(
    Placement(transformation(origin = {-48, 51}, extent = {{-6, -7}, {6, 7}})));
  
  Real N "Ratio (%) of angular velocity over reference velocity";
  Real beta "Interpolation variable for data in table";

  parameter Modelica.Units.SI.Temperature T0 = 293.15;
  parameter Modelica.Units.SI.Pressure P0 = 1e5;
  parameter Modelica.Units.SI.AngularVelocity omega0 = 120000*Modelica.Constants.pi/30;
  Modelica.Units.SI.HeatCapacity r;
  Modelica.Units.SI.Temperature Tin;
  Modelica.Blocks.Interfaces.RealOutput sen_Air_comp_speed annotation(
    Placement(transformation(origin = {57, -5}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-15, -77}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));

parameter Real flow_scale_factor=1;

equation
sen_Air_comp_speed=omega;

  r=Modelica.Constants.R/Medium.molarMass(in_props.state);
  Tin=in_props.T;
  
// Mechanical boundary conditions
  flange_a.phi = phi;
  flange_b.phi = phi;
  flange_a.tau + flange_b.tau = tau;
  der(phi) = omega;
  inlet.m_flow = qm;
  outlet.p/inlet.p = PR;
  
  
  N = omega/omega0 / (Tin/T0)^0.5;
  table_eta.u1 = beta;
  table_eta.u2 = N;
  table_eta.y = eta;
  table_PR.u1 = beta;
  table_PR.u2 = N;
  max(table_PR.y,1) = PR;
  table_qm.u1 = beta;
  table_qm.u2 = N;
  table_qm.y*flow_scale_factor = qm*(T0/Tin)^0.5/(inlet.p/P0);
  
  
// Mass balance
  inlet.m_flow + outlet.m_flow = 0;
  
//Energy conservation
  qm*(hout - hin) = tau*omega;
  Pcomp = qm*(hout - hin);
  hin = inStream(inlet.h_outflow);
  h_out_is = Medium.isentropicEnthalpy(outlet.p, in_props.state);
  eta = (h_out_is - hin)/(hout - hin);
 
  inlet.h_outflow = inStream(outlet.h_outflow);
  outlet.h_outflow = hout;
  
  in_props.p = inlet.p;
  in_props.h = inStream(inlet.h_outflow);
  in_props.Xi = inStream(inlet.Xi_outflow);
  
  out_props.p = outlet.p; 
  out_props.h = hout;
  out_props.Xi = actualStream(outlet.Xi_outflow);
  
  
  inlet.Xi_outflow = inStream(outlet.Xi_outflow);
  inlet.C_outflow = inStream(outlet.C_outflow);

  outlet.Xi_outflow = inStream(inlet.Xi_outflow);
  outlet.C_outflow = inStream(inlet.C_outflow);

  
  
  annotation(
    Icon(graphics = {Polygon(points = {{-80, 80}, {80, 40}, {80, -40}, {-80, -80}, {-80, 80}}), Text(origin = {0, 2}, extent = {{-52, 24}, {52, -24}}, textString = "Compressor")}),
    Diagram(graphics = {Polygon(points = {{-80, 80}, {80, 40}, {80, -40}, {-80, -80}, {-80, 80}})}),
  experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));

end CentrifugalCompressor;