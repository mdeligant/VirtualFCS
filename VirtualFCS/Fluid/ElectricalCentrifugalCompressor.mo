within VirtualFCS.Fluid;

model ElectricalCentrifugalCompressor
  // System
  outer Modelica.Fluid.System system "System properties";
  //*** DEFINE REPLACEABLE PACKAGES ***//
  
  
  parameter String file = Modelica.Utilities.Files.loadResource("modelica://VirtualFCS.Resources.CompressorData/TableCompresseurDimLessExtended.txt") "Compressor data file" annotation(
    Dialog(group = "Table data definition", enable = true, loadSelector(filter = "Text files (*.txt);;MATLAB MAT-files (*.mat)", caption = "Open file in which table is present")));
    
  parameter Modelica.Units.SI.AngularVelocity omega0 = 120000*Modelica.Constants.pi/30;
  parameter Modelica.Units.SI.Voltage Vmax=600;

  parameter Real Ke(unit="V.s/rad") = 0.04 "Constante f.e.m.";
  parameter Real Kt(unit="N.m/A")   = 0.04 "Constante de couple";
  
  parameter Modelica.Units.SI.Temperature T0 = 293.15;
  parameter Modelica.Units.SI.Pressure P0 = 1e5;  
  parameter Real flow_scale_factor=1;

  // Medium declaration
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component" annotation(
    choicesAllMatching = true);
  
  //*** INSTANTIATE COMPONENTS ***//
  // Interfaces and boundaries
  Modelica.Fluid.Interfaces.FluidPort_a Input(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-78, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b Output(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {84, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(
    Placement(transformation(origin = {-24, 96}, extent = {{10, -10}, {-10, 10}}), iconTransformation(origin = {20, 60}, extent = {{-6, -6}, {6, 6}})));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(
    Placement(transformation(origin = {-48, 96}, extent = {{10, -10}, {-10, 10}}), iconTransformation(origin = {-16, 60}, extent = {{-6, -6}, {6, 6}})));
  // Control components
  Modelica.Blocks.Interfaces.RealOutput sen_Air_comp_speed annotation(
    Placement(transformation(origin = {98, -94}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {0, -66}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput controlInterface annotation(
    Placement(transformation(origin = {-99, 41}, extent = {{-13, -13}, {13, 13}}), iconTransformation(origin = {-113, 61}, extent = {{-13, -13}, {13, 13}})));
  // Machines
  Modelica.Mechanics.Rotational.Components.Inertia inertia(J = 1e-5) annotation(
    Placement(transformation(origin = {-6, -16}, extent = {{-10, -10}, {10, 10}})));
  // Sensors
  // Other
  // Power & Efficiencies
  Real Power_Compressor(unit = "W") "The power consumed by the Compressor";
  CentrifugalCompressor centrifugalCompressor(redeclare package Medium = Medium, file = file, T0 = T0, P0 = P0, omega0 = omega0, flow_scale_factor = flow_scale_factor) annotation(
    Placement(transformation(origin = {44, -70}, extent = {{-20, -20}, {20, 20}})));
  VirtualFCS.Electrical.SimpleBLDC simpleBLDC(J = 2e-4, Ke = Ke, Kt = Kt, w(start = 10.471975511965978, fixed = true, displayUnit = "rpm")) annotation(
    Placement(transformation(origin = {-35, -16}, extent = {{-10, -10}, {10, 10}})));
  VirtualFCS.Control.PID pid(CSmax = Vmax, CSmin = 0, Kp = 1, PVmax = omega0*1.3, PVmin = 0, Ti = 0.1) annotation(
    Placement(transformation(origin = {-65, 38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Machines.Examples.ControlledDCDrives.Utilities.IdealDcDc idealDcDc(Td = 1e-6, Ti = 1e-6) annotation(
    Placement(transformation(origin = {-34, 38}, extent = {{-10, -10}, {10, 10}})));
equation
  Power_Compressor = pin_p.i*pin_p.v;
  connect(Input, centrifugalCompressor.inlet) annotation(
    Line(points = {{-78, -60}, {-36, -60}, {-36, -84}, {26, -84}}));
  connect(inertia.flange_b, centrifugalCompressor.flange_a) annotation(
    Line(points = {{4, -16}, {20, -16}, {20, -70}, {26, -70}}));
  connect(simpleBLDC.shaft, inertia.flange_a) annotation(
    Line(points = {{-26, -16}, {-16, -16}}));
  connect(controlInterface, pid.SP) annotation(
    Line(points = {{-99, 41}, {-76, 41}}, color = {0, 0, 127}));
  connect(pin_n, idealDcDc.pin_nBat) annotation(
    Line(points = {{-48, 96}, {-48, 60}, {-44, 60}, {-44, 48}}, color = {0, 0, 255}));
  connect(pin_p, idealDcDc.pin_pBat) annotation(
    Line(points = {{-24, 96}, {-24, 48}}, color = {0, 0, 255}));
  connect(pid.CS, idealDcDc.vRef) annotation(
    Line(points = {{-55, 38}, {-47, 38}}, color = {0, 0, 127}));
  connect(idealDcDc.pin_nMot, simpleBLDC.n) annotation(
    Line(points = {{-44, 28}, {-44, 10}, {-40, 10}, {-40, -6}}, color = {0, 0, 255}));
  connect(idealDcDc.pin_pMot, simpleBLDC.p) annotation(
    Line(points = {{-24, 27.8}, {-24, 9.8}, {-30, 9.8}, {-30, -6.2}}, color = {0, 0, 255}));
  connect(centrifugalCompressor.sen_Air_comp_speed, sen_Air_comp_speed) annotation(
    Line(points = {{42, -86}, {42, -94}, {98, -94}}, color = {0, 0, 127}));
  connect(centrifugalCompressor.sen_Air_comp_speed, pid.PV) annotation(
    Line(points = {{42, -86}, {42, -94}, {-96, -94}, {-96, 34}, {-75, 34}}, color = {0, 0, 127}));
  connect(centrifugalCompressor.outlet, Output) annotation(
    Line(points = {{58, -60}, {84, -60}}, color = {0, 127, 255}));
  annotation(
    Icon(graphics = {Polygon(visible = false, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{20, -75}, {50, -85}, {20, -95}, {20, -75}}), Line(visible = false, points = {{55, -85}, {-60, -85}}, color = {0, 128, 255}), Polygon(visible = false, lineColor = {0, 128, 255}, fillColor = {0, 128, 255}, fillPattern = FillPattern.Solid, points = {{20, -70}, {60, -85}, {20, -100}, {20, -70}}), Polygon(origin = {15, 0}, fillColor = {208, 208, 208}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, points = {{-115, 80}, {-115, -80}, {85, -50}, {85, 50}, {85, 50}, {-115, 80}}), Line(origin = {5.2, 0.61}, points = {{-57, 0}, {47, 0}}, thickness = 1.5, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 20, smooth = Smooth.Bezier)}, coordinateSystem(initialScale = 0.1)),
    Documentation(info = "<html><head></head><body>The Compressor model is designed to compress air from the ambient environment to the desired pressure and maintain a sufficient mass flow rate to support the needs of the <a href=\"modelica:77VirtualFCS.Electrochemical.Hydrogen.FuelCellStack\">Fuel Cell Stack</a>.&nbsp;<div><br></div><div><b>Description</b></div><div><b><br></b></div><div>The model features 5 interfaces: fluid ports in and out, electrical ports for positive and negative pins, and a control port. The fluid ports provide the upstream and downstream connections for the compressor, the electrical ports connect to the low-voltage power supply for the BoP components and the control port sets the rpm of the compressor. A <a href=\"modelica://Modelica.Electrical.Machines.BasicMachines.DCMachines.DC_PermanentMagnet\">DC Motor</a> drives the compressor and is connected to an inertia and torque source, which is linked to the resistance of the <a href=\"modelica://Modelica.Fluid.Machines.PrescribedPump\">Pump</a>.</div></body></html>"));
end ElectricalCentrifugalCompressor;