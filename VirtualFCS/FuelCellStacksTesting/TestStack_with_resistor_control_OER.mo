within VirtualFCS.FuelCellStacksTesting;

model TestStack_with_resistor_control_OER
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Media.MoistHydrogenThreeComponents;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack(redeclare package Cathode_Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {1, 5}, extent = {{-34, -34}, {34, 34}})));
  Modelica.Fluid.Sources.Boundary_pT SourcePressureAir(nPorts = 1, p = 2e5, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {52, 20}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceWater(nPorts = 1, m_flow = 2, T = 341.15, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-26, -46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkWater(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {24, -40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Fluid.Valves.ValveLinear valveLinearAir(redeclare package Medium = Cathode_Medium, dp_nominal = 5000, m_flow_nominal = 0.2, dp(start = 1e4)) annotation(
    Placement(transformation(origin = {58, -10}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(redeclare package Medium = Cathode_Medium, nPorts = 1, p = 1e5) annotation(
    Placement(transformation(origin = {90, -10}, extent = {{10, -10}, {-10, 10}})));
  VirtualFCS.Control.PID pid_OER_air(CSmax = 1, CSmin = 0, Kp = 1, PVmax = 4, PVmin = 0.1, Ti = 0.1, CSs(start = 1, fixed = true)) annotation(
    Placement(transformation(origin = {76, -48}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp air_OER_ramp(duration = 5, height = -348, offset = 350, startTime = 10) annotation(
    Placement(transformation(origin = {123, -44}, extent = {{10, -10}, {-10, 10}})));
  inner Modelica.Fluid.System system(m_flow_start = 1e-4)  annotation(
    Placement(transformation(origin = {90, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureH2(redeclare package Medium = Anode_Medium, T = 293.15, nPorts = 1, p = 1.7e5) annotation(
    Placement(transformation(origin = {-82, -18}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceH2(redeclare package Medium = Anode_Medium, T = 293.15, X = {0.7, 0.29, 0.01}, m_flow = 0.004042092, nPorts = 1) annotation(
    Placement(transformation(origin = {-78, 18}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.RampCurrent rampCurrent(I = 440, duration = 80, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {2, 64}, extent = {{10, -10}, {-10, 10}})));
equation
  connect(fuelCellStack.port_a_Coolant, SourceWater.ports[1]) annotation(
    Line(points = {{-6, -15}, {-6, -46}, {-16, -46}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_Coolant, SinkWater.ports[1]) annotation(
    Line(points = {{8, -15}, {8, -40}, {14, -40}}, color = {0, 127, 255}));
  connect(valveLinearAir.opening, pid_OER_air.CS) annotation(
    Line(points = {{58, -18}, {58, -48}, {66, -48}}, color = {0, 0, 127}));
  connect(air_OER_ramp.y, pid_OER_air.SP) annotation(
    Line(points = {{112, -44}, {86, -44}}, color = {0, 0, 127}));
  connect(SourcePressureAir.ports[1], fuelCellStack.port_a_Air) annotation(
    Line(points = {{42, 20}, {26, 20}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_Air, valveLinearAir.port_a) annotation(
    Line(points = {{26, -10}, {48, -10}}, color = {0, 127, 255}));
  connect(valveLinearAir.port_b, AirSink.ports[1]) annotation(
    Line(points = {{68, -10}, {80, -10}}, color = {0, 127, 255}));
  connect(fuelCellStack.OER_Cathode, pid_OER_air.PV) annotation(
    Line(points = {{26, -16}, {42, -16}, {42, -92}, {104, -92}, {104, -52}, {86, -52}}, color = {0, 0, 127}));
  connect(SourceH2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-68, 18}, {-24, 18}, {-24, 20}}, color = {0, 127, 255}));
  connect(SinkPressureH2.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-72, -18}, {-52, -18}, {-52, -10}, {-24, -10}}, color = {0, 127, 255}));
  connect(rampCurrent.p, fuelCellStack.pin_p) annotation(
    Line(points = {{12, 64}, {14, 64}, {14, 40}}, color = {0, 0, 255}));
  connect(rampCurrent.n, fuelCellStack.pin_n) annotation(
    Line(points = {{-8, 64}, {-12, 64}, {-12, 40}}, color = {0, 0, 255}));
  annotation(
    experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.02));
end TestStack_with_resistor_control_OER;