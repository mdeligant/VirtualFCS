within VirtualFCS.FuelCellStacksTesting;

model TestStack_with_resistor
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Media.MoistHydrogenThreeComponents;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack(redeclare package Cathode_Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {-3, 9}, extent = {{-34, -34}, {34, 34}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureH2(nPorts = 1, p = 2e5, T = 293.15, redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-66, -6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceAir(nPorts = 1, m_flow = 0.25, T = 313.15, redeclare package Medium = Cathode_Medium, X = {0.78, 0.196, 0.024}) annotation(
    Placement(transformation(origin = {48, 24}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureAir(nPorts = 1, p = 2e5, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {50, -10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceWater(nPorts = 1, m_flow = 2, T = 341.15, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-26, -46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkWater(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {24, -40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Fluid.Sources.MassFlowSource_T SourceH2(redeclare package Medium = Anode_Medium, T = 293.15, m_flow = 0.005, nPorts = 1, X = {0.7, 0.29, 0.01}) annotation(
    Placement(transformation(origin = {-66, 24}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.VariableResistor resistor annotation(
    Placement(transformation(origin = {-6, 62}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(duration = 80, height = -9.95, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {-74, 74}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(SourceAir.ports[1], fuelCellStack.port_a_Air) annotation(
    Line(points = {{38, 24}, {38, 25}, {22, 25}}, color = {0, 127, 255}));
  connect(SinkPressureAir.ports[1], fuelCellStack.port_b_Air) annotation(
    Line(points = {{40, -10}, {40, -7}, {22, -7}}, color = {0, 127, 255}));
  connect(SinkPressureH2.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-56, -6}, {-22, -6}, {-22, -7}, {-28, -7}}, color = {0, 127, 255}));
  connect(SourceH2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-56, 24}, {-42, 24}, {-42, 25}, {-28, 25}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, SourceWater.ports[1]) annotation(
    Line(points = {{-10, -11}, {-10, -46}, {-16, -46}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_Coolant, SinkWater.ports[1]) annotation(
    Line(points = {{4, -11}, {4, -40}, {14, -40}}, color = {0, 127, 255}));
  connect(ramp.y, resistor.R) annotation(
    Line(points = {{-63, 74}, {-6, 74}}, color = {0, 0, 127}));
  connect(resistor.n, fuelCellStack.pin_n) annotation(
    Line(points = {{-16, 62}, {-16, 43}, {-17, 43}}, color = {0, 0, 255}));
  connect(resistor.p, fuelCellStack.pin_p) annotation(
    Line(points = {{4, 62}, {10, 62}, {10, 44}}, color = {0, 0, 255}));
  annotation(
    experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.02));
end TestStack_with_resistor;