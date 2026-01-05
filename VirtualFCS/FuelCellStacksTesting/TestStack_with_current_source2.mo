within VirtualFCS.FuelCellStacksTesting;

model TestStack_with_current_source2
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Media.MoistHydrogenThreeComponents;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack(redeclare package Cathode_Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {1, 5}, extent = {{-34, -34}, {34, 34}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureH2(nPorts = 1, p = 2e5, T = 293.15, redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-70, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceAir(nPorts = 1, m_flow = 0.1, T = 341.15, redeclare package Medium = Cathode_Medium, X = {0.1904, 0.048, 0.7616}) annotation(
    Placement(transformation(origin = {48, 24}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureAir(nPorts = 1, p = 2e5, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {50, -10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceWater(nPorts = 1, m_flow = 2, T = 341.15, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-26, -46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkWater(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {24, -40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Fluid.Sources.MassFlowSource_T SourceH2(redeclare package Medium = Anode_Medium, T = 341.15, m_flow = 0.01, nPorts = 1, X = {0.69, 0.3, 0.01}) annotation(
    Placement(transformation(origin = {-66, 24}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.VariableResistor resistor annotation(
    Placement(transformation(origin = {4, 86}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(duration = 80, height = -9.3, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {-74, 106}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(SourceAir.ports[1], fuelCellStack.port_a_Air) annotation(
    Line(points = {{38, 24}, {38, 21}, {26, 21}}, color = {0, 127, 255}));
  connect(SinkPressureAir.ports[1], fuelCellStack.port_b_Air) annotation(
    Line(points = {{40, -10}, {40, -11}, {26, -11}}, color = {0, 127, 255}));
  connect(SinkPressureH2.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-60, -14}, {-22, -14}, {-22, -11}, {-24, -11}}, color = {0, 127, 255}));
  connect(SourceH2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-56, 24}, {-42, 24}, {-42, 21}, {-24, 21}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, SourceWater.ports[1]) annotation(
    Line(points = {{-6, -15}, {-6, -46}, {-16, -46}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_Coolant, SinkWater.ports[1]) annotation(
    Line(points = {{8, -15}, {8, -40}, {14, -40}}, color = {0, 127, 255}));
  connect(ramp.y, resistor.R) annotation(
    Line(points = {{-62, 106}, {4, 106}, {4, 98}}, color = {0, 0, 127}));
  connect(resistor.p, fuelCellStack.pin_p) annotation(
    Line(points = {{14, 86}, {14, 40}}, color = {0, 0, 255}));
  connect(resistor.n, fuelCellStack.pin_n) annotation(
    Line(points = {{-6, 86}, {-12, 86}, {-12, 40}}, color = {0, 0, 255}));
  annotation(
    experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.02));
end TestStack_with_current_source2;
