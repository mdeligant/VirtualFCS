within VirtualFCS.FuelCellStacksTesting;

model TestStack_with_resistor
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Modelica.Media.IdealGases.SingleGases.H2 constrainedby Modelica.Media.Interfaces.PartialSimpleIdealGasMedium;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack(redeclare package Cathode_Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {1, 3}, extent = {{-34, -34}, {34, 34}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureH2(nPorts = 1, p = 2e5, T = 293.15, redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-70, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceAir(nPorts = 1, m_flow = 0.25, T = 313.15, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {46, 20}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureAir(nPorts = 1, p = 1e5, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {50, -10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceWater(nPorts = 1, m_flow = 0.7, T = 323.15, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-26, -36}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkWater(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {24, -40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Fluid.Sources.MassFlowSource_T SourceH2(redeclare package Medium = Anode_Medium, T = 293.15, m_flow = 0.02, nPorts = 1) annotation(
    Placement(transformation(origin = {-70, 18}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.VariableResistor resistor annotation(
    Placement(transformation(origin = {2, 66}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Sources.Ramp ramp(height = -9.3, duration = 80, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {-76, 86}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(SourceAir.ports[1], fuelCellStack.port_a_Air) annotation(
    Line(points = {{36, 20}, {26, 20}, {26, 18}}, color = {0, 127, 255}));
  connect(SinkPressureAir.ports[1], fuelCellStack.port_b_Air) annotation(
    Line(points = {{40, -10}, {40, -12}, {26, -12}}, color = {0, 127, 255}));
  connect(SinkPressureH2.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-60, -14}, {-22, -14}, {-22, -12}, {-24, -12}}, color = {0, 127, 255}));
  connect(SourceH2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-60, 18}, {-24, 18}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, SourceWater.ports[1]) annotation(
    Line(points = {{-6, -18}, {-6, -36}, {-16, -36}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_Coolant, SinkWater.ports[1]) annotation(
    Line(points = {{8, -18}, {8, -40}, {14, -40}}, color = {0, 127, 255}));
  connect(fuelCellStack.pin_n, resistor.n) annotation(
    Line(points = {{-12, 38}, {-12, 66}, {-8, 66}}, color = {0, 0, 255}));
  connect(resistor.p, fuelCellStack.pin_p) annotation(
    Line(points = {{12, 66}, {14, 66}, {14, 38}}, color = {0, 0, 255}));
  connect(ramp.y, resistor.R) annotation(
    Line(points = {{-64, 86}, {2, 86}, {2, 78}}, color = {0, 0, 127}));
  annotation(
    experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.02));
end TestStack_with_resistor;