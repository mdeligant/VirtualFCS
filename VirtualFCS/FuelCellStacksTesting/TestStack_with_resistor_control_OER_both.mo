within VirtualFCS.FuelCellStacksTesting;

model TestStack_with_resistor_control_OER_both
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
  Control.PID pid_OER_air(CSmax = 1, CSmin = 0.001, Kp = 1, PVmax = 100, PVmin = 0, Ti = 0.1, CSs(start = 0.01, fixed = true)) annotation(
    Placement(transformation(origin = {76, -48}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp air_OER_ramp(duration = 5, height = 1.5, offset = 2, startTime = 10) annotation(
    Placement(transformation(origin = {123, -44}, extent = {{10, -10}, {-10, 10}})));
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {86, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureH2(redeclare package Medium = Anode_Medium, T = 293.15, nPorts = 1, p = 1e5) annotation(
    Placement(transformation(origin = {-90, -12}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SourcePressureH2(redeclare package Medium = Anode_Medium, T = 293.15, nPorts = 1, p = 2e5) annotation(
    Placement(transformation(origin = {-84, 22}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Valves.ValveLinear valveLinearH2(redeclare package Medium = Anode_Medium, dp(start = 1e4), dp_nominal = 1e4, m_flow_nominal = 0.002) annotation(
    Placement(transformation(origin = {-52, -12}, extent = {{10, 10}, {-10, -10}})));
  VirtualFCS.Control.PID pid_OER_H2(CSmax = 1, CSmin = 0.001, CSs(start = 0.01, fixed = true), Kp = 1, PVmax = 100, PVmin = 0, Ti = 0.1) annotation(
    Placement(transformation(origin = {-82, -42}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp H2_OER_ramp(duration = 5, height = 2, offset = 1.1, startTime = 1) annotation(
    Placement(transformation(origin = {-139, -38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.VariableResistor resistor annotation(
    Placement(transformation(origin = {-2, 76}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(duration = 80, height = -9.3, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {-70, 88}, extent = {{-10, -10}, {10, 10}})));
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
  connect(SourcePressureH2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-74, 22}, {-50, 22}, {-50, 20}, {-24, 20}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_H2, valveLinearH2.port_a) annotation(
    Line(points = {{-24, -10}, {-42, -10}, {-42, -12}}, color = {0, 127, 255}));
  connect(valveLinearH2.port_b, SinkPressureH2.ports[1]) annotation(
    Line(points = {{-62, -12}, {-80, -12}}, color = {0, 127, 255}));
  connect(pid_OER_H2.CS, valveLinearH2.opening) annotation(
    Line(points = {{-72, -42}, {-52, -42}, {-52, -20}}, color = {0, 0, 127}));
  connect(H2_OER_ramp.y, pid_OER_H2.SP) annotation(
    Line(points = {{-128, -38}, {-92, -38}}, color = {0, 0, 127}));
  connect(fuelCellStack.OER_Anode, pid_OER_H2.PV) annotation(
    Line(points = {{-24, -16}, {-40, -16}, {-40, -78}, {-108, -78}, {-108, -46}, {-92, -46}}, color = {0, 0, 127}));
  connect(ramp.y, resistor.R) annotation(
    Line(points = {{-59, 88}, {-2, 88}}, color = {0, 0, 127}));
  connect(resistor.p, fuelCellStack.pin_p) annotation(
    Line(points = {{8, 76}, {14, 76}, {14, 40}}, color = {0, 0, 255}));
  connect(resistor.n, fuelCellStack.pin_n) annotation(
    Line(points = {{-12, 76}, {-12, 40}}, color = {0, 0, 255}));
  annotation(
    experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-06, Interval = 0.02),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts",
  __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end TestStack_with_resistor_control_OER_both;
