within VirtualFCS.ComponentTesting;

model TestStack  

  // Medium decleration
  replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Modelica.Media.IdealGases.SingleGases.H2 constrainedby Modelica.Media.Interfaces.PartialSimpleIdealGasMedium;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  
  Electrochemical.Hydrogen.FuelCellStackWithMassFlowCorrectionPressureDropHeatTransfer fuelCellStack annotation(
    Placement(transformation(origin = {2, -16}, extent = {{-43, -64.5}, {43, 64.5}})));
  Modelica.Fluid.Sources.Boundary_pT H2_pressure_out(nPorts = 1, redeclare package Medium = Anode_Medium, p = 2e5)  annotation(
    Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_H2(nPorts = 1, m_flow = 0.003, redeclare package Medium = Anode_Medium)  annotation(
    Placement(transformation(origin = {-68, 26}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_Air(nPorts = 1, redeclare package Medium = Cathode_Medium, m_flow = 0.25, T = 313.15)  annotation(
    Placement(transformation(origin = {64, 26}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Fluid.Sources.Boundary_pT Air_pressure_out(nPorts = 1, p = 2e5, redeclare package Medium = Cathode_Medium)  annotation(
    Placement(transformation(origin = {68, -12}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_water(nPorts = 1, m_flow = 1.5, redeclare package Medium = Coolant_Medium, T = 323.15)  annotation(
    Placement(transformation(origin = {-40, -46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT Water_pressure_out(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium)  annotation(
    Placement(transformation(origin = {10, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Electrical.Analog.Basic.Resistor resistor(R = 0.7)  annotation(
    Placement(transformation(origin = {0, 72}, extent = {{-10, -10}, {10, 10}})));

Real mass_balance;
equation
  mass_balance = fuelCellStack.H2_mflow.y - fuelCellStack.O2_mflow.y + fuelCellStack.H2O_mflow.y;
  connect(Source_massflow_H2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-58, 26}, {-30, 26}}, color = {0, 127, 255}));
  connect(H2_pressure_out.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-64, -14}, {-30, -14}}, color = {0, 127, 255}));
  connect(Source_massflow_Air.ports[1], fuelCellStack.port_a_Air) annotation(
    Line(points = {{54, 26}, {34, 26}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_b_Air, Air_pressure_out.ports[1]) annotation(
    Line(points = {{34, -14}, {58, -14}, {58, -12}}, color = {0, 127, 255}));
  connect(Water_pressure_out.ports[1], fuelCellStack.port_b_Coolant) annotation(
    Line(points = {{10, -60}, {10, -20}}, color = {0, 127, 255}));
  connect(resistor.n, fuelCellStack.pin_p) annotation(
    Line(points = {{10, 72}, {20, 72}, {20, 48}}, color = {0, 0, 255}));
  connect(resistor.p, fuelCellStack.pin_n) annotation(
    Line(points = {{-10, 72}, {-16, 72}, {-16, 48}}, color = {0, 0, 255}));
  connect(fuelCellStack.port_a_Coolant, Source_massflow_water.ports[1]) annotation(
    Line(points = {{-6, -20}, {-6, -46}, {-30, -46}}, color = {0, 127, 255}));

annotation(
    experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end TestStack;