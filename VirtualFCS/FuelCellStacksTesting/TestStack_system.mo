within VirtualFCS.FuelCellStacksTesting;

model TestStack_system
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Media.MoistHydrogenThreeComponents;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack(redeclare package Cathode_Medium = Cathode_Medium, N_el = 10, m_pt_an = 0.1/100, m_pt_ca = 0.1/100) annotation(
    Placement(transformation(origin = {8, 19}, extent = {{-60, -60}, {60, 60}})));
  Modelica.Fluid.Sources.MassFlowSource_T SourceWater(nPorts = 1, m_flow = 0.5, T = 343.15, redeclare package Medium = Coolant_Medium, use_m_flow_in = true, use_T_in = true) annotation(
    Placement(transformation(origin = {-16, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Fluid.Sources.Boundary_pT SinkWater(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {30, -83}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
  Modelica.Fluid.Valves.ValveLinear valve_pressure_air(redeclare package Medium = Cathode_Medium, dp_nominal = 1e4, m_flow_nominal = 0.2, dp(start = 1e4)) annotation(
    Placement(transformation(origin = {69, -31}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Modelica.Fluid.Sources.Boundary_pT AirSink(redeclare package Medium = Cathode_Medium, nPorts = 1, p = 1e5) annotation(
    Placement(transformation(origin = {128, -62}, extent = {{10, -10}, {-10, 10}})));
  Control.PID pid_presure_cat(CSmax = 1, CSmin = 0.001, Kp = -1, PVmax = 5e5, PVmin = 0.1e5, Ti = 0.1, CSs(start = 0.1, fixed = true)) annotation(
    Placement(transformation(origin = {100, -31}, extent = {{-7, -7}, {7, 7}}, rotation = 180)));
  Modelica.Blocks.Sources.Ramp OER_set(duration = 80, height = 0, offset = 1.43, startTime = 1000) annotation(
    Placement(transformation(origin = {166, 39}, extent = {{10, -10}, {-10, 10}})));
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {138, 84}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT SinkPressureH2(redeclare package Medium = Anode_Medium, T = 293.15, nPorts = 1, p = 1e5) annotation(
    Placement(transformation(origin = {-136.5, -8.5}, extent = {{-10.5, -10.5}, {10.5, 10.5}})));
  Modelica.Fluid.Sources.Boundary_pT SourcePressureH2(redeclare package Medium = Anode_Medium, T = 293.15, nPorts = 1, p = 7e7, X = {1, 0, 0}, use_p_in = false) annotation(
    Placement(transformation(origin = {-134, 48}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Valves.ValveLinear purge_Valve(redeclare package Medium = Anode_Medium, dp(start = 1e4), dp_nominal = 1e4, m_flow_nominal = 0.002) annotation(
    Placement(transformation(origin = {-75, -9}, extent = {{6, 6}, {-6, -6}})));
  Control.PID pid_purge_valve(CSmax = 1, CSmin = 0, CSs(start = 0.01, fixed = true), Kp = -0.1, PVmax = 2, PVmin = 0, Ti = 1) annotation(
    Placement(transformation(origin = {-75, -32}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  Modelica.Blocks.Sources.Ramp purge_valve_set(duration = 80, height = 0, offset = 0.9, startTime = 1000) annotation(
    Placement(transformation(origin = {-170, -63}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSource(redeclare package Medium = Cathode_Medium, nPorts = 1, X = {0.2, 0.02, 0.78}) annotation(
    Placement(transformation(origin = {112, 85}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Fluid.ElectricalCentrifugalCompressor air_Compressor(redeclare package Medium = Cathode_Medium, flow_scale_factor = 1) annotation(
    Placement(transformation(origin = {75, 48}, extent = {{11, 11}, {-11, -11}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {86, 12}, extent = {{-6, -6}, {6, 6}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600) annotation(
    Placement(transformation(origin = {75, 22}, extent = {{-7, -7}, {7, 7}})));
  Modelica.Blocks.Sources.Ramp pressure_Air_set(duration = 600, height = 0.2e5, offset = 1.1e5, startTime = 300) annotation(
    Placement(transformation(origin = {165, -34}, extent = {{9, -9}, {-9, 9}})));
  Modelica.Fluid.Sensors.Pressure pressure_Air(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {117, -9}, extent = {{-7, -7}, {7, 7}}, rotation = -90)));
  Control.PID pid_OER(CSmax = 150000*Modelica.Constants.pi/30, CSmin = 15000*Modelica.Constants.pi/30, Kp = 0.1, PVmax = 10, PVmin = 0, Ti = 0.1) annotation(
    Placement(transformation(origin = {122, 41}, extent = {{5, 5}, {-5, -5}})));
  Modelica.Electrical.Analog.Sources.RampVoltage rampVoltage(V = -90, duration = 800, offset = 340, startTime = 100) annotation(
    Placement(transformation(origin = {8, 79}, extent = {{10, -10}, {-10, 10}})));
  Fluid.ElectricalCentrifugalCompressor recirculating_pump(redeclare package Medium = Anode_Medium, flow_scale_factor = 0.01) annotation(
    Placement(transformation(origin = {-54, 20}, extent = {{-8, -8}, {8, 8}}, rotation = 90)));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = 600) annotation(
    Placement(transformation(origin = {-74, 24}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(transformation(origin = {-89, 19}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  Modelica.Blocks.Sources.Ramp pump_speed(duration = 800, height = 100000*Modelica.Constants.pi/30, offset = 20000*Modelica.Constants.pi/30, startTime = 100) annotation(
    Placement(transformation(origin = {-169.5, 10.5}, extent = {{-9.5, -9.5}, {9.5, 9.5}})));
  Modelica.Blocks.Sources.Ramp Pressure_H2_set(duration = 80, height = 0.4e5, offset = 1.1e5, startTime = 10) annotation(
    Placement(transformation(origin = {-169, 89}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal(redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-54, -9}, extent = {{-6, -6}, {6, 6}})));
  Modelica.Fluid.Fittings.TeeJunctionIdeal teeJunctionIdeal1(redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-54, 47}, extent = {{-6, -6}, {6, 6}}, rotation = 180)));
  Modelica.Blocks.Sources.Ramp Temperature_coolant(duration = 800, offset = 70 + 273, startTime = 100) annotation(
    Placement(transformation(origin = {-2, -119}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Fluid.Valves.ValveLinear valve_pressure_an(redeclare package Medium = Cathode_Medium, dp(start = 6.989e7, fixed = true), dp_nominal = 7e7, m_flow_nominal = 0.004) annotation(
    Placement(transformation(origin = {-100, 48}, extent = {{-6, -6}, {6, 6}})));
  Modelica.Fluid.Sensors.Pressure pressure_H2(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {-70, 88}, extent = {{7, -7}, {-7, 7}})));
  Control.PID pid_presure_H2(CSmax = 1, CSmin = 0.001, CSs(start = 0.1, fixed = true), Kp = 5, PVmax = 5e5, PVmin = 0.1e5, Ti = 0.1) annotation(
    Placement(transformation(origin = {-100, 69}, extent = {{7, -7}, {-7, 7}}, rotation = 90)));
  Modelica.Blocks.Sources.Ramp m_flow_coolant(duration = 800, height = 4, offset = 0.3, startTime = 100) annotation(
    Placement(transformation(origin = {-54, -119}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Real efficiency;
equation
  efficiency = if time > 100 then fuelCellStack.eta_FC_LHV*(fuelCellStack.H2_mass_flow + fuelCellStack.H2_mem_crossover)/(valve_pressure_an.m_flow*inStream(valve_pressure_an.port_a.Xi_outflow[1])) else 0;
  connect(purge_Valve.port_b, SinkPressureH2.ports[1]) annotation(
    Line(points = {{-81, -9}, {-105.5, -9}, {-105.5, -8.5}, {-126, -8.5}}, color = {0, 127, 255}));
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{86, 18}, {86, 22}, {82, 22}}, color = {0, 0, 255}));
  connect(rampVoltage.n, fuelCellStack.pin_n) annotation(
    Line(points = {{-2, 79}, {-16, 79}}, color = {0, 0, 255}));
  connect(pid_OER.SP, OER_set.y) annotation(
    Line(points = {{127, 39}, {155, 39}}, color = {0, 0, 127}));
  connect(teeJunctionIdeal.port_3, recirculating_pump.Input) annotation(
    Line(points = {{-54, -3}, {-54, 11}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal.port_2, fuelCellStack.port_b_H2) annotation(
    Line(points = {{-48, -9}, {-36, -9}}, color = {0, 127, 255}));
  connect(teeJunctionIdeal1.port_1, fuelCellStack.port_a_H2) annotation(
    Line(points = {{-48, 47}, {-36, 47}}, color = {0, 127, 255}));
  connect(recirculating_pump.Output, teeJunctionIdeal1.port_3) annotation(
    Line(points = {{-54, 29}, {-54, 40.8}}, color = {0, 127, 255}));
  connect(pid_purge_valve.CS, purge_Valve.opening) annotation(
    Line(points = {{-75, -26}, {-75, -14}}, color = {0, 0, 127}));
  connect(pressure_Air_set.y, pid_presure_cat.SP) annotation(
    Line(points = {{155, -34}, {107, -34}}, color = {0, 0, 127}));
  connect(AirSource.ports[1], air_Compressor.Input) annotation(
    Line(points = {{112, 75}, {112.5, 75}, {112.5, 48}, {87, 48}}, color = {0, 127, 255}));
  connect(pid_OER.CS, air_Compressor.controlInterface) annotation(
    Line(points = {{117, 41}, {87, 41}}, color = {0, 0, 127}));
  connect(purge_Valve.port_a, teeJunctionIdeal.port_1) annotation(
    Line(points = {{-69, -9}, {-60, -9}}, color = {0, 127, 255}));
  connect(constantVoltage1.p, recirculating_pump.pin_p) annotation(
    Line(points = {{-74, 29}, {-64, 29}, {-64, 22}, {-59, 22}}, color = {0, 0, 255}));
  connect(constantVoltage1.n, recirculating_pump.pin_n) annotation(
    Line(points = {{-74, 19}, {-59, 19}}, color = {0, 0, 255}));
  connect(ground1.p, constantVoltage1.n) annotation(
    Line(points = {{-83, 19}, {-74, 19}}, color = {0, 0, 255}));
  connect(SourcePressureH2.ports[1], valve_pressure_an.port_a) annotation(
    Line(points = {{-124, 48}, {-106, 48}}, color = {0, 127, 255}));
  connect(valve_pressure_an.port_b, teeJunctionIdeal1.port_2) annotation(
    Line(points = {{-94, 48}, {-60, 48}}, color = {0, 127, 255}));
  connect(pid_presure_H2.CS, valve_pressure_an.opening) annotation(
    Line(points = {{-100, 62}, {-100, 53}}, color = {0, 0, 127}));
  connect(m_flow_coolant.y, SourceWater.m_flow_in) annotation(
    Line(points = {{-54, -108}, {-54, -104}, {-24, -104}, {-24, -88}}, color = {0, 0, 127}));
  connect(air_Compressor.Output, fuelCellStack.port_a_Air) annotation(
    Line(points = {{63, 48}, {50.5, 48}, {50.5, 47}, {52, 47}}, color = {0, 127, 255}));
  connect(constantVoltage.p, air_Compressor.pin_p) annotation(
    Line(points = {{68, 22}, {64, 22}, {64, 34}, {73, 34}, {73, 41}}, color = {0, 0, 255}));
  connect(constantVoltage.n, air_Compressor.pin_n) annotation(
    Line(points = {{82, 22}, {86, 22}, {86, 34}, {77, 34}, {77, 41}}, color = {0, 0, 255}));
  connect(fuelCellStack.OER_Cathode, pid_OER.PV) annotation(
    Line(points = {{52, 62}, {139, 62}, {139, 43}, {127, 43}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_b_Air, valve_pressure_air.port_a) annotation(
    Line(points = {{52, -9}, {69, -9}, {69, -25}}, color = {0, 127, 255}));
  connect(valve_pressure_air.port_b, AirSink.ports[1]) annotation(
    Line(points = {{69, -37}, {69, -62}, {118, -62}}, color = {0, 127, 255}));
  connect(valve_pressure_air.opening, pid_presure_cat.CS) annotation(
    Line(points = {{74, -31}, {93, -31}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_b_Air, pressure_Air.port) annotation(
    Line(points = {{52, -9}, {110, -9}}, color = {0, 127, 255}));
  connect(pressure_Air.p, pid_presure_cat.PV) annotation(
    Line(points = {{117, -17}, {117, -28}, {107, -28}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_b_Coolant, SinkWater.ports[1]) annotation(
    Line(points = {{20, -17}, {20, -33}, {30, -33}, {30, -73}}, color = {0, 127, 255}));
  connect(Temperature_coolant.y, SourceWater.T_in) annotation(
    Line(points = {{-2, -108}, {-2, -104}, {-20, -104}, {-20, -90}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_a_Coolant, SourceWater.ports[1]) annotation(
    Line(points = {{-4, -17}, {-4, -51}, {-16, -51}, {-16, -68}}, color = {0, 127, 255}));
  connect(rampVoltage.p, fuelCellStack.pin_p) annotation(
    Line(points = {{18, 79}, {32, 79}}, color = {0, 0, 255}));
  connect(pump_speed.y, recirculating_pump.controlInterface) annotation(
    Line(points = {{-159, 10.5}, {-109.5, 10.5}, {-109.5, 11}, {-59, 11}}, color = {0, 0, 127}));
  connect(purge_valve_set.y, pid_purge_valve.SP) annotation(
    Line(points = {{-159, -63}, {-77, -63}, {-77, -38}}, color = {0, 0, 127}));
  connect(fuelCellStack.purge_valve_indication, pid_purge_valve.PV) annotation(
    Line(points = {{-26, -35}, {-26, -49}, {-73, -49}, {-73, -38}}, color = {0, 0, 127}));
  connect(pressure_H2.port, valve_pressure_an.port_b) annotation(
    Line(points = {{-70, 81}, {-70, 48}, {-94, 48}}, color = {0, 127, 255}));
  connect(pressure_H2.p, pid_presure_H2.PV) annotation(
    Line(points = {{-78, 88}, {-97, 88}, {-97, 76}}, color = {0, 0, 127}));
  connect(Pressure_H2_set.y, pid_presure_H2.SP) annotation(
    Line(points = {{-158, 89}, {-103, 89}, {-103, 76}}, color = {0, 0, 127}));
  annotation(
    experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 0.1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
    Diagram(coordinateSystem(extent = {{-150, -100}, {150, 100}}, grid = {1, 1})),
    Icon(coordinateSystem(extent = {{-150, -100}, {150, 100}}, grid = {1, 1})));
end TestStack_system;
