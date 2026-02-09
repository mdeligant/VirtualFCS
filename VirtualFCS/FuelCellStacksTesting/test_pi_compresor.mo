within VirtualFCS.FuelCellStacksTesting;

model test_pi_compresor
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Media.MoistHydrogenThreeComponents;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Modelica.Fluid.Sources.Boundary_pT AirSink(redeclare package Medium = Cathode_Medium, nPorts = 1, p = 1.1e5, use_p_in = true) annotation(
    Placement(transformation(origin = {60, -56}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp air_OER_ramp(duration = 0.01, height = 0.5, offset = 3.5, startTime = 50) annotation(
    Placement(transformation(origin = {157, 30}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSource(redeclare package Medium = Cathode_Medium, nPorts = 1, X = {0.2, 0.02, 0.78}) annotation(
    Placement(transformation(origin = {90, 50}, extent = {{10, -10}, {-10, 10}})));
  Fluid.ElectricalCentrifugalCompressor electricalCentrifugalCompressor(redeclare package Medium = Cathode_Medium, flow_scale_factor = 1) annotation(
    Placement(transformation(origin = {61, 47}, extent = {{11, 11}, {-11, -11}}, rotation = -0)));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {82, 8}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600) annotation(
    Placement(transformation(origin = {57, 19}, extent = {{-7, -7}, {7, 7}})));
  Control.PID pid_OER(CSmax = 150000*Modelica.Constants.pi/30, CSmin = 15000*Modelica.Constants.pi/30/10, Kp = 0.01, PVmax = 10, PVmin = 0, Ti = 0.001, CSs(start = 0.5, fixed = true)) annotation(
    Placement(transformation(origin = {91, 32}, extent = {{5, 5}, {-5, -5}}, rotation = -0)));
  Modelica.Fluid.Pipes.DynamicPipe channelCathode(redeclare model FlowModel = Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow, redeclare model HeatTransfer = Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer, redeclare package Medium = Cathode_Medium, T_start = 333.15, crossArea = 0.204/10^6, diameter = 0.0002, isCircular = false, length = 0.15, modelStructure = Modelica.Fluid.Types.ModelStructure.a_v_b, nNodes = 2, nParallel = 370*72, p_a_start = 101000, perimeter = 1.6/1000, use_HeatTransfer = true) annotation(
    Placement(transformation(origin = {12, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {12, 56}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain gain(k = 1/(0.032/(96485*4)*370*100)) annotation(
    Placement(transformation(origin = {84, 82}, extent = {{-10, -10}, {10, 10}})));
  VirtualFCS.Fluid.Membrane GDL_ca(redeclare package Medium = Cathode_Medium, T_start = 353.15, V = 370*72*0.15*0.204/10^6 + 370*237/1000*250/10^6/2*0.5, p_start = 102000) annotation(
    Placement(transformation(origin = {11, -1}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Blocks.Sources.Ramp air_OER_ramp1(duration = 5, height = 0.2e5, offset = 1e5, startTime = 10) annotation(
    Placement(transformation(origin = {110, -51}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Math.Add add annotation(
    Placement(transformation(origin = {118, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
  Modelica.Blocks.Sources.Ramp air_OER_ramp2(duration = 0.01, height = -0.5, offset = 0, startTime = 150) annotation(
    Placement(transformation(origin = {165, -2}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor[2](each C = 10, each T(fixed = true, start = 333.15), each der_T(fixed = false)) annotation(
    Placement(transformation(origin = {-25, -21}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Pipes.DynamicPipe channelCathode1(redeclare model HeatTransfer = Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer, redeclare package Medium = Cathode_Medium, T_start = 333.15, crossArea = 0.204/10^6, diameter = 0.0002, isCircular = false, length = 0.15, modelStructure = Modelica.Fluid.Types.ModelStructure.a_v_b, nNodes = 1, nParallel = 370*72, p_a_start = 103000, perimeter = 1.6/1000, use_HeatTransfer = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, massDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, momentumDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, p_b_start = 103000, m_flow_start = 0.2) annotation(
    Placement(transformation(origin = {13, 27}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor1[1](each C = 10, each T(fixed = true, start = 333.15), each der_T(fixed = false)) annotation(
    Placement(transformation(origin = {-25, 36}, extent = {{-10, -10}, {10, 10}})));
equation
  GDL_ca.mass_flow_b = {0, 0, 0};
  GDL_ca.mass_flow_a = {0, 0, 0};
  GDL_ca.Q_flow = 0;
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{82, 18}, {73, 18}, {73, 19}, {64, 19}}, color = {0, 0, 255}));
  connect(AirSource.ports[1], electricalCentrifugalCompressor.Input) annotation(
    Line(points = {{80, 50}, {74, 50}, {74, 48}}, color = {0, 127, 255}));
  connect(constantVoltage.p, electricalCentrifugalCompressor.pin_p) annotation(
    Line(points = {{50, 19}, {50, 34}, {60, 34}, {60, 40}, {58, 40}}, color = {0, 0, 255}));
  connect(constantVoltage.n, electricalCentrifugalCompressor.pin_n) annotation(
    Line(points = {{64, 19}, {66, 19}, {66, 36}, {62, 36}, {62, 40}}, color = {0, 0, 255}));
  connect(pid_OER.CS, electricalCentrifugalCompressor.controlInterface) annotation(
    Line(points = {{86, 32}, {80, 32}, {80, 40}, {74, 40}}, color = {0, 0, 127}));
  connect(gain.u, massFlowRate.m_flow) annotation(
    Line(points = {{72, 82}, {72, 84}, {1, 84}, {1, 56}}, color = {0, 0, 127}));
  connect(gain.y, pid_OER.PV) annotation(
    Line(points = {{96, 82}, {114, 82}, {114, 34}, {96, 34}}, color = {0, 0, 127}));
  connect(electricalCentrifugalCompressor.Output, massFlowRate.port_a) annotation(
    Line(points = {{48, 48}, {48, 66}, {12, 66}}, color = {0, 127, 255}));
  connect(GDL_ca.port_2, channelCathode.port_a) annotation(
    Line(points = {{12, -10}, {12, -20}}, color = {0, 127, 255}));
  connect(channelCathode.port_b, AirSink.ports[1]) annotation(
    Line(points = {{12, -40}, {12, -56}, {50, -56}}, color = {0, 127, 255}));
  connect(air_OER_ramp1.y, AirSink.p_in) annotation(
    Line(points = {{100, -50}, {72, -50}, {72, -48}}, color = {0, 0, 127}));
  connect(air_OER_ramp.y, add.u2) annotation(
    Line(points = {{146, 30}, {132, 30}, {132, 18}, {130, 18}}, color = {0, 0, 127}));
  connect(add.y, pid_OER.SP) annotation(
    Line(points = {{108, 12}, {102, 12}, {102, 30}, {96, 30}}, color = {0, 0, 127}));
  connect(air_OER_ramp2.y, add.u1) annotation(
    Line(points = {{154, -2}, {130, -2}, {130, 6}}, color = {0, 0, 127}));
  connect(heatCapacitor.port, channelCathode.heatPorts) annotation(
    Line(points = {{-24, -30}, {8, -30}}, color = {191, 0, 0}, thickness = 0.5));
  connect(channelCathode1.port_b, GDL_ca.port_1) annotation(
    Line(points = {{14, 18}, {12, 18}, {12, 10}}, color = {0, 127, 255}));
  connect(channelCathode1.port_a, massFlowRate.port_b) annotation(
    Line(points = {{14, 38}, {12, 38}, {12, 46}}, color = {0, 127, 255}));
  connect(heatCapacitor1.port, channelCathode1.heatPorts) annotation(
    Line(points = {{-24, 26}, {8, 26}}, color = {191, 0, 0}, thickness = 0.5));
  annotation(
    experiment(StartTime = 0, StopTime = 200, Tolerance = 1e-06, Interval = 0.01),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end test_pi_compresor;
