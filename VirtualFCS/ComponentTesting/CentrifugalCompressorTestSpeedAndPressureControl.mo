within VirtualFCS.ComponentTesting;

model CentrifugalCompressorTestSpeedAndPressureControl "Model to test the electrical centrifugal compressor with speed and pressure control"
  extends Modelica.Icons.Example;
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {-90, -90}, extent = {{-10, -10}, {10, 10}})));
  replaceable package Medium = Modelica.Media.Air.MoistAir(Temperature(start = system.T_start), AbsolutePressure(start = system.p_start));
  Modelica.Fluid.Sources.Boundary_pT AirSource(redeclare package Medium = Medium, nPorts = 1, p = 1e5) annotation(
    Placement(transformation(origin = {-75, -11}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp1(height = 100000*Modelica.Constants.pi/30, duration = 50, offset = 30000*Modelica.Constants.pi/30, startTime = 100) annotation(
    Placement(transformation(origin = {-78, 31}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(redeclare package Medium = Medium, nPorts = 1, p = 1e5) annotation(
    Placement(transformation(origin = {91, -10}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Valves.ValveLinear valveLinear(dp_nominal = 5000, m_flow_nominal = 0.2, redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {64, -10}, extent = {{-10, -10}, {10, 10}})));
  Fluid.ElectricalCentrifugalCompressor electricalCentrifugalCompressor(redeclare package Medium = Medium)  annotation(
    Placement(transformation(origin = {-17, -11}, extent = {{-17, -17}, {17, 17}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {-50, 60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600) annotation(
    Placement(transformation(origin = {-18, 70}, extent = {{10, -10}, {-10, 10}})));
  VirtualFCS.Control.PID pid1(CSmax = 0, CSmin = 1, Kp = 100, PVmax = 3e5, PVmin = 1e5, Ti = 0.1) annotation(
    Placement(transformation(origin = {49, 32}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Medium) annotation(
    Placement(transformation(origin = {16, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp2(duration = 5, height = 0.5e5, offset = 1.1e5, startTime = 200) annotation(
    Placement(transformation(origin = {14, 36}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(AirSource.ports[1], electricalCentrifugalCompressor.Input) annotation(
    Line(points = {{-65, -11}, {-36, -11}}, color = {0, 127, 255}));
  connect(ramp1.y, electricalCentrifugalCompressor.controlInterface) annotation(
    Line(points = {{-67, 31}, {-44, 31}, {-44, -1}, {-36, -1}}, color = {0, 0, 127}));
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{-50, 70}, {-28, 70}}, color = {0, 0, 255}));
  connect(constantVoltage.n, electricalCentrifugalCompressor.pin_n) annotation(
    Line(points = {{-28, 70}, {-28, 26.5}, {-20, 26.5}, {-20, -1}}, color = {0, 0, 255}));
  connect(electricalCentrifugalCompressor.pin_p, constantVoltage.p) annotation(
    Line(points = {{-14, -1}, {-14, 30.45}, {-7.6, 30.45}, {-7.6, 70.2}}, color = {0, 0, 255}));
  connect(ramp2.y, pid1.SP) annotation(
    Line(points = {{25, 36}, {39, 36}}, color = {0, 0, 127}));
  connect(electricalCentrifugalCompressor.Output, pressure.port) annotation(
    Line(points = {{2, -11}, {10.7, -11}, {10.7, -10}, {16, -10}}, color = {0, 127, 255}));
  connect(pressure.port, valveLinear.port_a) annotation(
    Line(points = {{16, -10}, {54, -10}}, color = {0, 127, 255}));
  connect(pressure.p, pid1.PV) annotation(
    Line(points = {{27, 0}, {35, 0}, {35, 28}, {39, 28}}, color = {0, 0, 127}));
  connect(pid1.CS, valveLinear.opening) annotation(
    Line(points = {{59, 32}, {64, 32}, {64, -2}}, color = {0, 0, 127}));
  connect(valveLinear.port_b, AirSink.ports[1]) annotation(
    Line(points = {{74, -10}, {82, -10}}, color = {0, 127, 255}));
  annotation(
    experiment(StopTime = 500, Interval = 0.5, Tolerance = 1e-06, StartTime = 0),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end CentrifugalCompressorTestSpeedAndPressureControl;