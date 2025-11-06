within VirtualFCS.ComponentTesting;

model CentrifugalCompressorTest "Model to test the electrical centrifugal compressor"
  extends Modelica.Icons.Example;
  inner Modelica.Fluid.System system annotation(
    Placement(visible = true, transformation(origin = {90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  replaceable package Medium = Modelica.Media.Air.MoistAir(Temperature(start = system.T_start), AbsolutePressure(start = system.p_start));
  Fluid.ElectricalCentrifugalCompressor compressor(flow_scale_factor = 1, redeclare package Medium = Medium)   annotation(
    Placement(transformation(origin = {-4, 4}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Fluid.Sources.Boundary_pT AirSource(nPorts = 1, redeclare package Medium = Modelica.Media.Air.MoistAir, p = 1e5)  annotation(
    Placement(transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(nPorts = 1, redeclare package Medium = Modelica.Media.Air.MoistAir, use_p_in = true)  annotation(
    Placement(transformation(origin = {60, 0}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(height = 0.5e5, duration = 200, offset = 1.1e5, startTime = 120)  annotation(
    Placement(transformation(origin = {90, 8}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 110000*3.141592/30, offset = 20000*3.141592/30, startTime = 0) annotation(
    Placement(transformation(origin = {-86, 62}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {-32, 64}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600)  annotation(
    Placement(transformation(origin = {-4, 76}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
equation
  connect(AirSource.ports[1], compressor.Input) annotation(
    Line(points = {{-74, 0}, {-50, 0}, {-50, 4}, {-26, 4}}, color = {0, 127, 255}));
  connect(compressor.Output, AirSink.ports[1]) annotation(
    Line(points = {{18, 4}, {34, 4}, {34, 0}, {50, 0}}, color = {0, 127, 255}));
  connect(ramp.y, AirSink.p_in) annotation(
    Line(points = {{79, 8}, {72, 8}}, color = {0, 0, 127}));
  connect(ramp1.y, compressor.controlInterface) annotation(
    Line(points = {{-75, 62}, {-48, 62}, {-48, 16}, {-27, 16}}, color = {0, 0, 127}));
  connect(constantVoltage.p, compressor.pin_p) annotation(
    Line(points = {{6, 76}, {6, 44}, {0, 44}, {0, 16}}, color = {0, 0, 255}));
  connect(constantVoltage.n, compressor.pin_n) annotation(
    Line(points = {{-14, 76}, {-14, 45}, {-7, 45}, {-7, 16}}, color = {0, 0, 255}));
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{-32, 74}, {-32, 76}, {-14, 76}}, color = {0, 0, 255}));
  annotation(
    experiment(StopTime = 12000, Interval = 0.5, Tolerance = 1e-06, StartTime = 0),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
  __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end CentrifugalCompressorTest;