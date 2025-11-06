within VirtualFCS.ComponentTesting;

model TestSimpleBLDC
  import SI = Modelica.Units.SI;
  extends Modelica.Icons.Example;

  parameter SI.Voltage Va=500 "Actual armature voltage";
  parameter SI.Time tStart=0.2
    "Start of armature voltage ramp";
  parameter SI.Time tRamp=2 "Armature voltage ramp";
  parameter SI.Torque TLoad=1.8 "Nominal load torque";
  parameter SI.Time tStep=3 "Time of load torque step";
  parameter SI.Inertia JLoad=2e-4
    "Load's moment of inertia";
  Electrical.SimpleBLDC simpleBLDC(Ke = 0.04, Kt = 0.04)  annotation(
    Placement(transformation(origin = {-8, -28}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp ramp(duration = tRamp, height = Va, startTime = tStart) annotation(
    Placement(transformation(origin = {-7, 13}, extent = {{-80, 60}, {-60, 80}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation(
    Placement(transformation(origin = {1, 11}, extent = {{0, 30}, {-20, 50}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {-67, 53}, extent = {{-10, -10}, {10, 10}}, rotation = 270)));
  Modelica.Mechanics.Rotational.Components.Inertia loadInertia(J = JLoad) annotation(
    Placement(transformation(origin = {3, 13}, extent = {{40, -50}, {60, -30}})));
  Modelica.Mechanics.Rotational.Sources.TorqueStep loadTorqueStep(offsetTorque = 0, startTime = tStep, stepTorque = -TLoad, useSupport = false) annotation(
    Placement(transformation(origin = {3, 13}, extent = {{90, -50}, {70, -30}})));
equation
  connect(ramp.y, signalVoltage.v) annotation(
    Line(points = {{-66, 83}, {-66, 63}, {-9, 63}}, color = {0, 0, 255}));
  connect(signalVoltage.n, ground.p) annotation(
    Line(points = {{-19, 51}, {-38, 51}, {-38, 53}, {-57, 53}}, color = {0, 0, 255}));
  connect(loadInertia.flange_b, loadTorqueStep.flange) annotation(
    Line(points = {{63, -27}, {73, -27}}));
  connect(simpleBLDC.shaft, loadInertia.flange_a) annotation(
    Line(points = {{2, -28}, {23, -28}, {23, -26}, {44, -26}}));
  connect(signalVoltage.n, simpleBLDC.n) annotation(
    Line(points = {{-19, 51}, {-12, 51}, {-12, -18}}, color = {0, 0, 255}));
  connect(signalVoltage.p, simpleBLDC.p) annotation(
    Line(points = {{1, 51}, {-4, 51}, {-4, -18}, {-2, -18}}, color = {0, 0, 255}));
annotation(
    experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end TestSimpleBLDC;