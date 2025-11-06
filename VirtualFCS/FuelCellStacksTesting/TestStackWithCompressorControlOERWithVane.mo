within VirtualFCS.FuelCellStacksTesting;

model TestStackWithCompressorControlOERWithVane

  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Modelica.Media.IdealGases.SingleGases.H2 constrainedby Modelica.Media.Interfaces.PartialSimpleIdealGasMedium;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;

  extends Modelica.Icons.Example;
  // Medium declaration
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack annotation(
    Placement(transformation(origin = {-26, -12.5}, extent = {{-35, -35}, {35, 35}})));
  Modelica.Fluid.Sources.Boundary_pT H2_pressure_out(nPorts = 1, redeclare package Medium = Anode_Medium, p = 2e5) annotation(
    Placement(transformation(origin = {-84, -28}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_H2(nPorts = 1, m_flow = 0.02, redeclare package Medium = Anode_Medium, T = 293.15) annotation(
    Placement(transformation(origin = {-82, 4}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_water(nPorts = 1, m_flow = 0.7, redeclare package Medium = Coolant_Medium, T = 323.15) annotation(
    Placement(transformation(origin = {-56, -66}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT Water_pressure_out(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-18, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {-90, -89}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp pressure_ramp(duration = 5, height = 0.3e5, offset = 1.3e5, startTime = 10) annotation(
    Placement(transformation(origin = {-70, 85}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 0.11) annotation(
    Placement(transformation(origin = {446, 17}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Cathode_Medium, dp_nominal = 5000, m_flow_nominal = 0.2, allowFlowReversal = false) annotation(
    Placement(transformation(origin = {30, -28}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(nPorts = 1, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {62, -28}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {10, 72}, extent = {{-10, -10}, {10, 10}})));
  Control.PID pid_OER(CSmax = 1, CSmin = 0, Kp = 1, PVmax = 4, PVmin = 0, Ti = 0.1) annotation(
    Placement(transformation(origin = {50, -72}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp OER_ramp(duration = 5, height = 2, offset = 1.5, startTime = 1) annotation(
    Placement(transformation(origin = {87, -68}, extent = {{10, -10}, {-10, 10}})));
  Control.PID pid_pressure(CSmax = 150000*Modelica.Constants.pi/30, CSmin = 0, Kp = 1, PVmax = 3e5, PVmin = 1e5, Ti = 0.1) annotation(
    Placement(transformation(origin = {50, 81}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSource(nPorts = 1, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {90, 4}, extent = {{10, -10}, {-10, 10}})));
  Fluid.ElectricalCentrifugalCompressor electricalCentrifugalCompressor(flow_scale_factor = 3, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {36, 4}, extent = {{16, -16}, {-16, 16}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {56, 38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600) annotation(
    Placement(transformation(origin = {38, 48}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.RampCurrent rampCurrent(I = 390, duration = 80, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {-26, 50}, extent = {{10, -10}, {-10, 10}})));
equation
  connect(Source_massflow_H2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-72, 4}, {-52, 4}}, color = {0, 127, 255}));
  connect(H2_pressure_out.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-74, -28}, {-63, -28}, {-63, -29}, {-52, -29}}, color = {0, 127, 255}));
  connect(Water_pressure_out.ports[1], fuelCellStack.port_b_Coolant) annotation(
    Line(points = {{-18, -74}, {-18, -53}, {-19, -53}, {-19, -33.5}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, Source_massflow_water.ports[1]) annotation(
    Line(points = {{-33, -33.5}, {-33, -65.9571}, {-45.6, -65.9571}}, color = {0, 127, 255}));
  connect(pressure_ramp.y, pid_pressure.SP) annotation(
    Line(points = {{-59, 85}, {40, 85}}, color = {0, 0, 127}));
  connect(valveLinear.port_b, AirSink.ports[1]) annotation(
    Line(points = {{40, -28}, {52, -28}}, color = {0, 127, 255}));
  connect(valveLinear.opening, pid_OER.CS) annotation(
    Line(points = {{30, -36}, {30, -72}, {40, -72}}, color = {0, 0, 127}));
  connect(OER_ramp.y, pid_OER.SP) annotation(
    Line(points = {{76, -68}, {60, -68}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_a_Air, electricalCentrifugalCompressor.Output) annotation(
    Line(points = {{0, 4}, {18, 4}}, color = {0, 127, 255}));
  connect(electricalCentrifugalCompressor.Input, AirSource.ports[1]) annotation(
    Line(points = {{54, 4}, {80, 4}}, color = {0, 127, 255}));
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{56, 48}, {48, 48}}, color = {0, 0, 255}));
  connect(fuelCellStack.port_b_Air, valveLinear.port_a) annotation(
    Line(points = {{0, -29}, {15, -29}, {15, -28}, {20, -28}}, color = {0, 127, 255}));
  connect(fuelCellStack.pin_n, rampCurrent.n) annotation(
    Line(points = {{-40, 22.5}, {-40, 50}, {-36, 50}}, color = {0, 0, 255}));
  connect(fuelCellStack.OER, pid_OER.PV) annotation(
    Line(points = {{-5, -33.5}, {10, -33.5}, {10, -94}, {66, -94}, {66, -76}, {60, -76}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_a_Air, pressure.port) annotation(
    Line(points = {{0, 4}, {10, 4}, {10, 62}}, color = {0, 127, 255}));
  connect(pressure.p, pid_pressure.PV) annotation(
    Line(points = {{21, 72}, {29, 72}, {29, 77}, {39, 77}}, color = {0, 0, 127}));
  connect(electricalCentrifugalCompressor.pin_n, constantVoltage.n) annotation(
    Line(points = {{39, 14}, {40, 14}, {40, 34}, {48, 34}, {48, 48}}, color = {0, 0, 255}));
  connect(electricalCentrifugalCompressor.pin_p, constantVoltage.p) annotation(
    Line(points = {{33, 14}, {34, 14}, {34, 36}, {28, 36}, {28, 48}}, color = {0, 0, 255}));
  connect(electricalCentrifugalCompressor.controlInterface, pid_pressure.CS) annotation(
    Line(points = {{54, 14}, {70, 14}, {70, 81}, {60, 81}}, color = {0, 0, 127}));
  connect(fuelCellStack.pin_p, rampCurrent.p) annotation(
    Line(points = {{-12, 22.5}, {-12, 50}, {-16, 50}}, color = {0, 0, 255}));
  annotation(
    Documentation(info = "<html><head></head><body>
  
  This example tests a fuel cell stack with constant hydrogen, air and cooling supply. The inlet and outlet conditions are presented in the table below.
  
  
  <html>
  <head>
  
    <title>Simulation Parameters</title>
    <style>
      body {
        font-family: Arial, sans-serif;
        margin: 40px;
        background-color: #f9f9f9;
      }
      table {
        width: 80%;
        border-collapse: collapse;
        margin: 20px auto;
        background-color: #fff;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
      }
      caption {
        font-size: 1.4em;
        margin-bottom: 10px;
        font-weight: bold;
      }
      th, td {
        border: 1px solid #ccc;
        text-align: center;
        padding: 10px;
      }
      th {
        background-color: #e0f0ff;
      }
      tr:nth-child(even) {
        background-color: #f5faff;
      }
    </style>
  </head>
  <body>
  
    <table>
      <thead>
        <tr>
          <th>Circuit</th>
          <th>Inlet Temperature (°C)</th>
          <th>Inlet Flow Rate (kg/s)</th>
          <th>Outlet Pressure (bar)</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>Hydrogen</td>
          <td>20</td>
          <td>0.003</td>
          <td>2.0</td>
        </tr>
  
        <tr>
          <td>Air</td>
          <td>50</td>
          <td>0.25</td>
          <td>2.0</td>
        </tr>
        
  
        <tr>
          <td>Cooling water</td>
          <td>50</td>
          <td>1.5</td>
          <td>2.0</td>
        </tr>
  
      </tbody>
    </table>
  With these conditions, the fuel cell stack is connected to a fixed resistor of 0.7 Ohm, giving an output power of 124 kW with a stack temperature of 81.38 °C.
  
  </body>
  </html>
  
  
  
  
  
  </div></div></body></html>"),
    experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end TestStackWithCompressorControlOERWithVane;