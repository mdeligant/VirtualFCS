within VirtualFCS.FuelCellStacksTesting;

model TestStackWithTwoStageCompressorControlOER
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Modelica.Media.IdealGases.SingleGases.H2 constrainedby Modelica.Media.Interfaces.PartialSimpleIdealGasMedium;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack annotation(
    Placement(transformation(origin = {-16.6667, -2}, extent = {{-37.6667, -37.6667}, {37.6667, 37.6667}})));
  Modelica.Fluid.Sources.Boundary_pT H2_pressure_out(nPorts = 1, redeclare package Medium = Anode_Medium, p = 2e5) annotation(
    Placement(transformation(origin = {-86, -20}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_H2(nPorts = 1, m_flow = 0.003, redeclare package Medium = Anode_Medium, T = 293.15) annotation(
    Placement(transformation(origin = {-84, 16}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_water(nPorts = 1, m_flow = 1.5, redeclare package Medium = Coolant_Medium, T = 323.15) annotation(
    Placement(transformation(origin = {-56, -44}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT Water_pressure_out(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {-90, 89}, extent = {{-10, -10}, {10, 10}})));
  Control.PID pid(CSmax = 150000*Modelica.Constants.pi/30, CSmin = 0, Kp = 1, PVmax = 3e5, PVmin = 1e5, Ti = 0.1) annotation(
    Placement(transformation(origin = {74, 91}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp pressure_ramp(duration = 20, height = 0.8e5, offset = 1.2e5, startTime = 50) annotation(
    Placement(transformation(origin = {32, 95}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 0.11) annotation(
    Placement(transformation(origin = {446, 17}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSource(nPorts = 1, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {144, 34}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Cathode_Medium, dp_nominal = 5000, m_flow_nominal = 0.2) annotation(
    Placement(transformation(origin = {50, -20}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(nPorts = 1, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {82, -20}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {32, 62}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp resistor_ramp(height = 0, duration = 30, offset = 0.7, startTime = 10) annotation(
    Placement(transformation(origin = {-60, 90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.VariableResistor resistor1 annotation(
    Placement(transformation(origin = {-18, 66}, extent = {{10, -10}, {-10, 10}})));
  Control.PID pid1(CSmax = 1, CSmin = 0, Kp = 1, PVmax = 4, PVmin = 0.5, Ti = 0.1) annotation(
    Placement(transformation(origin = {72, -60}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp OER_ramp(duration = 5, height = 0.5, offset = 2.5, startTime = 200) annotation(
    Placement(transformation(origin = {119, -56}, extent = {{10, -10}, {-10, 10}})));
  Fluid.ElectricalCentrifugalCompressor electricalCentrifugalCompressor(redeclare package Medium = Cathode_Medium)  annotation(
    Placement(transformation(origin = {60, 34}, extent = {{10, 10}, {-10, -10}})));
  Fluid.ElectricalCentrifugalCompressor electricalCentrifugalCompressor1(redeclare package Medium = Cathode_Medium)  annotation(
    Placement(transformation(origin = {108, 34}, extent = {{10, 10}, {-10, -10}}, rotation = -0)));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {189, -2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600) annotation(
    Placement(transformation(origin = {125, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
equation
  connect(Source_massflow_H2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-74, 16}, {-44, 16}}, color = {0, 127, 255}));
  connect(H2_pressure_out.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-76, -20}, {-44, -20}}, color = {0, 127, 255}));
  connect(Water_pressure_out.ports[1], fuelCellStack.port_b_Coolant) annotation(
    Line(points = {{-12, -52}, {-12, -25}, {-9, -25}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, Source_massflow_water.ports[1]) annotation(
    Line(points = {{-24, -25}, {-24, -43.9571}, {-45.6, -43.9571}}, color = {0, 127, 255}));
  connect(pressure_ramp.y, pid.SP) annotation(
    Line(points = {{43, 95}, {64, 95}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_b_Air, valveLinear.port_a) annotation(
    Line(points = {{11, -20}, {40, -20}}, color = {0, 127, 255}));
  connect(valveLinear.port_b, AirSink.ports[1]) annotation(
    Line(points = {{60, -20}, {72, -20}}, color = {0, 127, 255}));
  connect(fuelCellStack.pin_n, resistor1.n) annotation(
    Line(points = {{-32, 36}, {-32, 65}, {-28, 65}, {-28, 66}}, color = {0, 0, 255}));
  connect(resistor1.p, fuelCellStack.pin_p) annotation(
    Line(points = {{-8, 66}, {-8, 36}, {-2, 36}}, color = {0, 0, 255}));
  connect(pressure.p, pid.PV) annotation(
    Line(points = {{43, 62}, {43, 87}, {64, 87}}, color = {0, 0, 127}));
  connect(valveLinear.opening, pid1.CS) annotation(
    Line(points = {{50, -28}, {50, -60}, {62, -60}}, color = {0, 0, 127}));
  connect(OER_ramp.y, pid1.SP) annotation(
    Line(points = {{108, -56}, {82, -56}}, color = {0, 0, 127}));
  connect(electricalCentrifugalCompressor.Input, electricalCentrifugalCompressor1.Output) annotation(
    Line(points = {{71, 34}, {97, 34}}, color = {0, 127, 255}));
  connect(electricalCentrifugalCompressor1.Input, AirSource.ports[1]) annotation(
    Line(points = {{119, 34}, {134, 34}}, color = {0, 127, 255}));
  connect(pid.CS, electricalCentrifugalCompressor1.controlInterface) annotation(
    Line(points = {{84, 92}, {128, 92}, {128, 28}, {119, 28}}, color = {0, 0, 127}));
  connect(pid.CS, electricalCentrifugalCompressor.controlInterface) annotation(
    Line(points = {{84, 92}, {92, 92}, {92, 28}, {71, 28}}, color = {0, 0, 127}));
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{189, 8}, {157, 8}, {157, 10}, {125, 10}}, color = {0, 0, 255}));
  connect(constantVoltage.n, electricalCentrifugalCompressor1.pin_n) annotation(
    Line(points = {{126, 10}, {126, 22}, {110, 22}, {110, 28}}, color = {0, 0, 255}));
  connect(constantVoltage.n, electricalCentrifugalCompressor.pin_n) annotation(
    Line(points = {{126, 10}, {126, 22}, {62, 22}, {62, 28}}, color = {0, 0, 255}));
  connect(constantVoltage.p, electricalCentrifugalCompressor1.pin_p) annotation(
    Line(points = {{126, -10}, {106, -10}, {106, 28}}, color = {0, 0, 255}));
  connect(constantVoltage.p, electricalCentrifugalCompressor.pin_p) annotation(
    Line(points = {{126, -10}, {106, -10}, {106, 18}, {58, 18}, {58, 28}}, color = {0, 0, 255}));
connect(fuelCellStack.OER, pid1.PV) annotation(
    Line(points = {{5, -25}, {34, -25}, {34, -78}, {92, -78}, {92, -64}, {82, -64}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_a_Air, pressure.port) annotation(
    Line(points = {{10, 16}, {32, 16}, {32, 52}}, color = {0, 127, 255}));
  connect(electricalCentrifugalCompressor.Output, fuelCellStack.port_a_Air) annotation(
    Line(points = {{49, 34}, {32, 34}, {32, 16}, {10, 16}}, color = {0, 127, 255}));
  connect(resistor_ramp.y, resistor1.R) annotation(
    Line(points = {{-48, 90}, {-18, 90}, {-18, 78}}, color = {0, 0, 127}));
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
    experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 2),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"),
  Diagram(coordinateSystem(extent = {{-100, -100}, {200, 150}})),
  Icon(coordinateSystem(extent = {{-100, -100}, {200, 150}})));
end TestStackWithTwoStageCompressorControlOER;