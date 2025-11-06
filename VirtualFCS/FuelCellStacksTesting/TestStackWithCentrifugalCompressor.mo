within VirtualFCS.FuelCellStacksTesting;

model TestStackWithCentrifugalCompressor
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Modelica.Media.IdealGases.SingleGases.H2 constrainedby Modelica.Media.Interfaces.PartialSimpleIdealGasMedium;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack annotation(
    Placement(transformation(origin = {-23.3333, 4}, extent = {{-36.3333, -54.5}, {36.3333, 54.5}})));
  Modelica.Fluid.Sources.Boundary_pT H2_pressure_out(nPorts = 1, redeclare package Medium = Anode_Medium, p = 2e5) annotation(
    Placement(transformation(origin = {-86, -16}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_H2(nPorts = 1, m_flow = 0.003, redeclare package Medium = Anode_Medium, T = 293.15) annotation(
    Placement(transformation(origin = {-84, 28}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_water(nPorts = 1, m_flow = 1.5, redeclare package Medium = Coolant_Medium, T = 323.15) annotation(
    Placement(transformation(origin = {-56, -44}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT Water_pressure_out(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-12, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {-90, -91}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Ramp speed_ramp(duration = 5, height = 60000*Modelica.Constants.pi/30, offset = 60000*Modelica.Constants.pi/30, startTime = 10) annotation(
    Placement(transformation(origin = {86, 71}, extent = {{10, -10}, {-10, 10}})));
  Control.PID pid1(CSmax = 0, CSmin = 1, Kp = 100, PVmax = 3e5, PVmin = 1e5, Ti = 0.1) annotation(
    Placement(transformation(origin = {53, -49}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp pressure_ramp(duration = 5, height = 0.5e5, offset = 1e5, startTime = 200) annotation(
    Placement(transformation(origin = {86, -43}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 0.11) annotation(
    Placement(transformation(origin = {446, 17}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSource(nPorts = 1, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {90, 30}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Valves.ValveLinear valveLinear(redeclare package Medium = Cathode_Medium, dp_nominal = 5000, m_flow_nominal = 0.2) annotation(
    Placement(transformation(origin = {38, -18}, extent = {{-10, 10}, {10, -10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(nPorts = 1, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {70, -18}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {18, -72}, extent = {{-10, 10}, {10, -10}})));
  Fluid.ElectricalCentrifugalCompressor electricalCentrifugalCompressor(flow_scale_factor = 3, redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {38, 30}, extent = {{16, -16}, {-16, 16}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(transformation(origin = {62, 58}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V = 600) annotation(
    Placement(transformation(origin = {38, 68}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Sources.RampCurrent rampCurrent(I = 390, duration = 80, offset = 10, startTime = 10) annotation(
    Placement(transformation(origin = {-23, 82}, extent = {{10, -10}, {-10, 10}})));
equation
  connect(Source_massflow_H2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-74, 28}, {-60.5, 28}, {-60.5, 29}, {-50, 29}}, color = {0, 127, 255}));
  connect(H2_pressure_out.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-76, -16}, {-63, -16}, {-63, -21}, {-50, -21}}, color = {0, 127, 255}));
  connect(Water_pressure_out.ports[1], fuelCellStack.port_b_Coolant) annotation(
    Line(points = {{-12, -52}, {-12, -39.5}, {-16, -39.5}, {-16, -29}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, Source_massflow_water.ports[1]) annotation(
    Line(points = {{-31, -29}, {-31, -43.9571}, {-45.6, -43.9571}}, color = {0, 127, 255}));
  connect(pressure_ramp.y, pid1.SP) annotation(
    Line(points = {{75, -43}, {71, -43}, {71, -45}, {63, -45}}, color = {0, 0, 127}));
  connect(valveLinear.port_b, AirSink.ports[1]) annotation(
    Line(points = {{48, -18}, {60, -18}}, color = {0, 127, 255}));
  connect(valveLinear.opening, pid1.CS) annotation(
    Line(points = {{38, -26}, {38, -49}, {43, -49}}, color = {0, 0, 127}));
  connect(fuelCellStack.port_a_Air, electricalCentrifugalCompressor.Output) annotation(
    Line(points = {{3, 29}, {18, 29}, {18, 30}, {20, 30}}, color = {0, 127, 255}));
  connect(electricalCentrifugalCompressor.Input, AirSource.ports[1]) annotation(
    Line(points = {{55.6, 30}, {79.6, 30}}, color = {0, 127, 255}));
  connect(speed_ramp.y, electricalCentrifugalCompressor.controlInterface) annotation(
    Line(points = {{75, 71}, {72, 71}, {72, 40}, {56, 40}}, color = {0, 0, 127}));
  connect(constantVoltage.p, electricalCentrifugalCompressor.pin_p) annotation(
    Line(points = {{28, 68}, {28, 52}, {35, 52}, {35, 40}}, color = {0, 0, 255}));
  connect(electricalCentrifugalCompressor.pin_n, constantVoltage.n) annotation(
    Line(points = {{41, 40}, {41, 52}, {48, 52}, {48, 68}}, color = {0, 0, 255}));
  connect(ground.p, constantVoltage.n) annotation(
    Line(points = {{62, 68}, {48, 68}}, color = {0, 0, 255}));
  connect(valveLinear.port_a, fuelCellStack.port_b_Air) annotation(
    Line(points = {{28, -18}, {19, -18}, {19, -21}, {3, -21}}, color = {0, 127, 255}));
  connect(pressure.p, pid1.PV) annotation(
    Line(points = {{30, -72}, {68, -72}, {68, -52}, {64, -52}}, color = {0, 0, 127}));
  connect(pressure.port, fuelCellStack.port_b_Air) annotation(
    Line(points = {{18, -62}, {18, -21}, {3, -21}}, color = {0, 127, 255}));
  connect(rampCurrent.p, fuelCellStack.pin_p) annotation(
    Line(points = {{-12, 82}, {-8, 82}, {-8, 58}}, color = {0, 0, 255}));
  connect(rampCurrent.n, fuelCellStack.pin_n) annotation(
    Line(points = {{-32, 82}, {-38, 82}, {-38, 58}}, color = {0, 0, 255}));
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
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end TestStackWithCentrifugalCompressor;