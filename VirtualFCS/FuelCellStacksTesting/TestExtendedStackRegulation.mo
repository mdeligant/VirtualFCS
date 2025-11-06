within VirtualFCS.FuelCellStacksTesting;

model TestExtendedStackRegulation
  extends Modelica.Icons.Example;
  // Medium models
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  //    replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Anode_Medium = Modelica.Media.IdealGases.SingleGases.H2 constrainedby Modelica.Media.Interfaces.PartialSimpleIdealGasMedium;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  Electrochemical.Hydrogen.FuelCellStack fuelCellStack annotation(
    Placement(transformation(origin = {-34, 6}, extent = {{-43, -43}, {43, 43}})));
  Modelica.Fluid.Sources.Boundary_pT H2_pressure_out(nPorts = 1, redeclare package Medium = Anode_Medium, p = 2e5) annotation(
    Placement(transformation(origin = {-92, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_H2(nPorts = 1, m_flow = 0.003, redeclare package Medium = Anode_Medium, T = 293.15) annotation(
    Placement(transformation(origin = {-88, 28}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_water(nPorts = 1, m_flow = 1.5, redeclare package Medium = Coolant_Medium, T = 323.15) annotation(
    Placement(transformation(origin = {-60, -44}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT Water_pressure_out(nPorts = 1, p = 2e5, redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-24, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  inner Modelica.Fluid.System system annotation(
    Placement(transformation(origin = {-90, -91}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sources.Boundary_pT AirSink(nPorts = 1, redeclare package Medium = Cathode_Medium, use_p_in = false, p = 2e5) annotation(
    Placement(transformation(origin = {24, -6}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Ramp resistor_ramp(height = -99.3, duration = 30, offset = 100, startTime = 0) annotation(
    Placement(transformation(origin = {-90, 88}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.VariableResistor resistor1 annotation(
    Placement(transformation(origin = {-34, 64}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sources.MassFlowSource_T Source_massflow_Air(redeclare package Medium = Cathode_Medium, T = 333.15, m_flow = 0.2, nPorts = 1, use_m_flow_in = true, use_T_in = true, X = {0.78, 0.21, 0.01}) annotation(
    Placement(transformation(origin = {17, 26}, extent = {{10, -10}, {-10, 10}})));
  Control.PID pid_OER(CSmax = 2, CSmin = 0.01, Kp = 1, PVmax = 5, PVmin = 0, Ti = 0.1) annotation(
    Placement(transformation(origin = {59, 51}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Constant const1(k = 1.5) annotation(
    Placement(transformation(origin = {90, 56}, extent = {{10, -10}, {-10, 10}})));
  Control.PID pid_phi(CSmax = 274, CSmin = 400, Kp = 1, PVmax = 2, PVmin = 0, Ti = 1, CSs(start = 0.6315, fixed = true)) annotation(
    Placement(transformation(origin = {52, 19}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.Constant const11(k = 0.80) annotation(
    Placement(transformation(origin = {89, 22}, extent = {{10, -10}, {-10, 10}})));
equation
  connect(Source_massflow_H2.ports[1], fuelCellStack.port_a_H2) annotation(
    Line(points = {{-78, 28}, {-68, 28}, {-68, 26}, {-66, 26}}, color = {0, 127, 255}));
  connect(H2_pressure_out.ports[1], fuelCellStack.port_b_H2) annotation(
    Line(points = {{-82, -14}, {-66, -14}}, color = {0, 127, 255}));
  connect(Water_pressure_out.ports[1], fuelCellStack.port_b_Coolant) annotation(
    Line(points = {{-24, -52}, {-24, -38}, {-25, -38}, {-25, -20}}, color = {0, 127, 255}));
  connect(fuelCellStack.port_a_Coolant, Source_massflow_water.ports[1]) annotation(
    Line(points = {{-43, -20}, {-43, -43.7571}, {-49.2, -43.7571}}, color = {0, 127, 255}));
  connect(fuelCellStack.pin_n, resistor1.n) annotation(
    Line(points = {{-51, 49}, {-51, 64}, {-44, 64}}, color = {0, 0, 255}));
  connect(fuelCellStack.port_a_Air, Source_massflow_Air.ports[1]) annotation(
    Line(points = {{-2.46667, 26.0667}, {3.03331, 26.0667}, {3.03331, 26}, {7, 26}}, color = {0, 127, 255}));
  connect(AirSink.ports[1], fuelCellStack.port_b_Air) annotation(
    Line(points = {{14, -6}, {6, -6}, {6, -14}, {-2, -14}}, color = {0, 127, 255}));
  connect(pid_OER.CS, Source_massflow_Air.m_flow_in) annotation(
    Line(points = {{49, 51}, {41, 51}, {41, 34}, {27, 34}}, color = {0, 0, 127}));
  connect(pid_phi.CS, Source_massflow_Air.T_in) annotation(
    Line(points = {{42, 19}, {40, 19}, {40, 30}, {29, 30}}, color = {0, 0, 127}));
  connect(pid_phi.SP, const11.y) annotation(
    Line(points = {{62, 23}, {78, 23}, {78, 22}}, color = {0, 0, 127}));
  connect(fuelCellStack.OER, pid_OER.PV) annotation(
    Line(points = {{-4, -20}, {116, -20}, {116, 40}, {76, 40}, {76, 48}, {70, 48}}, color = {0, 0, 127}));
  connect(fuelCellStack.phi, pid_phi.PV) annotation(
    Line(points = {{-4, -26}, {74, -26}, {74, 15}, {62, 15}}, color = {0, 0, 127}));
  connect(pid_OER.SP, const1.y) annotation(
    Line(points = {{70, 56}, {79, 56}}, color = {0, 0, 127}));
  connect(resistor1.p, fuelCellStack.pin_p) annotation(
    Line(points = {{-24, 64}, {-16, 64}, {-16, 50}}, color = {0, 0, 255}));
  connect(resistor_ramp.y, resistor1.R) annotation(
    Line(points = {{-78, 88}, {-34, 88}, {-34, 76}}, color = {0, 0, 127}));
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
end TestExtendedStackRegulation;