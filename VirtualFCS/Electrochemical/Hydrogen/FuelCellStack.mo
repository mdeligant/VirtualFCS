within VirtualFCS.Electrochemical.Hydrogen;

model FuelCellStack
  //*** DEFINE REPLACEABLE PACKAGES ***//
  // System
  outer Modelica.Fluid.System system "System properties";
  // Medium models
  //  replaceable package Cathode_Medium = Modelica.Media.Air.MoistAir;
  replaceable package Cathode_Medium = Media.MoistAirThreeComponents;
  replaceable package Anode_Medium = Media.MoistHydrogenThreeComponents;
  replaceable package Coolant_Medium = Modelica.Media.Water.ConstantPropertyLiquidWater constrainedby Modelica.Media.Interfaces.PartialMedium;
  // Component indice
  parameter Integer i_O2 = 1;
  parameter Integer i_H2 = 1;
  parameter Integer i_H2O = 2;
  parameter Integer i_N2 = 3;
  //*** DECLARE PARAMETERS ***//
  // Physical parameters
  // Fuel Cell Stack Paramters
  parameter Real m_FC_stack(unit = "kg") = 42 "FC stack mass";
  parameter Real L_FC_stack(unit = "m") = 0.420 "FC stack length";
  parameter Real W_FC_stack(unit = "m") = 0.582 "FC stack width";
  parameter Real H_FC_stack(unit = "m") = 0.156 "FC stack height";
  parameter Real vol_FC_stack(unit = "m3") = L_FC_stack*W_FC_stack*H_FC_stack "FC stack volume";
  parameter Real I_rated_FC_stack(unit = "A") = 450 "FC stack rated current";
  parameter Real i_L_FC_stack(unit = "A") = 760 "FC stack cell maximum limiting current";
  parameter Real N_FC_stack(unit = "1") = 455 "FC stack number of cells";
  parameter Real A_FC_surf(unit = "m2") = 2*(L_FC_stack*W_FC_stack) + 2*(L_FC_stack*H_FC_stack) + 2*(W_FC_stack*H_FC_stack) "FC stack surface area";
  // Electrochemical parameters
  parameter Real i_0_FC_stack(unit = "A") = 0.0091 "FC stack cell exchange current";
  parameter Real i_x_FC_stack(unit = "A") = 0.001 "FC stack cell cross-over current";
  parameter Real b_1_FC_stack(unit = "V/dec") = 0.0985 "FC stack cell Tafel slope";
  parameter Real b_2_FC_stack(unit = "V/dec") = 0.0985 "FC stack cell trasport limitation factor";
  parameter Real R_O_FC_stack(unit = "Ohm") = 0.00022*N_FC_stack "FC stack cell ohmic resistance";
  // Thermal parameters
  parameter Real Cp_FC_stack(unit = "J/(kg.K)") = 110.0 "FC stack specific heat capacity";
  //*** DECLARE VARIABLES ***//
  // Physical constants
  import Modelica.Constants.R;
  import Modelica.Constants.F;
  // Fuel cell variables
  Modelica.Units.SI.Voltage V_cell;
  Modelica.Units.SI.Power P_th;
  Modelica.Units.SI.Pressure p_O2_log(min = 0);
  Modelica.Units.SI.Pressure p_0 = 100000;  
  Modelica.Units.SI.Pressure P_GDL_an_in[3];
  Modelica.Units.SI.Pressure P_GDL_an_out[3];
  Modelica.Units.SI.Pressure P_GDL_an_ave[3];
  Modelica.Units.SI.Pressure P_GDL_ca_in[3];
  Modelica.Units.SI.Pressure P_GDL_ca_out[3];
  Modelica.Units.SI.Pressure P_GDL_ca_ave[3];
  // mass flow
  Modelica.Units.SI.MassFlowRate H2_mass_flow;
  Modelica.Units.SI.MassFlowRate O2_mass_flow;
  //*** INSTANTIATE COMPONENTS ***//
  // Efficiencies
  Real eta_FC_LHV(unit = "100") "Lower heating value efficiency of fuel cell stack";
  // Electrical Components
  // Fluid Components
  Modelica.Fluid.Interfaces.FluidPort_b port_b_H2(redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-150, -100}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-110, -70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Interfaces.FluidPort_a port_a_H2(redeclare package Medium = Anode_Medium) annotation(
    Placement(visible = true, transformation(origin = {-150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a port_a_Coolant(redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-140, -126}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-30, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_b_Coolant(redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {140, -126}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {30, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Pipes.DynamicPipe pipeCoolant(redeclare package Medium = Coolant_Medium, T_start = 293.15, diameter = 0.003, length = 1, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = 1, nParallel = 500, p_a_start = 102502, use_HeatTransfer = true) annotation(
    Placement(transformation(origin = {0, -126}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  // Thermal Components
  // Other Components
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(
    Placement(visible = true, transformation(origin = {60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(
    Placement(visible = true, transformation(origin = {-60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Resistor R_ohmic(R = R_O_FC_stack) annotation(
    Placement(visible = true, transformation(origin = {60, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Electrical.Analog.Sources.SignalVoltage potentialSource annotation(
    Placement(visible = true, transformation(origin = {-60, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
    Placement(transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 10000) annotation(
    Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
    Placement(transformation(origin = {36, -60}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(visible = true, transformation(origin = {-100, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor annotation(
    Placement(transformation(origin = {-36, -42}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
  Modelica.Fluid.Interfaces.FluidPort_a port_a_Air(redeclare package Medium = Cathode_Medium) annotation(
    Placement(visible = true, transformation(origin = {150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b port_b_Air(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {150, -102}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Pipes.DynamicPipe channelAnode(nParallel = 500, length = 10, diameter = 0.003, redeclare model FlowModel = Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow, use_HeatTransfer = true, redeclare model HeatTransfer = Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer, redeclare package Medium = Anode_Medium, nNodes = 10, modelStructure = Modelica.Fluid.Types.ModelStructure.a_v_b) annotation(
    Placement(transformation(origin = {-120, -74}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Fluid.Pipes.DynamicPipe channelCathode(nParallel = 500, length = 1, diameter = 0.004, redeclare model FlowModel = Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow, use_HeatTransfer = true, redeclare model HeatTransfer = Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer, redeclare package Medium = Cathode_Medium, nNodes = 10, modelStructure = Modelica.Fluid.Types.ModelStructure.a_v_b) annotation(
    Placement(transformation(origin = {120, -74}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_H2(G = 50000) annotation(
    Placement(transformation(origin = {-92, -74}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_air(G = 50000) annotation(
    Placement(transformation(origin = {88, -74}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Units.SI.Power Power_stack;
  Modelica.Blocks.Math.Add water_prod annotation(
    Placement(transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Math.Gain gain(k = -1) annotation(
    Placement(transformation(origin = {42, 14}, extent = {{-10, -10}, {10, 10}})));
  Real mass_balance;
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = Cp_FC_stack*m_FC_stack, T(start = 333.15, fixed = true), der_T(fixed = false)) annotation(
    Placement(transformation(origin = {0, -39}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sensors.Temperature inlet_air_temperature(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {106, 104}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput OER_Cathode annotation(
    Placement(transformation(origin = {112, -114}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {108, -92}, extent = {{-10, -10}, {10, 10}})));
  //Modelica.Media.Water.StandardWater. WaterHumidityProps;
  Modelica.Blocks.Interfaces.RealOutput phi_Cathode annotation(
    Placement(transformation(origin = {112, -140}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {108, -112}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sensors.Temperature outlet_air_temperature(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {59, -91}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation(
    Placement(transformation(origin = {29, -34}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-1, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput OER_Anode annotation(
    Placement(transformation(origin = {-90, -122}, extent = {{10, -10}, {-10, 10}}), iconTransformation(origin = {-110, -90}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput phi_Anode annotation(
    Placement(transformation(origin = {-92, -144}, extent = {{10, -10}, {-10, 10}}), iconTransformation(origin = {-110, -110}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sensors.Temperature inlet_hydrogen_temperature(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {-120, 98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Units.SI.Current i_lim_react "Limiting current from reactive flows";
  Fluid.Membrane GDL_an(redeclare package Medium = Anode_Medium, V = 0.027, T_start = 353.15) annotation(
    Placement(transformation(origin = {-120, 14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Fluid.Membrane GDL_ca(redeclare package Medium = Cathode_Medium, V = 0.027, T_start = 353.15) annotation(
    Placement(transformation(origin = {120, 14}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain nitogen_cross_over(k = -1) annotation(
    Placement(transformation(origin = {1, -9}, extent = {{-10, -10}, {10, 10}})));
equation
  i_lim_react = min(port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_H2])/(Anode_Medium.MMX[i_H2]/(F*2)*N_FC_stack), port_a_Air.m_flow*inStream(port_a_Air.Xi_outflow[i_O2])/(Cathode_Medium.MMX[i_O2]/(F*4)*N_FC_stack));
  phi_Cathode = Cathode_Medium.relativeHumidity(Cathode_Medium.setState_phX(channelCathode.port_b.p, channelCathode.port_b.h_outflow, channelCathode.port_b.Xi_outflow));
  phi_Anode = Anode_Medium.relativeHumidity(Anode_Medium.setState_phX(channelAnode.port_b.p, channelAnode.port_b.h_outflow, channelAnode.port_b.Xi_outflow));
  OER_Cathode = port_a_Air.m_flow*inStream(port_a_Air.Xi_outflow[i_O2])/O2_mass_flow;
  OER_Anode = port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_H2])/H2_mass_flow;
//transfer in the membrane
  H2_mass_flow = Anode_Medium.MMX[i_H2]/(F*2)*N_FC_stack*min(i_lim_react, currentSensor.i);
  O2_mass_flow = Cathode_Medium.MMX[i_O2]/(F*4)*N_FC_stack*min(i_lim_react, currentSensor.i);
  GDL_an.mass_flow_a = {-H2_mass_flow, 0, 0};
  GDL_ca.mass_flow_a[i_O2] = -O2_mass_flow;
  GDL_an.mass_flow_b = {0, 0, 0};
  GDL_ca.mass_flow_b = {0, 0, 0};
  mass_balance = sum(GDL_an.mass_flow_a + GDL_ca.mass_flow_a);
  Power_stack = V_cell*N_FC_stack*pin_n.i;
//*** DEFINE EQUATIONS ***//
// activities
  P_GDL_an_in = Anode_Medium.partialPressure(actualStream(GDL_an.port_1.Xi_outflow), port_a_H2.p);
  P_GDL_an_out = Anode_Medium.partialPressure(actualStream(GDL_an.port_2.Xi_outflow), port_b_H2.p);
  P_GDL_an_ave = (P_GDL_an_in + P_GDL_an_ave)/2;
  P_GDL_ca_in = Cathode_Medium.partialPressure(actualStream(GDL_ca.port_1.Xi_outflow), port_a_Air.p);
  P_GDL_ca_out = Cathode_Medium.partialPressure(actualStream(GDL_ca.port_2.Xi_outflow), port_b_Air.p);
  P_GDL_ca_ave =  (P_GDL_ca_in + P_GDL_ca_out)/2;
//average activities of H2 and H2O based on logarithmic average of O2 partial pressure
  p_O2_log = (P_GDL_ca_in[i_O2] - P_GDL_ca_out[i_O2])/(log(P_GDL_ca_in[i_O2]) - log(P_GDL_ca_out[i_O2]));
// ELECTROCHEMICAL EQUATIONS //
// Calculate the stack voltage
  potentialSource.v = N_FC_stack*(1.229 + R*temperatureSensor.T/(2*F)*log((P_GDL_an_ave[i_H2]/p_0*(p_O2_log/p_0)^0.5)/(P_GDL_ca_ave[i_H2O]/p_0)) - b_1_FC_stack*log10((abs(currentSensor.i) + i_x_FC_stack)/i_0_FC_stack) + b_2_FC_stack*log10(1 - (abs(currentSensor.i) + i_x_FC_stack)/min(i_lim_react, i_L_FC_stack)));
// Calculate the voltage of the cell
  V_cell = pin_p.v/N_FC_stack;
// THERMAL EQUATIONS //
  P_th = (1.481 - V_cell)*abs(currentSensor.i)*N_FC_stack;
  GDL_an.Q_flow = if time < 1 then 100000 else 0;
  GDL_ca.Q_flow = if time < 1 then 100000 else 0;
// Assign the thermal power value to the heat flow component
  prescribedHeatFlow.Q_flow = P_th;
// Efficiencies
  eta_FC_LHV = (V_cell/1.253)*100;
//*** DEFINE CONNECTIONS ***//
  connect(pipeCoolant.port_b, port_b_Coolant) annotation(
    Line(points = {{10, -126}, {140, -126}}, color = {255, 0, 0}, thickness = 1));
  connect(pipeCoolant.port_a, port_a_Coolant) annotation(
    Line(points = {{-10, -126}, {-140, -126}}, color = {0, 0, 255}, thickness = 1));
  connect(R_ohmic.n, pin_p) annotation(
    Line(points = {{60, 130}, {60, 150}}, color = {0, 0, 255}));
  connect(pin_n, potentialSource.n) annotation(
    Line(points = {{-60, 150}, {-60, 150}, {-60, 130}, {-60, 130}}, color = {0, 0, 255}));
  connect(pin_n, ground.p) annotation(
    Line(points = {{-60, 150}, {-100, 150}, {-100, 138}, {-100, 138}}, color = {0, 0, 255}));
  connect(channelCathode.port_b, port_b_Air) annotation(
    Line(points = {{120, -84}, {120, -102}, {150, -102}}, color = {0, 127, 255}));
  connect(thermalConductor_air.port_b, channelCathode.heatPorts[1]) annotation(
    Line(points = {{98, -74}, {116, -74}}, color = {191, 0, 0}));
  connect(channelAnode.heatPorts[1], thermalConductor_H2.port_a) annotation(
    Line(points = {{-115.6, -74.1}, {-101.6, -74.1}}, color = {191, 0, 0}));
  connect(water_prod.y, gain.u) annotation(
    Line(points = {{0, 39}, {0, 14}, {30, 14}}, color = {0, 0, 127}));
  connect(pipeCoolant.heatPorts[1], thermalConductor.port_a) annotation(
    Line(points = {{0, -122}, {0, -110}}, color = {191, 0, 0}));
  connect(thermalConductor.port_b, heatCapacitor.port) annotation(
    Line(points = {{0, -90}, {0, -49}}, color = {191, 0, 0}));
  connect(thermalConductor_air.port_a, heatCapacitor.port) annotation(
    Line(points = {{78, -74}, {0, -74}, {0, -49}}, color = {191, 0, 0}));
  connect(prescribedHeatFlow.port, heatCapacitor.port) annotation(
    Line(points = {{26, -60}, {0, -60}, {0, -48}}, color = {191, 0, 0}));
  connect(temperatureSensor.port, heatCapacitor.port) annotation(
    Line(points = {{-36, -52}, {-36, -60}, {0, -60}, {0, -48}}, color = {191, 0, 0}));
  connect(port_a_Air, inlet_air_temperature.port) annotation(
    Line(points = {{150, 80}, {106, 80}, {106, 94}}));
  connect(outlet_air_temperature.port, port_b_Air) annotation(
    Line(points = {{59, -101}, {150, -101}, {150, -102}}, color = {0, 127, 255}));
  connect(thermalConductor_H2.port_b, heatCapacitor.port) annotation(
    Line(points = {{-82, -74}, {0, -74}, {0, -48}}, color = {191, 0, 0}));
  connect(heatCapacitor.port, heatPort) annotation(
    Line(points = {{0, -48}, {30, -48}, {30, -34}}, color = {191, 0, 0}));
  connect(channelAnode.port_b, port_b_H2) annotation(
    Line(points = {{-120, -84}, {-120, -100}, {-150, -100}}, color = {0, 127, 255}));
  connect(port_a_H2, inlet_hydrogen_temperature.port) annotation(
    Line(points = {{-150, 80}, {-120, 80}, {-120, 88}}));
  connect(port_a_H2, GDL_an.port_1) annotation(
    Line(points = {{-150, 80}, {-120, 80}, {-120, 24}}));
  connect(GDL_an.port_2, channelAnode.port_a) annotation(
    Line(points = {{-120, 4}, {-120, -64}}, color = {0, 127, 255}));
  connect(port_a_Air, GDL_ca.port_1) annotation(
    Line(points = {{150, 80}, {120, 80}, {120, 24}}));
  connect(GDL_ca.port_2, channelCathode.port_a) annotation(
    Line(points = {{120, 4}, {120, -64}}, color = {0, 127, 255}));
  connect(GDL_an.mass_flow_a[i_H2], water_prod.u2) annotation(
    Line(points = {{-116, 14}, {-72, 14}, {-72, 78}, {-6, 78}, {-6, 62}}, color = {0, 0, 127}));
  connect(GDL_ca.mass_flow_a[i_O2], water_prod.u1) annotation(
    Line(points = {{116, 14}, {80, 14}, {80, 76}, {6, 76}, {6, 62}}, color = {0, 0, 127}));
  connect(gain.y, GDL_ca.mass_flow_a[i_H2O]) annotation(
    Line(points = {{54, 14}, {116, 14}}, color = {0, 0, 127}));
  connect(GDL_an.mass_flow_a[i_N2], nitogen_cross_over.u) annotation(
    Line(points = {{-116, 14}, {-72, 14}, {-72, -8}, {-10, -8}}, color = {0, 0, 127}));
  connect(nitogen_cross_over.y, GDL_ca.mass_flow_a[i_N2]) annotation(
    Line(points = {{12, -8}, {80, -8}, {80, 14}, {116, 14}}, color = {0, 0, 127}));
  connect(potentialSource.p, currentSensor.p) annotation(
    Line(points = {{-60, 110}, {-60, 100}, {-10, 100}}, color = {0, 0, 255}));
  connect(currentSensor.n, R_ohmic.p) annotation(
    Line(points = {{10, 100}, {60, 100}, {60, 110}}, color = {0, 0, 255}));
  annotation(
    Diagram(coordinateSystem(extent = {{-150, -150}, {150, 150}}, initialScale = 0.1)),
    Icon(coordinateSystem(extent = {{-150, -150}, {150, 150}}, initialScale = 0.1), graphics = {Line(origin = {20.1754, 1.92106}, points = {{0, 78}, {0, -80}, {0, -82}}), Rectangle(origin = {80, 0}, fillColor = {0, 178, 227}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-20, 100}, {20, -100}}), Line(origin = {40.1315, 2}, points = {{0, 78}, {0, -80}, {0, -82}}), Line(origin = {0.219199, 1.92106}, points = {{0, 78}, {0, -80}, {0, -82}}), Line(origin = {-40.0001, 1.61404}, points = {{0, 78}, {0, -80}, {0, -82}}), Rectangle(origin = {-80, 0}, fillColor = {170, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-20, 100}, {20, -100}}), Text(origin = {10, -54}, textColor = {255, 0, 0}, extent = {{-11, 6}, {11, -6}}, textString = "K"), Line(origin = {-20.0439, -0.307018}, points = {{0, 80}, {0, -80}, {0, -80}}), Rectangle(origin = {35, 54}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Vertical, extent = {{-95, 26}, {25, -134}}), Text(origin = {-80, 6}, extent = {{-26, 24}, {26, -24}}, textString = "A"), Text(origin = {80, 6}, extent = {{-26, 24}, {26, -24}}, textString = "C"), Text(origin = {74, -94}, extent = {{-22, 10}, {22, -10}}, textString = "OER"), Text(origin = {75, -113}, extent = {{-21, 11}, {21, -11}}, textString = "phi"), Text(origin = {-86, -94}, extent = {{-22, 10}, {22, -10}}, textString = "OER"), Text(origin = {-85, -113}, extent = {{-21, 11}, {21, -11}}, textString = "phi")}),
    Documentation(info = "<html><head></head><body>This model describes the dynamic behaviour of a proton exchange membrane fuel cell (PEMFC) stack. The model includes components describing the electrical, fluidic, and thermal properties of the stack.&nbsp;<div><br></div><div>The electrical performance is modelled using a 0-D polarization curve model , which incorporates Nernstian thermodynamic effects due to hydrogen and oxygen pressure changes, Tafel kinetics to calculate activation overpotentials, and an empirical relationship to calculate mass-transport overpotentials. These effects are combined in&nbsp;<span style=\"font-family: 'Courier New';\">potentialSource.v</span><span style=\"font-family: 'Courier New'; font-size: 12pt;\">,</span>which calculates the open-circuit voltage for a single cell, adjusts for hydrogen and oxygen partial pressures, subtracts the activation and mass-transport overpotentials, and finally multiplies by the number of cells in the stack. A simple resistor is included after the potential source to cover all Ohmic resistive losses in the fuel cell. Default parameters fit the polarization curve given by Powercell in their Powercellution data sheet, available <a href=\"https://powercellution.com/p-stack\">here</a>.</div><div><br></div><div>The fluidic performance is modelled using simple ideal flow components for the air and hydrogen gas lines, connected to mass sink boundary conditions. The magnitude of the mass sink is coupled to the electrical current in the stack using Faraday's law.&nbsp;&nbsp;
</div><div><br></div><div>The thermal performance is considered by coupling a model describing the flow of liquid coolant to a thermal heat source. The magnitude of the heat source is calculated using the higher heating value of hydrogen and the calculated electrical voltage of the cell.<div><br></div><div>The hydrogen, air, and coolant ports can be connected to their respective subsystems, either by using the <a href=\"modelica://VirtualFCS.SubSystems.FuelCellSubSystems\">FuelCellSubSystems</a> block, or individual <a href=\"modelica://VirtualFCS.SubSystems.Hydrogen.SubSystemHydrogen\">SubSystemHydrogen</a>, <a href=\"modelica://VirtualFCS.SubSystems.Air.SubSystemAir\">SubSystemAir</a>, and <a href=\"modelica://VirtualFCS.SubSystems.Cooling.SubSystemCooling\">SubSystemCooling</a> blocks.<br>&nbsp; 

<table border=\"0.9\"><caption style=\"text-align: left;\" align=\"Left\"><strong><u>Default Parameters</u></strong></caption><caption style=\"text-align: left;\" align=\"Left\"><strong><u><br></u></strong></caption>
<tbody>
<tr>
<th>Parameter name</th>
<th>Value</th>
<th>Unit</th>
</tr>
<tr>
<td align=\"Left\">m_FC_stack</td>
<td>=42</td>
<td align=\"Right\">kg</td>
</tr>
<tr>
<td align=\"Left\">L_FC_stack</td>
<td>=0.42</td>
<td align=\"Right\">m</td>
</tr>
<tr>
<td align=\"Left\">W_FC_stack</td>
<td>=0.582</td>
<td align=\"Right\">m</td>
</tr>
<tr>
<td align=\"Left\">H_FC_stack</td>
<td>=0.156</td>
<td align=\"Right\">m</td>
</tr>
<tr>
<td align=\"Left\">I_rated_FC_stack</td>
<td>=450</td>
<td align=\"Right\">A</td>
</tr>
<tr>
<td align=\"Left\">i_L_FC_stack</td>
<td>=760</td>
<td align=\"Right\">A</td>
</tr>

<tr>
<td align=\"Left\">N_FC_stack</td>
<td>=455</td>
<td align=\"Right\">-</td>
</tr>

<tr>
<td align=\"Left\">i_0_FC_stack</td>
<td>=0.0091</td>
<td align=\"Right\">A</td>
</tr>
<tr>
<td align=\"Left\">i_x_FC_stack</td>
<td>=0.001</td>
<td align=\"Right\">A</td>
</tr>
<tr>
<td align=\"Left\">b_1_FC_stack</td>
<td>=0.0985</td>
<td align=\"Right\">V/dec</td>
</tr>
<tr>
<td align=\"Left\">b_2_FC_stack</td>
<td>=0.0985</td>
<td align=\"Right\">V/dec</td>
</tr>
<tr>
<td align=\"Left\">R_0</td>
<td>=0.00022*N_FC_stack</td>
<td align=\"Right\">Ohm</td>
</tr>
<tr>
<td align=\"Left\">Cp</td>
<td>=1100</td>
<td align=\"Left\">J/(kg.K)</td>
</tr>

</tbody>
</table><br><br><br>



<div><span style=\"text-decoration: underline;\"><strong>Electrochemical equations: </strong></span></div><div>In the equations below, i<sub>stack</sub>&nbsp;represents the current flowing through the stack, accessible in the code as <font face=\"Courier New\">currentSensor.i</font>.</div>
<p><i><u>The Nernst equilibrium potential, per cell</u>&nbsp;</i></p>
<p>U<sub>FC</sub><sup>Nernst </sup>= (U<sup>0</sup> -((RT)/(2F) ln( 1/(p<sub>H2</sub> (p<sub>O2</sub><sup>0.5</sup>))), U<sup>0 </sup>= 1.229 V</p>
<p><span style=\"text-decoration: underline;\"><i>Activation overpotential, per cell</i></span></p>
<p>η<sup>act </sup>= b<sub>1 </sub>ln( 1-(i<sub>stack&nbsp;</sub>+ i<sub>x</sub>) / i<sub>0</sub>)</p>
<p><u><i>Concentration overpotential, per cell</i></u></p>
<p>η<sup>con </sup>= -b<sub>2 </sub>ln( 1-(i<sub>stack&nbsp;</sub>+ i<sub>x</sub>) / i<sub>L</sub>)</p><p><u><i>Stack voltage</i></u></p><p>V<sub>stack</sub> = N<sub>cell</sub> (U<sub>FC</sub><sup>Nernst</sup>&nbsp;- η<sup>act&nbsp;</sup>&nbsp;- i<sub>FC</sub>R<sub>0</sub> - η<sup>con</sup>)</p>
<p><span style=\"text-decoration: underline;\"><strong>Thermal equations:</strong> </span></p>
<p><i><u>Electrochemical heat generation</u></i></p>
<p>Q<sub>gen</sub><sup>&nbsp;</sup>= (V<sub>TN</sub> - V<sub><font size=\"2\">stack</font></sub>)i<sub>stack</sub>, V<sub>TN</sub> = 1.481 V</p>
<p><br></p>
<p>&nbsp;</p>


</div></div></body></html>"),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
end FuelCellStack;
