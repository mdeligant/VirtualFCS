within VirtualFCS.Electrochemical.Hydrogen;

model FuelCellStack
  // sources
  // for water transport in the membrane
  //[1] K. Jiao et X. Li, « Water transport in polymer electrolyte membrane fuel cells », Progress in Energy and Combustion Science, vol. 37, nᵒ 28, p. 221‑291, juin 2011, doi: 10.1016/j.pecs.2010.06.002.
  // for calcul of potential
  //[2] B. Xie, G. Zhang, J. Xuan, et K. Jiao, « Three-dimensional multi-phase model of PEM fuel cell coupled with improved agglomerate sub-model of catalyst layer », Energy Conversion and Management, vol. 199, p. 112051, nov. 2019, doi: 10.1016/j.enconman.2019.112051.
  // for diffusion in GDL
  // [3] B. Xie et al., « Validation methodology for PEM fuel cell three-dimensional simulation », International Journal of Heat and Mass Transfer, vol. 189, nᵒ 25, p. 122705, juin 2022, doi: 10.1016/j.ijheatmasstransfer.2022.122705.
  // for diffusion of H2 and N2 trought membrane
  //[4] R. Omrani et B. Shabani, « An analytical model for hydrogen and nitrogen crossover rates in proton exchange membrane fuel cells », International Journal of Hydrogen Energy, vol. 45, nᵒ 55, p. 31041‑31055, nov. 2020, doi: 10.1016/j.ijhydene.2020.08.089.
  // for heat exchange
  //[5] D. Musse et D. Lee, « Computational evaluation of PEMFC performance based on bipolar plate material types », Energy Reports, vol. 11, nᵒ 18, p. 4886‑4903, juin 2024, doi: 10.1016/j.egyr.2024.04.052.
  // for dimention off cell and bipolar plate
  // [6] [1] P. von Tettau et al., « Laboratory assessments applied to mass-produced automotive fuel cells », International Journal of Hydrogen Energy, vol. 52, p. 1127‑1136, janv. 2024, doi: 10.1016/j.ijhydene.2023.10.309.
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
  parameter Integer N_el = 10 "discretization of heat exchange";
  //*** DECLARE PARAMETERS ***//
  // Physical parameters
  // Fuel Cell Stack Paramters
  parameter Real m_FC_stack(unit = "kg") = 42 "FC stack mass";
  parameter Real L_FC_stack(unit = "m") = 0.420 "FC stack length";
  parameter Real W_FC_stack(unit = "m") = 0.582 "FC stack width";
  parameter Real H_FC_stack(unit = "m") = 0.156 "FC stack height";
  parameter Real vol_FC_stack(unit = "m3") = L_FC_stack*W_FC_stack*H_FC_stack "FC stack volume";
  parameter Real N_FC_stack(unit = "1") = 370 "FC stack number of cells";
  parameter Real A_FC_surf(unit = "m2") = 2*(L_FC_stack*W_FC_stack) + 2*(L_FC_stack*H_FC_stack) + 2*(W_FC_stack*H_FC_stack) "FC stack surface area";
  // Thermal parameters
  parameter Real Cp_FC_stack(unit = "J/(kg.K)") = 800 "FC stack specific heat capacity";
  // cell parameter
  parameter Modelica.Units.SI.Density rhoMem = 1980;
  //density of the memebrane
  parameter Modelica.Units.SI.MolarMass EW = 1.100;
  // molar mass of the membrane
  parameter Real GDL_porosity = 0.55;
  //[2]
  parameter Modelica.Units.SI.Length t_CL_an = 3*10^(-6);
  //thickness of catalist layer in anode [2]
  parameter Modelica.Units.SI.Length t_CL_ca = 10*10^(-6);
  //thickness of catalist layer in cathode [2]
  parameter Modelica.Units.SI.Length t_GDL_an = 240*10^(-6);
  //thickness of gaz diffusion layer in anode [6]
  parameter Modelica.Units.SI.Length t_GDL_ca = 240*10^(-6);
  //thickness of gaz diffusion layer in cathode [6]
  parameter Modelica.Units.SI.Length t_mem = 20*10^(-6);
  // thickness of memebrane [6]
  parameter Modelica.Units.SI.Length t_bp = 0.13/100;
  // thickness of bipolar plate betwwen channel and coollant (toyota review 66)
  parameter Real number_channel = 72 "number of channel by cell [6]";
  parameter Modelica.Units.SI.Length channel_length = Cell_Area/(w_channel_an + w_land_an)/number_channel;
  parameter Modelica.Units.SI.Area channel_area_an = 0.121/10^6;
  // area of channel anode [6]
  parameter Modelica.Units.SI.Length channel_perimeter_an = 1.5/1000;
  // perimeter of channel anode [6]
  parameter Modelica.Units.SI.Area channel_area_ca = 0.204/10^6;
  // area of channel cathode [6]
  parameter Modelica.Units.SI.Length channel_perimeter_ca = 1.6/1000;
  // perimeter of channel cathode [6]
  parameter Modelica.Units.SI.Length w_channel_an = 0.56/1000;
  // witdh of channel [6]
  parameter Modelica.Units.SI.Length w_land_an = 0.84/1000;
  // width of land part next to channel [6]
  parameter Modelica.Units.SI.Length w_channel_ca = 0.61/1000;
  // witdh of channel [6]
  parameter Modelica.Units.SI.Length w_land_ca = 0.79/1000;
  // width of land part next to channel [6]
  //thickness of the membrane
  parameter Modelica.Units.SI.Area Cell_Area = 0.0237;
  // active area of a single cell
  // activation overvoltage parameter [2]
  parameter Modelica.Units.SI.CurrentDensity i_0_ref_an = 3.5;
  //Reference exchange current density
  parameter Modelica.Units.SI.CurrentDensity i_0_ref_ca = 3.5*10^(-4);
  //Reference exchange current density
  parameter Real alpha_an(unit = "1") = 0.5;
  //Transfer coefficient
  parameter Real alpha_ca(unit = "1") = 0.5;
  //Transfer coefficient
  parameter Real m_pt_an(unit = "kg/m2") = 0.05/10^6*10000;
  //Pt loading [2]
  parameter Real m_pt_ca(unit = "kg/m2") = 0.05/10^6*10000;
  //Pt loading [2]
  parameter Real ECSA(unit = "m2/kg") = 70*1000;
  //electrochemically active surface area [2]
  parameter Real r_pt_c(unit = "1") = m_pt_ca/(0.05/10^6*10000)*0.1;
  //platinum weight percentage of Pt/carbon catalyst
  parameter Real r_im_c(unit = "1") = 0.95;
  //mass ratio of ionomer to carbon
  parameter Modelica.Units.SI.MolalConcentration C_H2_ref = 56.4;
  // reference concentration of H2 in the ionomer
  parameter Modelica.Units.SI.MolalConcentration C_O2_ref = 3.39;
  // reference concentration of O2 in the ionomer
  parameter Modelica.Units.SI.MolarEnergy w = 10000;
  // energy parameter link to theta_Pt_o_cov
  // parameter for heat exchange
  // heat capacity and exchange [5]
  parameter Modelica.Units.SI.ThermalConductivity T_conductivity_GDL = 1.6 "ThermalConductivity of gas difusion layer";
  parameter Modelica.Units.SI.ThermalConductivity T_conductivity_CL = 8 "ThermalConductivity of catalyst layer";
  parameter Modelica.Units.SI.ThermalConductivity T_conductivity_mem = 0.67 "ThermalConductivity of PEM";
  parameter Modelica.Units.SI.ThermalConductivity T_conductivity_bp = 21.9 "ThermalConductivity of bipolar plate";
  // Conductivity of titanium [6]
  parameter Modelica.Units.SI.ThermalConductance Thermal_conductance_mem = T_conductivity_mem/t_mem*Cell_Area*N_FC_stack;
  parameter Modelica.Units.SI.ThermalConductance Thermal_conductance_an = 1/(t_CL_an/T_conductivity_CL + t_GDL_an/T_conductivity_GDL)*Cell_Area*N_FC_stack;
  parameter Modelica.Units.SI.ThermalConductance Thermal_conductance_ca = 1/(t_CL_ca/T_conductivity_CL + t_GDL_ca/T_conductivity_GDL)*Cell_Area*N_FC_stack;
  parameter Modelica.Units.SI.ThermalConductance Thermal_conductance_coollant = T_conductivity_bp/t_bp*Cell_Area*N_FC_stack*(w_land_an/(w_channel_an + w_land_an) + w_land_ca/(w_channel_ca + w_land_ca));
  //thermic transfer with only the bp consider and no contribution of gas in the channel
  //*** DECLARE VARIABLES ***//
  // Physical constants
  import Modelica.Constants.R;
  import Modelica.Constants.F;
  // Fuel cell variables
  Modelica.Units.SI.Voltage V_cell;
  Modelica.Units.SI.CurrentDensity j_cell;
  Modelica.Units.SI.Voltage V_rev;
  Modelica.Units.SI.Voltage V_nernst;
  Modelica.Units.SI.Voltage V_act_an;
  Modelica.Units.SI.Voltage V_act_ca;
  Modelica.Units.SI.Voltage V_act;
  Modelica.Units.SI.Voltage V_ohm;
  Modelica.Units.SI.Power P_th;
  Modelica.Units.SI.Power P_in;
  Modelica.Units.SI.Power P_out;
  // partial pressure
  Modelica.Units.SI.Pressure p_O2_log(min = 0);
  Modelica.Units.SI.Pressure p_0 = 100000;
  Modelica.Units.SI.Pressure P_CL_an_in[3];
  Modelica.Units.SI.Pressure P_CL_an_out[3];
  Modelica.Units.SI.Pressure P_CL_an_ave[3];
  Modelica.Units.SI.Pressure P_CL_ca_in[3];
  Modelica.Units.SI.Pressure P_CL_ca_out[3];
  Modelica.Units.SI.Pressure P_CL_ca_ave[3];
  Modelica.Units.SI.Pressure P_GDL_an_in[3];
  Modelica.Units.SI.Pressure P_GDL_an_out[3];
  Modelica.Units.SI.Pressure P_GDL_an_ave[3];
  Modelica.Units.SI.Pressure P_GDL_ca_in[3];
  Modelica.Units.SI.Pressure P_GDL_ca_out[3];
  Modelica.Units.SI.Pressure P_GDL_ca_ave[3];
  // temperature
  Modelica.Units.SI.Temperature T_an;
  Modelica.Units.SI.Temperature T_ca;
  Modelica.Units.SI.Temperature T_ave;
  // mass flow due to electrochemical reaction
  Modelica.Units.SI.MassFlowRate H2_mass_flow;
  Modelica.Units.SI.MassFlowRate O2_mass_flow;
  // water mass flow trought membrane based on [1]
  Modelica.Units.SI.MassFlowRate EOD_mass_flow;
  // water mass flow due to electro osmoic drag
  Real Nd;
  //Characterize the number of water molecules transported through the membrane for hydrogen transport
  Modelica.Units.SI.DiffusionArea Dnmw;
  // diffusivity of water in the membrane
  Modelica.Units.SI.MassFlowRate diff_mass_flow;
  // water mass flow due to diffusion
  Modelica.Units.SI.MassFlowRate water_mass_flow_billan;
  Modelica.Units.SI.MassFlowRate water_mass_flow_billan_theo;
  // the water mass flow is limited to avoid convergence issue during initialisation. after initialisation water_mass_flow_billan and  water_mass_flow_billan_theo should be always equal
  // diffusion of N2 and H2 trought membrane [4]
  Modelica.Units.SI.MassFlowRate H2_mem_crossover;
  Modelica.Units.SI.MassFlowRate N2_mem_crossover;
  Modelica.Units.SI.MassFlowRate O2_mem_crossover_reac;
  Real mem_N2_diff(unit = "mol/m/s/Pa");
  Real mem_N2_diff_0(unit = "mol/m/s/Pa");
  Modelica.Units.SI.MolarEnergy mem_N2_diff_Ea;
  Real mem_H2_diff(unit = "mol/m/s/Pa");
  Real mem_H2_diff_0(unit = "mol/m/s/Pa");
  Modelica.Units.SI.MolarEnergy mem_H2_diff_Ea;
  Modelica.Units.SI.CurrentDensity j_H2_crossover;
  //humidity in the cell
  Real phi_Anode_v[N_el];
  Real phi_Cathode_v[N_el];
  Real lambda_an(unit = "1");
  // water content of the membrane anode side
  Real lambda_ca(unit = "1");
  // water content of the membrane cathode side
  Real lambda_m(unit = "1");
  // water content of the membrane average
  // diffusion in GDL
  Modelica.Units.SI.DiffusionArea D_an[3];
  Modelica.Units.SI.DiffusionArea D_ca[3];
  // ohmic overvoltage variable [2]
  Modelica.Units.SI.Conductivity Conductivity;
  // dependant of water content of the membrane
  // activation overvoltage variable [2]
  Real A_eff_pt_an(unit = "/m");
  //specific surface area of pt anode side
  Real A_eff_pt_ca(unit = "/m");
  //specific surface area of pt cathode side
  Real theta_T_an(unit = "1");
  // temperature correction coefficient anode side
  Real theta_T_ca(unit = "1");
  // temperature correction coefficient cathode side
  Real theta_Pt_o_cov(unit = "1");
  // Pt-oxide coverage-dependent correction coefficient
  Real theta_Pt_im_cov(unit = "1");
  // Pt-ionomer coverage-dependent correction coefficient
  Real H_H2(unit = "Pa.m3/mol");
  //henry constant to have the concentration in the ionomer next to platinium
  Real H_O2(unit = "Pa.m3/mol");
  //henry constant to have the concentration in the ionomer next to platinium
  //*** INSTANTIATE COMPONENTS ***//
  // Efficiencies
  Real eta_FC_LHV(unit = "100") "Lower heating value efficiency of fuel cell stack";
   Modelica.Units.SI.MassFlowRate mass_balance[4];
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
  Modelica.Fluid.Pipes.DynamicPipe pipeCoolant(redeclare package Medium = Coolant_Medium, T_start = 333.15, diameter = 0.0004, length = channel_length, modelStructure = Modelica.Fluid.Types.ModelStructure.a_vb, nNodes = N_el, nParallel = N_FC_stack*number_channel, p_a_start = 102502, use_HeatTransfer = true) annotation(
    Placement(transformation(origin = {0, -126}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  // Thermal Components
  // Other Components
  Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation(
    Placement(visible = true, transformation(origin = {60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation(
    Placement(visible = true, transformation(origin = {-60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-60, 150}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sources.SignalVoltage potentialSource annotation(
    Placement(visible = true, transformation(origin = {-60, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Electrical.Analog.Sensors.CurrentSensor currentSensor annotation(
    Placement(transformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor[N_el](each G = Thermal_conductance_coollant/N_el) annotation(
    Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow_an[N_el] annotation(
    Placement(transformation(origin = {-84, -96}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
    Placement(visible = true, transformation(origin = {-100, 128}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a port_a_Air(redeclare package Medium = Cathode_Medium) annotation(
    Placement(visible = true, transformation(origin = {150, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {110, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b port_b_Air(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {150, -102}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Pipes.DynamicPipe channelAnode(nParallel = N_FC_stack*number_channel, length = channel_length, diameter = 0.0002, redeclare model FlowModel = Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow, use_HeatTransfer = true, redeclare model HeatTransfer = Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer, redeclare package Medium = Anode_Medium, nNodes = N_el, modelStructure = Modelica.Fluid.Types.ModelStructure.a_v_b, T_start = 333.15, p_a_start = 105000, isCircular = false, crossArea = channel_area_an, perimeter = channel_perimeter_an) annotation(
    Placement(transformation(origin = {-120, -74}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Fluid.Pipes.DynamicPipe channelCathode(nParallel = N_FC_stack*number_channel, length = channel_length, diameter = 0.0002, redeclare model FlowModel = Modelica.Fluid.Pipes.BaseClasses.FlowModels.DetailedPipeFlow, use_HeatTransfer = true, redeclare model HeatTransfer = Modelica.Fluid.Pipes.BaseClasses.HeatTransfer.LocalPipeFlowHeatTransfer, redeclare package Medium = Cathode_Medium, nNodes = N_el, modelStructure = Modelica.Fluid.Types.ModelStructure.a_v_b, T_start = 333.15, p_a_start = 1.2e5, isCircular = false, crossArea = channel_area_ca, perimeter = channel_perimeter_ca) annotation(
    Placement(transformation(origin = {130, -74}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_an[N_el](each G = Thermal_conductance_an/N_el) annotation(
    Placement(transformation(origin = {-50, -74}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_ca[N_el](each G = Thermal_conductance_ca/N_el) annotation(
    Placement(transformation(origin = {50, -74}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Units.SI.Power Power_stack;
  Modelica.Blocks.Math.Add water_prod annotation(
    Placement(transformation(origin = {2, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));

  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor[N_el](each C = Cp_FC_stack*m_FC_stack/N_el, each T(start = 333.15, fixed = true), each der_T(fixed = false)) annotation(
    Placement(transformation(origin = {0, -61}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sensors.Temperature inlet_air_temperature(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {106, 104}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput OER_Cathode annotation(
    Placement(transformation(origin = {116, -140}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 108}, extent = {{-10, -10}, {10, 10}})));
  //Modelica.Media.Water.StandardWater. WaterHumidityProps;
  Modelica.Blocks.Interfaces.RealOutput phi_Cathode annotation(
    Placement(transformation(origin = {80, -140}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {108, -112}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput purge_valve_indication annotation(
    Placement(transformation(origin = {-54, -154}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-84, -136}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput Phi_cathode_out annotation(
    Placement(transformation(origin = {54, -146}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {76, -138}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
 
  Modelica.Fluid.Sensors.Temperature outlet_air_temperature(redeclare package Medium = Cathode_Medium) annotation(
    Placement(transformation(origin = {205, -91}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort[N_el] annotation(
    Placement(transformation(origin = {47, -50}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-1, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput OER_Anode annotation(
    Placement(transformation(origin = {-84, -142}, extent = {{10, -10}, {-10, 10}}), iconTransformation(origin = {-110, -90}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput phi_Anode annotation(
    Placement(transformation(origin = {-122, -142}, extent = {{10, -10}, {-10, 10}}), iconTransformation(origin = {-110, -110}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sensors.Temperature inlet_hydrogen_temperature(redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-120, 98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Units.SI.Current i_lim_react "Limiting current from reactive flows";
  Fluid.Membrane CL_an(redeclare package Medium = Anode_Medium, V = N_FC_stack*Cell_Area*(t_GDL_an/2 + t_CL_an)*GDL_porosity, T_start = 353.15, X_start = {0.5, 0.4, 0.1}, p_start = 1.1e5) annotation(
    Placement(transformation(origin = {-60, 14}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Fluid.Membrane CL_ca(redeclare package Medium = Cathode_Medium, V = N_FC_stack*Cell_Area*(t_GDL_ca/2 + t_CL_an)*GDL_porosity, T_start = 353.15, p_start = 1.3e5) annotation(
    Placement(transformation(origin = {80, 14}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain nitogen_cross_over(k = -1) annotation(
    Placement(transformation(origin = {3, -7}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add water_prod_and_transfer(k1 = -1, k2 = -1) annotation(
    Placement(transformation(origin = {36, 14}, extent = {{-10, -10}, {10, 10}})));
  VirtualFCS.Fluid.Membrane GDL_an(redeclare package Medium = Anode_Medium, T_start = 353.15, V = N_FC_stack*number_channel*channel_length*channel_area_an + N_FC_stack*Cell_Area*t_GDL_an/2*GDL_porosity, X_start = {0.5, 0.4, 0.1}, p_start = 1.1e5) annotation(
    Placement(transformation(origin = {-120, 13}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  VirtualFCS.Fluid.Membrane GDL_ca(redeclare package Medium = Cathode_Medium, T_start = 353.15, V = N_FC_stack*number_channel*channel_length*channel_area_ca + N_FC_stack*Cell_Area*t_GDL_ca/2*GDL_porosity, p_start = 1.3e5) annotation(
    Placement(transformation(origin = {130, 14}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Gain diff_trought_GDL_an[3](each k = -1) annotation(
    Placement(transformation(origin = {-89, 14}, extent = {{-6, -6}, {6, 6}})));
  Modelica.Blocks.Math.Gain diff_trought_GDL_ca[3](each k = -1) annotation(
    Placement(transformation(origin = {104, 14}, extent = {{6, -6}, {-6, 6}}, rotation = -0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow_ca[N_el] annotation(
    Placement(transformation(origin = {92, -102}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_mem[N_el](each G = Thermal_conductance_mem/N_el) annotation(
    Placement(transformation(origin = {0, -36}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sensors.Temperature outlet_hydrogen_temperature(redeclare package Medium = Anode_Medium) annotation(
    Placement(transformation(origin = {-194, -89}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Fluid.Sensors.Temperature inlet_coolant_temperature(redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {-170, -141}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Fluid.Sensors.Temperature outlet_coolant_temperature(redeclare package Medium = Coolant_Medium) annotation(
    Placement(transformation(origin = {179, -135}, extent = {{10, -10}, {-10, 10}})));
equation
  i_lim_react = min(port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_H2])/(Anode_Medium.MMX[i_H2]/(F*2)*N_FC_stack), port_a_Air.m_flow*inStream(port_a_Air.Xi_outflow[i_O2])/(Cathode_Medium.MMX[i_O2]/(F*4)*N_FC_stack));
  j_cell*Cell_Area = currentSensor.i;
// temperature
  T_an = sum(channelAnode.heatPorts.T)/N_el;
  T_ca = sum(channelCathode.heatPorts.T)/N_el;
  T_ave = (T_an + T_ca)/2;
// humidity
  for i in 1:N_el loop
     phi_Anode_v[i] = ((P_CL_an_out[i_H2O] - P_CL_an_in[i_H2O])/(N_el-1)*(i-1)+ P_CL_an_in[i_H2O])/Modelica.Media.Water.StandardWater.saturationPressure(channelAnode.heatPorts[i].T);
  phi_Cathode_v[i] = ((P_CL_ca_out[i_H2O] - P_CL_ca_in[i_H2O])/(N_el-1)*(i-1)+ P_CL_ca_in[i_H2O])/Modelica.Media.Water.StandardWater.saturationPressure(channelCathode.heatPorts[i].T);
    end for;
  phi_Cathode = Cathode_Medium.relativeHumidity(Cathode_Medium.setState_pTX(CL_ca.medium.p, T_ca, CL_ca.medium.X));
  phi_Anode = Anode_Medium.relativeHumidity(Anode_Medium.setState_pTX(CL_an.medium.p, T_an, CL_an.medium.X));
  //Drain_valve_indication = max(max(phi_Anode_v[end], phi_Cathode_v[end]), port_b_H2.Xi_outflow[i_N2]*30);
  purge_valve_indication = max(phi_Anode_v[end], P_GDL_an_ave[i_N2]/port_a_H2.p*2);
  Phi_cathode_out = phi_Cathode_v[end];
  lambda_an = if phi_Anode < 1 then 0.043 + 17.81*phi_Anode - 39.85*phi_Anode^2 + 36*phi_Anode^3 else 14 + 1.4*(phi_Anode - 1)*2*1.27/1000;
// could be modified to have an expression dependant of temperature
  lambda_ca = if phi_Cathode < 1 then 0.043 + 17.81*phi_Cathode - 39.85*phi_Cathode^2 + 36*phi_Cathode^3 else 14 + 1.4*(phi_Cathode - 1)*2*1.27/1000;
  lambda_m = (lambda_an + lambda_ca)/2;
// hydrogen and oxygen excess ratio
  OER_Cathode = smooth(1, if (O2_mass_flow < -10^(-4) or O2_mass_flow > 10^(-4)) then channelCathode.m_flows[1]*inStream(port_a_Air.Xi_outflow[i_O2])/O2_mass_flow else 0);
  H2_mass_flow*OER_Anode = port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_H2]);
//transfer in the membrane
  H2_mass_flow = Anode_Medium.MMX[i_H2]/(F*2)*N_FC_stack*currentSensor.i;
  O2_mass_flow = Cathode_Medium.MMX[i_O2]/(F*4)*N_FC_stack*min(i_lim_react, currentSensor.i);
  H2_mem_crossover = mem_H2_diff/t_mem*P_CL_an_ave[i_H2]*Cell_Area*N_FC_stack*Anode_Medium.MMX[i_H2];
  N2_mem_crossover = mem_N2_diff/t_mem*(P_CL_ca_ave[i_N2] - P_CL_an_ave[i_N2])*Cell_Area*N_FC_stack*Anode_Medium.MMX[i_N2];
  O2_mem_crossover_reac = mem_H2_diff/t_mem*P_CL_an_ave[i_H2]*Cell_Area*N_FC_stack*Anode_Medium.MMX[i_O2]/2;
  j_H2_crossover = H2_mem_crossover*F*2/(Anode_Medium.MMX[i_H2]*N_FC_stack*Cell_Area);
  CL_an.mass_flow_a = {-H2_mass_flow - H2_mem_crossover, water_mass_flow_billan, N2_mem_crossover};
  CL_ca.mass_flow_a[i_O2] = -O2_mass_flow - O2_mem_crossover_reac;
  GDL_an.mass_flow_b = {0, 0, 0};
  GDL_ca.mass_flow_b = {0, 0, 0};
  CL_an.mass_flow_b = -Cell_Area*N_FC_stack/(t_GDL_an*R*T_an)*D_an.*(P_CL_an_ave .- P_GDL_an_ave).*Anode_Medium.MMX;
  CL_ca.mass_flow_b = -Cell_Area*N_FC_stack/(t_GDL_ca*R*T_ca)*D_ca.*(P_CL_ca_ave .- P_GDL_ca_ave).*Cathode_Medium.MMX;
  
  mass_balance[i_H2] = port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_H2]) + port_b_H2.m_flow*actualStream(port_b_H2.Xi_outflow[i_H2]) - H2_mass_flow - H2_mem_crossover;
  
  mass_balance[i_H2O] = port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_H2O]) + port_b_H2.m_flow*actualStream(port_b_H2.Xi_outflow[i_H2O]) + port_a_Air.m_flow*inStream(port_a_Air.Xi_outflow[i_H2O]) + port_b_Air.m_flow*actualStream(port_b_Air.Xi_outflow[i_H2O])+ H2_mass_flow +  O2_mass_flow + H2_mem_crossover + O2_mem_crossover_reac;
  
  mass_balance[i_N2] = port_a_H2.m_flow*inStream(port_a_H2.Xi_outflow[i_N2]) + port_b_H2.m_flow*actualStream(port_b_H2.Xi_outflow[i_N2]) +  port_a_Air.m_flow*inStream(port_a_Air.Xi_outflow[i_N2]) + port_b_Air.m_flow*actualStream(port_b_Air.Xi_outflow[i_N2]);
  
  mass_balance[4] = port_a_Air.m_flow*inStream(port_a_Air.Xi_outflow[i_O2]) + port_b_Air.m_flow*actualStream(port_b_Air.Xi_outflow[i_O2]) - O2_mass_flow -O2_mem_crossover_reac;
  Power_stack = V_cell*N_FC_stack*pin_n.i;
// water transfer
  Nd = 2.5*lambda_m/22;
  EOD_mass_flow = -Nd*Anode_Medium.MMX[i_H2O]*j_cell*Cell_Area/F*N_FC_stack;
  Dnmw = 4.17*10^(-8)*lambda_m*(161*exp(-lambda_m) + 1)*exp(-2346/T_ave);
// Diffusivity of water thought the memebrane is simplified and valid only for lambda > 3 to avoid convergence issues; this is acceptable since normal operation assumes higher membrane humidity.
  diff_mass_flow = Anode_Medium.MMX[i_H2O]*rhoMem/EW*Dnmw*(lambda_ca - lambda_an)/t_mem*Cell_Area*N_FC_stack;
  water_mass_flow_billan_theo = EOD_mass_flow + diff_mass_flow;
  water_mass_flow_billan = min(max(-inStream(GDL_an.port_1.Xi_outflow[i_H2O])*GDL_an.port_1.m_flow*2, water_mass_flow_billan_theo), inStream(GDL_ca.port_1.Xi_outflow[i_H2O])*GDL_ca.port_1.m_flow*2);
// hydrogen and nitrogen crossover trought the membrane
// H2 ans N2 diffusion trought membrane
  mem_N2_diff_0 = 1.163*10^(-12)*exp(3.772*(phi_Anode + phi_Cathode)/2);
  mem_N2_diff_Ea = 1000*11.39*exp(0.5511*(phi_Anode + phi_Cathode)/2);
  mem_N2_diff = mem_N2_diff_0*exp(-mem_N2_diff_Ea/(R*T_ave));
  mem_H2_diff_0 = 9.27*10^(-11)*exp(2.14*(phi_Anode + phi_Cathode)/2)*100/100000;
  mem_H2_diff_Ea = 1000*6.87*exp(0.467*(phi_Anode + phi_Cathode)/2);
  mem_H2_diff = mem_H2_diff_0*exp(-mem_H2_diff_Ea/(R*T_ave));
// diffusion between gaz channel and CL [3]
  D_an = GDL_porosity^1.5*{1.005*10^(-4), 1.005*10^(-4), 1.005*10^(-4)}*(T_an/332.15)^1.5*(101325/CL_an.medium.p);
  D_ca = GDL_porosity^1.5*{2.652*10^(-5), 2.982*10^(-5), 2.982*10^(-5)}*(T_ca/332.15)^1.5*(101325/CL_ca.medium.p);
//*** DEFINE EQUATIONS ***//
// activities
  P_GDL_an_in = Anode_Medium.partialPressure(actualStream(GDL_an.port_1.Xi_outflow), port_a_H2.p);
  P_GDL_an_out = Anode_Medium.partialPressure(actualStream(GDL_an.port_2.Xi_outflow), port_b_H2.p);
  P_GDL_an_ave = (P_GDL_an_in + P_GDL_an_out)/2;
  P_GDL_ca_in = Cathode_Medium.partialPressure(actualStream(GDL_ca.port_1.Xi_outflow), port_a_Air.p);
  P_GDL_ca_out = Cathode_Medium.partialPressure(actualStream(GDL_ca.port_2.Xi_outflow), port_b_Air.p);
  P_GDL_ca_ave = (P_GDL_ca_in + P_GDL_ca_out)/2;
  P_CL_an_ave = Anode_Medium.partialPressure(CL_an.medium.Xi, CL_an.port_2.p);
  P_CL_ca_ave = Cathode_Medium.partialPressure(CL_ca.medium.Xi, CL_ca.port_2.p);
  P_CL_an_in = P_GDL_an_in + (P_CL_an_ave - P_GDL_an_ave);
  P_CL_an_out = P_GDL_an_out + (P_CL_an_ave - P_GDL_an_ave);
  P_CL_ca_in = P_GDL_ca_in + (P_CL_ca_ave - P_GDL_ca_ave);
  P_CL_ca_out = P_GDL_ca_out + (P_CL_ca_ave - P_GDL_ca_ave);
//average activities of H2 and H2O based on logarithmic average of O2 partial pressure
  p_O2_log = (P_GDL_ca_in[i_O2] - P_GDL_ca_out[i_O2])/(log(P_GDL_ca_in[i_O2]) - log(P_GDL_ca_out[i_O2])) + (P_CL_ca_ave[i_O2] - P_GDL_ca_ave[i_O2]);
//(P_CL_ca_in[i_O2] - P_CL_ca_out[i_O2])/(log(P_CL_ca_in[i_O2]) - log(P_CL_ca_out[i_O2]));
// ELECTROCHEMICAL EQUATIONS //
// Calculate the stack voltage
// activation over voltage
  A_eff_pt_an = m_pt_an*ECSA/t_CL_an;
  A_eff_pt_ca = m_pt_ca*ECSA/t_CL_ca;
  theta_T_an = exp(-1400*(1/T_an - 1/353.15));
  theta_T_ca = exp(-7900*(1/T_ca - 1/353.15));
  theta_Pt_im_cov = 1 - exp(-30.6*(r_pt_c/(r_im_c*(1 - r_pt_c)))^2.6);
  theta_Pt_o_cov = 1/(1 + exp(22.4*(0.818 - (V_nernst - V_act_ca - R*T_ca/(2*F)*log((P_CL_an_ave[i_H2]/p_0))))));
  H_H2 = 2583.7875*exp(170/T_an);
  H_O2 = 101325/(4.408 - 0.09712*lambda_m);
  j_cell/t_CL_an = i_0_ref_an*A_eff_pt_an*theta_T_an*(P_CL_an_ave[i_H2]/(H_H2*C_H2_ref))^0.5*(exp((2*F*alpha_an*V_act_an)/(R*T_an)) - exp(-(2*F*alpha_an*V_act_an)/(R*T_an)));
  (j_cell + j_H2_crossover)/t_CL_ca = p_O2_log/H_O2/(C_O2_ref/(i_0_ref_ca*A_eff_pt_ca*theta_T_ca*theta_Pt_im_cov*(1 - theta_Pt_o_cov)*(exp((4*F*alpha_ca*V_act_ca)/(R*T_ca)) - exp(-(4*F*alpha_ca*V_act_ca + w*theta_Pt_o_cov)/(R*T_ca)))));
// stack voltage
  V_rev = 1.229 - 0.85*10^(-3)*(T_ave - 298.15);
  V_nernst = V_rev + R*T_ave/(2*F)*log((P_CL_an_ave[i_H2]/p_0*(p_O2_log/p_0)^0.5)/(P_CL_ca_ave[i_H2O]/p_0));
  V_act = V_act_an + V_act_ca;
  Conductivity = (0.5139*lambda_m - 0.326)*exp(1268*(1/303 - 1/T_ave));
  V_ohm = t_mem/Conductivity*j_cell;
  V_cell = V_nernst - V_act - V_ohm;
  potentialSource.v = N_FC_stack*V_cell;
// THERMAL EQUATIONS //
  P_th = (1.481 - V_cell)*currentSensor.i*N_FC_stack;
  P_in = inStream(port_a_H2.h_outflow)*port_a_H2.m_flow + inStream(port_a_Air.h_outflow)*port_a_Air.m_flow + inStream(port_a_Coolant.h_outflow)*port_a_Coolant.m_flow;
  P_out = actualStream(port_b_H2.h_outflow)*port_b_H2.m_flow + actualStream(port_b_Air.h_outflow)*port_b_Air.m_flow + actualStream(port_b_Coolant.h_outflow)*port_b_Coolant.m_flow;
  CL_an.Q_flow = 0;
  CL_ca.Q_flow = 0;
  GDL_an.Q_flow = 0;
  GDL_ca.Q_flow = 0;
// Assign the thermal power value to the heat flow component
  prescribedHeatFlow_an.Q_flow = fill(V_act_an*j_cell*Cell_Area*N_FC_stack/N_el, N_el);
  prescribedHeatFlow_ca.Q_flow = fill(((141.80*10^6*Anode_Medium.MMX[1]/(2*F) - (V_cell + V_act_an))*j_cell*Cell_Area*N_FC_stack + 119.96*10^6*H2_mem_crossover)/N_el, N_el);
// Efficiencies
  eta_FC_LHV = if currentSensor.i > 0 then potentialSource.v*currentSensor.i/(119.96*10^6*(H2_mass_flow + H2_mem_crossover)) else 0;
//*** DEFINE CONNECTIONS ***//
  connect(pipeCoolant.port_b, port_b_Coolant) annotation(
    Line(points = {{10, -126}, {140, -126}}, color = {255, 0, 0}, thickness = 1));
  connect(pipeCoolant.port_a, port_a_Coolant) annotation(
    Line(points = {{-10, -126}, {-140, -126}}, color = {0, 0, 255}, thickness = 1));
  connect(pin_n, potentialSource.n) annotation(
    Line(points = {{-60, 150}, {-60, 150}, {-60, 130}, {-60, 130}}, color = {0, 0, 255}));
  connect(pin_n, ground.p) annotation(
    Line(points = {{-60, 150}, {-100, 150}, {-100, 138}, {-100, 138}}, color = {0, 0, 255}));
  connect(channelCathode.port_b, port_b_Air) annotation(
    Line(points = {{130, -84}, {130, -102}, {150, -102}}, color = {0, 127, 255}));
  connect(port_a_Air, inlet_air_temperature.port) annotation(
    Line(points = {{150, 80}, {106, 80}, {106, 94}}));
  connect(outlet_air_temperature.port, port_b_Air) annotation(
    Line(points = {{205, -101}, {152, -101}, {152, -102}, {150, -102}}, color = {0, 127, 255}));
  connect(channelAnode.port_b, port_b_H2) annotation(
    Line(points = {{-120, -84}, {-120, -100}, {-150, -100}}, color = {0, 127, 255}));
  connect(port_a_H2, inlet_hydrogen_temperature.port) annotation(
    Line(points = {{-150, 80}, {-120, 80}, {-120, 88}}));
  connect(potentialSource.p, currentSensor.p) annotation(
    Line(points = {{-60, 110}, {-60, 100}, {-10, 100}}, color = {0, 0, 255}));
  connect(water_prod.y, water_prod_and_transfer.u1) annotation(
    Line(points = {{2, 39}, {2, 19}, {24, 19}}, color = {0, 0, 127}));
  connect(currentSensor.n, pin_p) annotation(
    Line(points = {{10, 100}, {60, 100}, {60, 150}}, color = {0, 0, 255}));
  connect(CL_an.mass_flow_a[1], water_prod.u2) annotation(
    Line(points = {{-55, 14}, {-49, 14}, {-49, 78}, {-4, 78}, {-4, 62}}, color = {0, 0, 127}));
  connect(CL_ca.mass_flow_a[1], water_prod.u1) annotation(
    Line(points = {{75, 14}, {47, 14}, {47, 78}, {8, 78}, {8, 62}}, color = {0, 0, 127}));
  connect(CL_an.mass_flow_a[2], water_prod_and_transfer.u2) annotation(
    Line(points = {{-55, 14}, {-38, 14}, {-38, 8}, {24, 8}}, color = {0, 0, 127}));
  connect(water_prod_and_transfer.y, CL_ca.mass_flow_a[2]) annotation(
    Line(points = {{47, 14}, {75, 14}}, color = {0, 0, 127}));
  connect(CL_an.mass_flow_a[3], nitogen_cross_over.u) annotation(
    Line(points = {{-55, 14}, {-52, 14}, {-52, -7}, {-9, -7}}, color = {0, 0, 127}));
  connect(nitogen_cross_over.y, CL_ca.mass_flow_a[3]) annotation(
    Line(points = {{14, -7}, {50, -7}, {50, 14}, {75, 14}}, color = {0, 0, 127}));
  connect(port_a_H2, GDL_an.port_1) annotation(
    Line(points = {{-150, 80}, {-120, 80}, {-120, 23}}));
  connect(GDL_an.port_2, channelAnode.port_a) annotation(
    Line(points = {{-120, 3}, {-120, -64}}, color = {0, 127, 255}));
  connect(port_a_Air, GDL_ca.port_1) annotation(
    Line(points = {{150, 80}, {130, 80}, {130, 24}}));
  connect(GDL_ca.port_2, channelCathode.port_a) annotation(
    Line(points = {{130, 4}, {130, -64}}, color = {0, 127, 255}));
  connect(GDL_an.mass_flow_a, diff_trought_GDL_an.u) annotation(
    Line(points = {{-116, 14}, {-96, 14}}, color = {0, 0, 127}, thickness = 0.5));
  connect(diff_trought_GDL_an.y, CL_an.mass_flow_b) annotation(
    Line(points = {{-82, 14}, {-64, 14}}, color = {0, 0, 127}, thickness = 0.5));
  connect(CL_ca.mass_flow_b, diff_trought_GDL_ca.y) annotation(
    Line(points = {{84, 14}, {98, 14}}, color = {0, 0, 127}, thickness = 0.5));
  connect(diff_trought_GDL_ca.u, GDL_ca.mass_flow_a) annotation(
    Line(points = {{112, 14}, {126, 14}}, color = {0, 0, 127}, thickness = 0.5));
  connect(heatCapacitor.port, thermalConductor.port_b) annotation(
    Line(points = {{0, -70}, {0, -90}}, color = {191, 0, 0}, thickness = 0.5));
  connect(thermalConductor.port_a, pipeCoolant.heatPorts) annotation(
    Line(points = {{0, -110}, {0, -122}}, color = {191, 0, 0}, thickness = 0.5));
  connect(thermalConductor_an.port_b, heatCapacitor.port) annotation(
    Line(points = {{-40, -74}, {0, -74}, {0, -70}}, color = {191, 0, 0}, thickness = 0.5));
  connect(thermalConductor_ca.port_a, heatCapacitor.port) annotation(
    Line(points = {{40, -74}, {0, -74}, {0, -70}}, color = {191, 0, 0}, thickness = 0.5));
  connect(thermalConductor_ca.port_b, channelCathode.heatPorts) annotation(
    Line(points = {{60, -74}, {126, -74}}, color = {191, 0, 0}, thickness = 0.5));
  connect(thermalConductor_an[end:-1:1].port_a, channelAnode.heatPorts) annotation(
    Line(points = {{-60, -74}, {-116, -74}}, color = {191, 0, 0}, thickness = 0.5));
  connect(channelAnode.heatPorts, thermalConductor_mem[end:-1:1].port_a) annotation(
    Line(points = {{-116, -74}, {-102, -74}, {-102, -36}, {-10, -36}}, color = {127, 0, 0}, thickness = 0.5));
  connect(thermalConductor_mem.port_b, channelCathode.heatPorts) annotation(
    Line(points = {{10, -36}, {100, -36}, {100, -74}, {126, -74}}, color = {191, 0, 0}, thickness = 0.5));
  connect(channelAnode.heatPorts, prescribedHeatFlow_an.port) annotation(
    Line(points = {{-116, -74}, {-84, -74}, {-84, -86}}, color = {127, 0, 0}, thickness = 0.5));
  connect(prescribedHeatFlow_ca.port, channelCathode.heatPorts) annotation(
    Line(points = {{92, -92}, {92, -74}, {126, -74}}, color = {191, 0, 0}, thickness = 0.5));
  connect(heatPort, heatCapacitor.port) annotation(
    Line(points = {{48, -50}, {22, -50}, {22, -74}, {0, -74}, {0, -70}}, color = {191, 0, 0}, thickness = 0.5));
  connect(outlet_hydrogen_temperature.port, port_b_H2) annotation(
    Line(points = {{-194, -98}, {-150, -98}, {-150, -100}}, color = {0, 127, 255}));
  connect(outlet_coolant_temperature.port, port_b_Coolant) annotation(
    Line(points = {{180, -144}, {140, -144}, {140, -126}}, color = {0, 127, 255}));
  connect(inlet_coolant_temperature.port, port_a_Coolant) annotation(
    Line(points = {{-170, -150}, {-140, -150}, {-140, -126}}, color = {0, 127, 255}));
  annotation(
    Diagram(coordinateSystem(extent = {{-150, -150}, {150, 150}}, initialScale = 0.1)),
    Icon(coordinateSystem(extent = {{-150, -150}, {150, 150}}, initialScale = 0.1), graphics = {Line(origin = {20.1754, 1.92106}, points = {{0, 78}, {0, -80}, {0, -82}}), Rectangle(origin = {80, 0}, fillColor = {0, 178, 227}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-20, 100}, {20, -100}}), Line(origin = {40.1315, 2}, points = {{0, 78}, {0, -80}, {0, -82}}), Line(origin = {0.219199, 1.92106}, points = {{0, 78}, {0, -80}, {0, -82}}), Line(origin = {-40.0001, 1.61404}, points = {{0, 78}, {0, -80}, {0, -82}}), Rectangle(origin = {-80, 0}, fillColor = {170, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-20, 100}, {20, -100}}), Text(origin = {10, -54}, textColor = {255, 0, 0}, extent = {{-11, 6}, {11, -6}}, textString = "K"), Line(origin = {-20.0439, -0.307018}, points = {{0, 80}, {0, -80}, {0, -80}}), Rectangle(origin = {35, 54}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Vertical, extent = {{-95, 26}, {25, -134}}), Text(origin = {-80, 6}, extent = {{-26, 24}, {26, -24}}, textString = "A"), Text(origin = {80, 6}, extent = {{-26, 24}, {26, -24}}, textString = "C"), Text(origin = {82, 110}, extent = {{-22, 10}, {22, -10}}, textString = "OER"), Text(origin = {75, -113}, extent = {{-21, 11}, {21, -11}}, textString = "phi"), Text(origin = {-86, -94}, extent = {{-22, 10}, {22, -10}}, textString = "OER"), Text(origin = {-85, -113}, extent = {{-21, 11}, {21, -11}}, textString = "phi")}),
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
    experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 1),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts -d=aliasConflicts",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end FuelCellStack;
