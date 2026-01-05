within VirtualFCS.Fluid;

model Membrane "Splitting/joining component with static balances for a dynamic control volume"
  extends Fluid.PartialMembrane;
  extends Modelica.Fluid.Interfaces.PartialLumpedVolume(final fluidVolume = V);
  //replaceable package Medium=Modelica.Media.Interfaces.PartialMedium   ;
  parameter Modelica.Units.SI.Volume V "Mixing volume inside junction";
  Modelica.Units.SI.HeatFlowRate Q_flow;
  Modelica.Units.SI.MolarMass M_in;
  Modelica.Units.SI.MolarMass M_out;
  Modelica.Blocks.Interfaces.RealInput mass_flow_a[Medium.nXi](each unit = "kg/s") annotation(
    Placement(transformation(origin = {-50, 114}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {-1, 47}, extent = {{-13, -13}, {13, 13}}, rotation = -90)));   
  Modelica.Blocks.Interfaces.RealInput mass_flow_b[Medium.nXi](each unit = "kg/s") annotation(
    Placement(transformation(origin = {-48, -103}, extent = {{-20, -20}, {20, 20}}, rotation = 90), iconTransformation(origin = {1, -46}, extent = {{-13, -13}, {13, 13}}, rotation = 90)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation(
    Placement(transformation(origin = {28, -98}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {48, 50}, extent = {{-10, -10}, {10, 10}})));
equation
// Only one connection allowed to a port to avoid unwanted ideal mixing
  assert(cardinality(port_1) <= 1, "
port_1 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
  ");
  assert(cardinality(port_2) <= 1, "
port_2 of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
  ");
// Boundary conditions
  port_1.h_outflow = medium.h;
  port_2.h_outflow = medium.h;
  port_1.Xi_outflow = medium.Xi;
  port_2.Xi_outflow = medium.Xi;
  port_1.C_outflow = C;
  port_2.C_outflow = C;
// Mass balances
  mb_flow = port_1.m_flow + port_2.m_flow + sum(mass_flow_a + mass_flow_b) "Mass balance";
  mbXi_flow = port_1.m_flow*actualStream(port_1.Xi_outflow) + port_2.m_flow*actualStream(port_2.Xi_outflow) +  mass_flow_a + mass_flow_b "Component mass balances";
  mbC_flow = port_1.m_flow*actualStream(port_1.C_outflow) + port_2.m_flow*actualStream(port_2.C_outflow) "Trace substance mass balances";
// Momentum balance (suitable for compressible media)
  port_1.p = medium.p;
  port_2.p = medium.p;
// Energy balance
  Hb_flow = port_1.m_flow*actualStream(port_1.h_outflow) + port_2.m_flow*actualStream(port_2.h_outflow) +      mass_flow_a[1]*Medium.specificEnthalpy_pTX(p = medium.p, T = medium.T, X = {1,0,0}) + mass_flow_b[1]*Medium.specificEnthalpy_pTX(p = medium.p, T = medium.T, X = {1,0,0}) +
mass_flow_a[2]*Medium.specificEnthalpy_pTX(p = medium.p, T = medium.T, X = {0,1,0}) + mass_flow_b[2]*Medium.specificEnthalpy_pTX(p = medium.p, T = medium.T, X = {0,1,0}) +
mass_flow_a[3]*Medium.specificEnthalpy_pTX(p = medium.p, T = medium.T, X = {0,0,1}) + mass_flow_b[3]*Medium.specificEnthalpy_pTX(p = medium.p, T = medium.T, X = {0,0,1})  ; 
  port_a.T = medium.T;
  Qb_flow = port_a.Q_flow + Q_flow;
  Wb_flow = 0;
  M_in = 1/(sum(actualStream(port_1.Xi_outflow)./Medium.MMX));
  M_out = 1/(sum(actualStream(port_2.Xi_outflow)./Medium.MMX));
  annotation(
    Placement(transformation(origin = {-18, 102}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-37, 101}, extent = {{-21, -21}, {21, 21}})),
    Icon);
end Membrane;
