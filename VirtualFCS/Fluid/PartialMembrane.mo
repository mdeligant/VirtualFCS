within VirtualFCS.Fluid;

partial model PartialMembrane
  import Modelica.Fluid.Types;
  import Modelica.Fluid.Types.PortFlowDirection;
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component" annotation(
    choicesAllMatching = true);
  Modelica.Fluid.Interfaces.FluidPort_a port_1(redeclare package Medium = Medium, m_flow(min = if (portFlowDirection_1 == PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf, max = if (portFlowDirection_1 == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) annotation(
    Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
  Modelica.Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = Medium, m_flow(min = if (portFlowDirection_2 == PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf, max = if (portFlowDirection_2 == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) annotation(
    Placement(transformation(extent = {{90, -10}, {110, 10}})));
protected
  parameter PortFlowDirection portFlowDirection_1 = PortFlowDirection.Bidirectional "Flow direction for port_1" annotation(
    Dialog(tab = "Advanced"));
  parameter PortFlowDirection portFlowDirection_2 = PortFlowDirection.Bidirectional "Flow direction for port_2" annotation(
    Dialog(tab = "Advanced"));
  public  annotation(
    Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics = {Text(textColor = {0, 0, 255}, extent = {{-150, -89}, {150, -129}}, textString = "%name"), Rectangle(origin = {0, -23}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -23}, {100, 23}}), Rectangle(origin = {0, 23}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Vertical, extent = {{-100, -23}, {100, 23}})}),
    Diagram(graphics));
end PartialMembrane;
