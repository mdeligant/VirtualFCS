within VirtualFCS.Control;

model PID "PID controller with anti-windup"
  import SI = Modelica.Units.SI;

  parameter Real Kp "Proportional gain (normalised units)";
  parameter SI.Time Ti "Integral time";
  parameter Boolean integralAction = true "Use integral action";
  parameter Real PVmin "Minimum value of process variable for scaling";
  parameter Real PVmax "Maximum value of process variable for scaling";
  parameter Real CSmin "Minimum value of control signal for scaling";
  parameter Real CSmax "Maximum value of control signal for scaling";
  parameter SI.Time Td = 0 "Derivative time";
  parameter Real Nd = 1 "Derivative action up to Nd / Td rad/s";
  parameter Real Ni = 1 "Ni*Ti is the time constant of anti-windup compensation";
  parameter Real b = 1 "Setpoint weight on proportional action";
  parameter Real c = 0 "Setpoint weight on derivative action";
  parameter Real PVstart = 0.5 "Start value of PV (scaled)";
  parameter Real CSstart = 0.5 "Start value of CS (scaled)";
  parameter Boolean holdWhenSimplified = false "Hold CSs at start value when homotopy=simplified";
  parameter Boolean steadyStateInit = false "Initialize in steady state";
  Real CSs_hom "Control signal scaled in per units, used when homotopy=simplified";
  Real P "Proportional action / Kp";
  Real I(start = CSstart / Kp) "Integral action / Kp";
  Real D "Derivative action / Kp";
  Real Dx(start = c * PVstart - PVstart) "State of approximated derivator";
  Real PVs "Process variable scaled in per unit";
  Real SPs "Setpoint variable scaled in per unit";
  Real CSs(start = CSstart) "Control signal scaled in per unit";
  Real CSbs(start = CSstart) "Control signal scaled in per unit before saturation";
  Real track "Tracking signal for anti-windup integral action";
  Modelica.Blocks.Interfaces.RealInput PV "Process variable signal" annotation(
    Placement(transformation(extent = {{-112, -52}, {-88, -28}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput CS "Control signal" annotation(
    Placement(transformation(extent = {{88, -12}, {112, 12}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput SP "Set point signal" annotation(
    Placement(transformation(extent = {{-112, 28}, {-88, 52}}, rotation = 0)));
equation
// Scaling
  SPs = (SP - PVmin) / (PVmax - PVmin);
  PVs = (PV - PVmin) / (PVmax - PVmin);
  CS = CSmin + CSs * (CSmax - CSmin);
// Controller actions
  P = b * SPs - PVs;
  if integralAction then
    assert(Ti > 0, "Integral time must be positive");
    Ti * der(I) = SPs - PVs + track;
  else
    I = 0;
  end if;
  if Td > 0 then
    Td / Nd * der(Dx) + Dx = c * SPs - PVs "State equation of approximated derivator";
    D = Nd * (c * SPs - PVs - Dx) "Output equation of approximated derivator";
  else
    Dx = 0;
    D = 0;
  end if;
  if holdWhenSimplified then
    CSs_hom = CSstart;
  else
    CSs_hom = CSbs;
  end if;
  CSbs = Kp * (P + I + D) "Control signal before saturation";
  CSs = homotopy(smooth(0, if CSbs > 1 then 1 else if CSbs < 0 then 0 else CSbs), CSs_hom) "Saturated control signal";
  track = (CSs - CSbs) / (Kp * Ni);
initial equation
  if steadyStateInit then
    if Ti > 0 then
      der(I) = 0;
    end if;
    if Td > 0 then
      D = 0;
    end if;
  end if;
  annotation(
    Documentation(info = "<html>
<p><i>Cet élément est fonctionnel.</i></p>
<p>Ce contrôleur PID est issu de la librairie <a href=https://build.openmodelica.org/Documentation/ThermoPower.html>Thermopower</a>. Par défaut, c'est un PI.</p>
<p> Les principaux paramètres à ajuster sont (voir l'exemple dans le modèle CAES_Eleves.Tests.test_CombustionChamber_with_PID) :
<ul>
<li> Le gain Kp (normalisé);</li>
<li> Le temps d'intégration Ti;</li>
<li> Le booléen integralAction pour activer l'intégrateur (actif par défaut);</li>
<li> Les valeurs minimales et maximales attendues de la grandeur physique à réguler (PVmin et PVmax), pour la normalisation;</li>
<li> Les valeurs minimales et maximales du signal physique à fournir en sortie (CSmin et CSmax), pour la normalisation;</li>
<li> Le temps du dérivateur Td: par défaut Td=0, ce qui désactive le dérivateur.</li>
</ul>
</p>
<p> Ce contrôleur est lui-même inspiré par l'élément <a href=https://build.openmodelica.org/Documentation/Modelica.Blocks.Continuous.LimPID.html> Modelica.Blocks.Continuous.LimPID</a> de la librairie standard.</p>
<p> Voir par exemple <a href=https://fr.wikipedia.org/wiki/R%C3%A9gulateur_PID>PID</a> pour des rappels.</p>
</html>"));
  annotation(
    Diagram(graphics),
    Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-54, 40}, {52, -34}}, lineColor = {0, 0, 255}, textString = "PID"), Text(extent = {{-110, -108}, {110, -142}}, lineColor = {0, 0, 255}, lineThickness = 0.5, textString = "%name")}));
end PID;