within VirtualFCS;

package Media
  package MoistAirThreeComponents
    "Moist air with detailed concentration of Nitrogen and Oxygen"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "FlueGasSixComponents", data = {Modelica.Media.IdealGases.Common.SingleGasesData.N2, Modelica.Media.IdealGases.Common.SingleGasesData.O2, Modelica.Media.IdealGases.Common.SingleGasesData.H2O}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.N2, Modelica.Media.IdealGases.Common.FluidData.O2, Modelica.Media.IdealGases.Common.FluidData.H2O}, substanceNames = {"Nitrogen", "Oxygen", "Water" }, reference_X = {0.78, 0.21, 0.01});
  
  
  
  
  
  function relativeHumidity
    "Return relative humidity as a function of the thermodynamic state record"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "Thermodynamic state";
    output Real phi "Relative humidity";
  algorithm
    phi := state.X[3]*state.p / Modelica.Media.Water.StandardWater.saturationPressure(state.T);
  
  
  end relativeHumidity;
  
  
  
  
  annotation(
      Documentation(info = "<html>
  
  </html>"));
end MoistAirThreeComponents;
  extends Modelica.Icons.Package;

end Media;