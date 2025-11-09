within VirtualFCS;

package Media
  package MoistAirThreeComponents
    "Moist air with detailed concentration of Nitrogen and Oxygen"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "MoistAirThreeComponents", data = {Modelica.Media.IdealGases.Common.SingleGasesData.N2, Modelica.Media.IdealGases.Common.SingleGasesData.O2, Modelica.Media.IdealGases.Common.SingleGasesData.H2O}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.N2, Modelica.Media.IdealGases.Common.FluidData.O2, Modelica.Media.IdealGases.Common.FluidData.H2O}, substanceNames = {"Nitrogen", "Oxygen", "Water" }, reference_X = {0.78, 0.21, 0.01});
  
  
  
  
  
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
  
  package MoistHydrogenThreeComponents
    "Moist hydrogen with detailed concentration of Nitrogen and steam"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "MoistHydrogenThreeComponents", data = {Modelica.Media.IdealGases.Common.SingleGasesData.H2, Modelica.Media.IdealGases.Common.SingleGasesData.N2, Modelica.Media.IdealGases.Common.SingleGasesData.H2O}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.H2, Modelica.Media.IdealGases.Common.FluidData.N2, Modelica.Media.IdealGases.Common.FluidData.H2O}, substanceNames = {"Hydrogen", "Nitrogen", "Water"}, reference_X = {0.7, 0.3, 0.0});
  
    function relativeHumidity "Return relative humidity as a function of the thermodynamic state record"
      extends Modelica.Icons.Function;
      input ThermodynamicState state "Thermodynamic state";
      output Real phi "Relative humidity";
    algorithm
      phi := state.X[3]*state.p/Modelica.Media.Water.StandardWater.saturationPressure(state.T);
    end relativeHumidity;
    annotation(
      Documentation(info = "<html>
    
    </html>"));
  end MoistHydrogenThreeComponents;
  extends Modelica.Icons.Package;

end Media;