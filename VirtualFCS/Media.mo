within VirtualFCS;

package Media
  package MoistAirThreeComponents
    "Moist air with detailed concentration of Nitrogen and Oxygen"
   extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "MoistAirThreeComponents", data = {Modelica.Media.IdealGases.Common.SingleGasesData.O2, Modelica.Media.IdealGases.Common.SingleGasesData.H2O, Modelica.Media.IdealGases.Common.SingleGasesData.N2}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.O2, Modelica.Media.IdealGases.Common.FluidData.H2O, Modelica.Media.IdealGases.Common.FluidData.N2}, substanceNames = {"Oxygen", "Water", "Nitrogen"}, reference_X = {0.15, 0.07, 0.78});
  
  
  
  
  
  function relativeHumidity
    "Return relative humidity as a function of the thermodynamic state record"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "Thermodynamic state";
    output Real phi "Relative humidity";
    protected Real phi_v[MoistAirThreeComponents.nX];
    algorithm
      phi_v := MoistAirThreeComponents.massToMoleFractions(state.X, MoistAirThreeComponents.MMX)*state.p/ Modelica.Media.Water.StandardWater.saturationPressure(state.T);
      phi := phi_v[2];
  end relativeHumidity;
  
  function partialPressure
    extends Modelica.Icons.Function;
    input MassFraction X[MoistAirThreeComponents.nX];
    input Modelica.Units.SI.Pressure P;
    output Modelica.Units.SI.Pressure Px[MoistAirThreeComponents.nX];
    algorithm
      Px := P*X./MoistAirThreeComponents.MMX/(sum(X./MoistAirThreeComponents.MMX));
   end partialPressure;
  
  
  
  annotation(
      Documentation(info = "<html>
  
  </html>"));
end MoistAirThreeComponents;
  
  package MoistHydrogenThreeComponents
    "Moist hydrogen with detailed concentration of Nitrogen and steam"
    extends Modelica.Media.IdealGases.Common.MixtureGasNasa(mediumName = "MoistHydrogenThreeComponents", data = {Modelica.Media.IdealGases.Common.SingleGasesData.H2, Modelica.Media.IdealGases.Common.SingleGasesData.H2O, Modelica.Media.IdealGases.Common.SingleGasesData.N2}, fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.H2, Modelica.Media.IdealGases.Common.FluidData.H2O, Modelica.Media.IdealGases.Common.FluidData.N2}, substanceNames = {"Hydrogen", "Water", "Nitrogen"}, reference_X = {0.7, 0.2, 0.1});
    
  function relativeHumidity "Return relative humidity as a function of the thermodynamic state record"
    extends Modelica.Icons.Function;
    input ThermodynamicState state "Thermodynamic state";
    output Real phi "Relative humidity";
    protected Real phi_v[MoistAirThreeComponents.nX];
    algorithm
      phi_v := MoistHydrogenThreeComponents.massToMoleFractions(state.X, MoistHydrogenThreeComponents.MMX)*state.p/Modelica.Media.Water.StandardWater.saturationPressure(state.T);
      phi := phi_v[2];
  end relativeHumidity;
    
  function partialPressure
    extends Modelica.Icons.Function;
    input MassFraction X[MoistHydrogenThreeComponents.nX];
    input Modelica.Units.SI.Pressure P;
    output Modelica.Units.SI.Pressure Px[MoistHydrogenThreeComponents.nX];
    algorithm
      Px := P*X./MoistHydrogenThreeComponents.MMX/(sum(X./MoistHydrogenThreeComponents.MMX));
  end partialPressure; 
    
    annotation(
      Documentation(info = "<html>
    
    </html>"));
  end MoistHydrogenThreeComponents;
  extends Modelica.Icons.Package;

end Media;
