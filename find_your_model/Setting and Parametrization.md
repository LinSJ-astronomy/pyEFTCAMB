### Setting Parameters in the EFT Structure

There are two methods to set parameters in the EFT structure:

1. **Case-by-Case Definition**

   This method is suitable for specific models and functions in the structure. Parameters have physical significance within their corresponding models.

   For example, consider the designer \( f(R) \) and the background dark energy equation of state.

   Exact examples will be provided in the examples folder, and we will list them here in the future.

2. **Default Definition**

   If there is no specific definition for the model, a default name will be used. This name is constructed as "Function name" + index, where "Function name" is defined in the model flag, and the index ranges from 0 to \( N-1 \), where \( N \) is the number of parameters. The LaTeX format will be: 

   \[
   \text{"FunctionName"}_{index}
   \]

### Parametraizations for EFT functions

    In the flowchart, you will notice several functions labeled with "Parametrizations." We have implemented a number of parametrization methods within EFTCAMB. The flag numbers and corresponding parameter names are listed here. Still, if there is no specific definition for a parameter name, the default definition will be used.

   **We use $\Omega$ in pureEFT as an example
      - Flags for $\Omega$ parametrizations: **PureEFTmodelOmega** 
      - Zero: **PureEFTmodelOmega** = 0&emsp; ->&emsp; $\Omega(a) = 0$
      - Constant: **PureEFTmodelOmega** = 1&emsp; ->&emsp; $\Omega(a) = \Omega_0$
      - Linear: **PureEFTmodelOmega** = 2&emsp; ->&emsp; $\Omega(a) = \Omega_0 a$
      - Power Law: **PureEFTmodelOmega** = 3&emsp; ->&emsp; $\Omega(a) = \Omega_0 a^s$. Parameter Names: $\Omega_0$-> **EFTOmega0** , $s$ -> **EFTOmegaExp**. Latex format label: ($\Omega_0$, $\Omega_n$).
      - Exponential: **PureEFTmodelOmega** = 4&emsp; ->&emsp; $\Omega(a) = \exp(\Omega_0 a^s) -1$. Parameter Names: $\Omega_0$-> **EFTOmega0** , $s$ -> **EFTOmegaExp**. Latex format label: ($\Omega_0$, $\Omega_n$).
      - Taylor series: **PureEFTmodelOmega** = 5&emsp; ->&emsp; 
      - Pade series: **PureEFTmodelOmega** = 6&emsp; ->&emsp; 
      - Fourier: **PureEFTmodelOmega** = 7&emsp; ->&emsp; 
      - Steplog: **PureEFTmodelOmega** = 8&emsp; ->&emsp; 
      - Spline: **PureEFTmodelOmega** = 9&emsp; ->&emsp;
      - Spline5: **PureEFTmodelOmega** = 10&emsp; ->&emsp;
      - Exponential_Parametrization_2_1D: **PureEFTmodelOmega** = 11&emsp; ->&emsp;
    
