### Setting Parameters in the EFT Structure

There are two methods to set parameters in the EFT structure:

1. **Case-by-Case Definition**

   This method is suitable for specific models and functions in the structure. Parameters have physical significance within their corresponding models.

   Exact examples will also be provided in the examples folder, and we list their names here.

   * Flag for different $w_{DE}$ parametrizations (LCDM, wCDM, CPL, JBL, Taylor or TurningPoint) in PureEFTmode: **EFTwDE**.
        - **EFTwDE**&emsp;  = 0&emsp; ->&emsp; $w_{DE} = -1$
        - **EFTwDE**&emsp;  = 1&emsp; ->&emsp; $w_{DE} = w_0$
        - **EFTwDE**&emsp;  = 2&emsp; ->&emsp; $w_{DE} = w_0 + w_a(1-a)$
        - **EFTwDE**&emsp;  = 3&emsp; ->&emsp; $w_{DE} = w_0 + w_a (1-a) a^{(n-1)}$
        - **EFTwDE**&emsp;  = 5&emsp; ->&emsp; $w_{DE} = w_0 + w_a a + \frac12 w_2 a^2 + \frac16 w_3 a^3$
        - **EFTwDE**&emsp;  = 4&emsp; ->&emsp; $w_{DE} = w_0 + w_a (a_t-a)^2$
    
   The parameters above can be fixed with the flags:
     - $w_0$->**EFTw0** ,  $w_a$->**EFTwa** , $n$->**EFTwn**, $a_t$->**EFTwat** , $w_2$->**EFTw2** , $w_3$->**EFTw3**

   * Parameter names of **AltParEFTmodel**:
   
   **AltParEFTmodel** = 1
   - Flag for different $w_{DE}$ parametrizations (LCDM, wCDM, CPL, JBL, Taylor or TurningPoint): **RPHwDE**. The parameter names are the same with **EFTwDE**, except changeing "EFT" instead of "RPH".
   - Flag for parametrizations of RPH functions: **RPHHmassPmode, RPHkineticitymodel, RPHbraidingmodel, RPHtensormodel**.
   
   Parametrizations name rule see the text below. Their function names are **RPHalphaM, RPHmassP, RPHkineticity, RPHbraiding, RPHtenso**

   and Latex format are $\alpha^{\rm M},\tilde{M},\alpha^{\rm K},{\alpha^{\rm B}},{\alpha^{\rm T}}$, respectively.

   **AltParEFTmodel** = 2
   
   Flags: **OLLambdamodel, OLOmegamodel, OLGamma1model, OLGamma2model, OLGamma3model**

   Names: **OLLambda, OLOmega, OLGamma1, OLGamma2, OLGamma3**
   
   Latex format: $\tilde{\lambda}, \Omega, \gamma_1, \gamma_2, \gamma_3$

   **AltParEFTmodel** = 3
   
   Flags: **OLLambdamodel, OLmassPmodel, OLkineticitymodel, OLbraidingmodel, OLtensormodel**
   
   Names: **OLLambda, OLmass, OLkineticity, OLbraiding, OLtensor**
   
   Latex format: $\tilde{\lambda}, \tilde{M}, \alpha^{\rm K}, \alpha^{\rm B}, \alpha^{\rm T}$   

   **AltParEFTmodel** = 4
   
   - Flag for different $w_{DE}$ parametrizations: **EFTwDE**.





   

3. **Default Definition**

   If there is no specific definition for the model, a default name will be used. This name is constructed as "Function name" + index, where "Function name" is defined in the model flag, and the index ranges from 0 to \( N-1 \), where \( N \) is the number of parameters. The LaTeX format will be: \[\text{"FunctionName"}_{index}\]

### Parametrizations for EFT functions

In the flowchart, you will notice several functions labeled with "Parametrizations." We have implemented a number of parametrization methods within EFTCAMB. The flag numbers and corresponding parameter names are listed here. Still, if there is no specific definition for a parameter name, the default definition will be used.

**We use $\Omega$ in pureEFT as an example**

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
    
