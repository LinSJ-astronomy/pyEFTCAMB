===================
EFTCAMB MODELS
===================

This file lists all the models and model selection flags that can be used with EFTCAMB:

- EFTflag = 0: standard CAMB with GR.

- EFTflag = 1: EFTCAMB in Pure EFT model. 
  
  - PureEFTmodel = 1: standard formulation with Omega and Gammas

- EFTflag = 2: EFTCAMB with alternative base parametrizations.

  - AltParEFTmodel = 1: alternative parametrization with Omega and alphas
  - AltParEFTmodel = 2: alternative parametrization with Omega, Lambda and Gammas
  - AltParEFTmodel = 3: alternative parametrization with Omega, Lambda and alphas
  - AltParEFTmodel = 4: shift symmetric alpha B

- EFTflag = 3: EFTCAMB designer models.

  - DesignerEFTmodel = 1: designer f(R)
  - DesignerEFTmodel = 2: designer minimally coupled quintessence
  - DesignerEFTmodel = 3: designer DGP

- EFTflag = 4: EFTCAMB full mapping models.

  - FullMappingEFTmodel = 1: Horava
  - FullMappingEFTmodel = 2: Acoustic Dark Energy
  - FullMappingEFTmodel = 3: K-mouflage
  - FullMappingEFTmodel = 4: Quintessence
  - FullMappingEFTmodel = 5: Beyond Honrdeski
  - FullMappingEFTmodel = 6: Scaling Cubic Galileon
