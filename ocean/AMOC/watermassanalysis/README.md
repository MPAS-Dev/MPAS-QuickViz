# Water Mass Analysis

Development space for code and visualization related to water mass analysis in MPAS-Ocean.

The primary method used in this analysis is the Walin air-sea flux framework ([Walin 1982, Tellus](https://doi.org/10.3402/tellusa.v34i2.10801); [Speer and Tziperman 1992](https://doi.org/10.1175/1520-0485(1992)022<0093:ROWMFI>2.0.CO;2)). Briefly, the mode water transformation $F$ and formation $M$ are defined in terms of the net air/sea density flux $\Phi$ area-integrated over and isopycnal outcrop.

$$\underbrace{F(\rho) = \frac{\partial}{\partial\rho}\iint_{\rho_0}^{\rho}\Phi_{\rho}ds}_{\text{Transformation}} \hspace{1cm} \underbrace{M(\rho) = - \frac{\partial F}{\partial\rho}}_{\text{Formation}}$$