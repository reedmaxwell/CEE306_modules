# %%
import math
# %% [markdown]
#Reference ET calculation for FAO-56
#===================================
# This notebook follows the [FAO-56] (http://www.fao.org/3/X0490E/X0490E00.htm) methodology which implements the Penman
# approach to ET, to get a reference ET0 (often called potential ET) value
# it then uses crop coefficients to scale ET0 to some actual ET
#
# Overall Game Plan
#------------------
# This document outlines the FAO-56 approach which is useful when there are limited observations.
# We will calculate a theoretical solar radation budget, we will calculate a vapor pressure based 
# only on temperature (min), then we will calculate the slope of the vapor pressure curve ($\delta$)
# and psychrometric constant ($\gamma$) which are foundational components to the Penman equation.
# Each of the supporting calculations are arranged in sections.  If you have observations for (e.g.)
# long and shortwave net radiation, you would generally use those in place of the theoretical 
# calculations I show here.
#
#Calculate dew-point temperature, $T_{dew}$
#-------------------------------------
# This derives $T_{dew}$ from the relative humidity using FAO-56 Eq(10). This is done
# when there are no observations of relative humidity and is an estimate based on the 
# minimum daily temperature, $T_min$ (this one reason is why many products, e.g. DayMET, report 
# min/max temps).  $T_{dew} = T_{min} + K_0$ 
# where the empirical dewpoint offset, $K_0 = 2$ \[deg C].
# %%
Tmin = 15  #deg C
K0=2 #dewpoint offset, empirical
# %%
Tdew=Tmin-K0
print('Tdew=',Tdew,' deg C')
# %% [markdown]
#Calculate actual vapor pressure, $e_a$
#-----------------------------------
# The actual vapor pressure is the pressure when air is cooled to be saturated.  This is the 
# saturation vapor pressure at the dewpoint temp, $T_{dew}$.  We use FAO-56 Eq(14) to calculate
# this, it is $e_a=e_0(T_{dew})$ \[kPa] or $e_a=0.6108*e^{(17.27*T_{dew})/(T_{dew}+237.3)}$
# %%
ea=0.6108*math.exp((17.27*Tdew)/(Tdew+237.3))  #[kPa]
print('ea=',ea)
# %% [markdown]
#Calculate the net solar radiation
#---------------------------------
# A majority of the calculations in FAO-56 and in codes like PRMS are to estimate the net solar
# radiation.  This uses the equations we developed earlier in the class for theoretical net solar short 
# and long wave radiation (e.g. 2.6.1 in Brutsaert).ß    
# %%
