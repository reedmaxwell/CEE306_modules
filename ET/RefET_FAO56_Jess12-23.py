#!/usr/bin/env python
# coding: utf-8
#Draft for CEE 306 module
#Reference ET calculation for FAO-56
#Editing and expanding on ref_et_calc.py from 12/22/2021

#Still need to incorporate crop functions. Ended at ET0 today (12/23)
# In[1]:


import math

# This notebook follows the [FAO-56] (http://www.fao.org/3/X0490E/X0490E00.htm) methodology which implements the Penman
# approach to ET, to get a reference ET0 (often called potential ET) value
# it then uses crop coefficients to scale ET0 to some actual ET
# Overall Game Plan
# This document outlines the FAO-56 approach which is useful when there are limited observations.
# We will calculate a theoretical solar radation budget, we will calculate a vapor pressure based 
# only on temperature (minimum and maximum), then we will calculate the slope of the vapor pressure curve ($\delta$)
# and psychrometric constant ($\gamma$) which are foundational components to the Penman equation.
# Each of the supporting calculations are arranged in sections.  If you have observations for (e.g.)
# long and shortwave net radiation, you would generally use those in place of the theoretical 
# calculations I show here.
# Calculate dew-point temperature, $T_{dew}$
# Dew-point temperature is the temperature at which a parcel of air is fully saturated 
# with water so that droplets (dew) start to form.
# $T_{dew}$ can be derived from the relative humidity using FAO-56 Eq(10). 
# When there are no observations of relative humidity, we can use an estimate based on the 
# minimum daily temperature, $T_min$ (this is one reason is why many products, e.g. DayMET, report 
# min/max temps).  $T_{dew} = T_{min} + K_0$ 
# where the empirical dewpoint offset, $K_0 = 2$ \[deg C].
# In[2]:


Tmin = 15  #deg C
K0=2 #dewpoint offset, empirical


# In[3]:


Tdew=Tmin+K0
print('Tdew=',Tdew,' deg C')


# Calculate actual vapor pressure, $e_a$
# Vapor pressure is the partial pressure of water vapor within air. Saturation vapor pressure
# is the vapor pressure of an air parcel at some temperature when it has plenty of water. As seen
# in FAO-56 Eq(10), actual vapor pressure is the saturation vapor pressure times relative humidity.
# Since dew-point temperature is the temp at which air is saturated, then the actual vapor pressure
# at the dew-point temperature is equal to the saturation vapor pressure at that temperature.
# We use FAO-56 Eq(14) to calculate this as
# $e_a=e_0(T_{dew})$ \[kPa] or $e_a=0.6108*e^{(17.27*T_{dew})/(T_{dew}+237.3)}$
# In[4]:


ea=0.6108*math.exp((17.27*Tdew)/(Tdew+237.3))  #[kPa]
print('ea=',ea,' kPa')


# Calculate the saturation vapor pressure, $e_{s}$, and then the vapor pressure deficit, $e_{s}-e_{a}$
# Because the saturation vapor pressure equation is non-linear, the average should first calculate
# the saturation vapor pressure at the minimum and maximum temperature, and then average the results.
# $e_s=(e_0(T_{max})+e_0(T_{min}))/2$ \[kPa]
# Vapor pressure deficit is then the difference between the saturated and actual vapor pressures.
# In[5]:


Tmax = 25 #[C] maximum temperature
e0_Tmax=0.6108*math.exp((17.27*Tmax)/(Tmax+237.3))  #[kPa]
e0_Tmin=0.6108*math.exp((17.27*Tmin)/(Tmin+237.3))  #[kPa]
es=(e0_Tmax+e0_Tmin)/2 #[kPa]
print('es=',es,' kPa')
print('VPD=',es-ea,' kPa')


# Calculate the slope of the vapor pressure curve, $\Delta$
# The slope of the vapor pressure curve, $\Delta$, or D, is quite simply the derivative
# or rise/run of the FAO-56 Eq(11) for saturation vapor pressure as a function of temperature.
# D in units of kPa per degree C is listed in FAO-56 Eq(13) as
# $D=4096*(0.6108*math.exp(17.27*T/(T+237.3)))/(T+237.3)**2$
# We can take the slope at a mean temperature.
# In[6]:


Tmean = (Tmax+Tmin)/2 #[C]
D=4096*(0.6108*math.exp(17.27*Tmean/(Tmean+237.3)))/(Tmean+237.3)**2
print('$\Delta$=',D,' kPa/C')


# Calculate the net solar radiation
# A majority of the calculations in FAO-56 and in codes like PRMS are to estimate the net solar
# radiation.  This uses the equations we developed earlier in the class for theoretical net solar short 
# and long wave radiation (e.g. 2.6.1 in Brutsaert).ß   
# For daily periods, extraterrestrial radiation can be calculated using the solar constant, the 
# relative Earth-Sun distance, the sunset hour angle, the latitude, and the solar declination angle
# FAO-56 Eq(21)
# In[7]:


#Extraterrestrial radiation = Ra
Gsc = 0.0820 #[MJ/(m^2 min)], Solar Constant
J = 200 #day of year after Jan 1 (less than 366)
dr = 1+0.033*math.cos(2*math.pi*J/365) #relative Earth-Sun distance
sd = 0.409*math.sin(2*math.pi*J/365 - 1.39) #[rad], solar declination angle
j = math.pi/4 #[rad], latitude
ws = math.acos(-math.tan(j)*math.tan(sd)) #[rad], sunset hour angle
Ra = (24*60/math.pi)*Gsc*dr*(ws*math.sin(j)*math.sin(sd)+math.cos(j)*math.cos(sd)*math.sin(ws))
print('Ra=',Ra,' MJ/(m^2*day)')

# So that is how much solar radiation is reaching the top of atmosphere, but what is actually
# touching the ground? We need to estimate how much is lost in the atmosphere.
#FAO-56 Eq(35)
# In[8]:


N = 24*ws/math.pi #maximum daylight hours given sunshine hour angle
print('N = ',N,' hours') #check N before deciding n


# In[9]:


n = 12 #actual sunshine hours (must be less than N)
a_s = 0.25 #recommended Angstrom constant
b_s = 0.55 #recommended Angstrom constant
Rs = (a_s+b_s*(n/N))*Ra #Solar Radiaton at ground
print('Rs=',Rs,' [MJ/(m^2*day)]')

# if n=N (and it is clear sky all day) then
# Rso = (a_s + b_s)*Ra
# In[10]:


#Clear sky solar radiation
Rso = (a_s+b_s)*Ra
print('Rso=',Rso,' [MJ/(m^2*day)]')

# Finally, net solar radiation is what comes to the ground minus what gets reflected.
# Albedo tells us what fraction gets reflected, so net short waave radiation is 
# $R_{ns} = Rs*(1-albedo)$
# In[11]:


albedo = 0.3
Rns = Rs*(1-albedo)
print('Rns=',Rns,' [MJ/(m^2*day)]')


# Calculate the net longwave radiation
# FAO-56 Eq(39) details an approach to estimating net long-wave radiation out. Long-wave radiation is
# energy emitted at Earth temperatures (much cooler than the Sun). Net long-wave is 
# radiation emitted up by the surface minus radiation emitted down by the atmosphere. The atmospheric
# emission depends strongly on water content, so the equation uses $e_a$.
# In[12]:


s = 4.903E-9 #MJ/(K^4*m^2) Stefan Boltzmann Constant
Rnl = s*((Tmean+273.16)**4)*(0.34-0.14*math.sqrt(ea))*(1.35*Rs/Rso-0.35)
print('Rnl=',Rnl,' [MJ/(m^2*day)]')


# Calculate net radiation
# Net radiation is simply the balance of short and long-wave radiation. Since our short-wave was
# positive into the ground and our long-wave was positive out of the ground, 
# FAO-56 Eq(40) writes Rn = Rns - Rnl
# In[13]:


Rn = Rns - Rnl
print('Rn=',Rn,'[MJ/(m^2*day)]')

# G is soil heat flux. For this exercise, we will take it to be zero. 
# Calculate the wind speed 2 meters above the ground, $u_2$
# Wind tends to have a logarithmic profile with zero wind exactly at the ground, and faster
# winds higher up. Our equation requires us to have the wind 2 meters above the ground.
#FAO-56 Eq(47) gives the wind at 2 meters given u_z and z
# In[14]:


uz = 4 #[m/s], wind at not 2 m
z = 10 #[m], height at which uz was measured
u2 = uz*4.87/math.log(67.8*z-5.42)
print('u2=',u2,' m/s')


# Calculate reference evapotranspiration, $ET_0$
# FAO-56 Eq(6) uses the same units for the above terms as given here to get a reference ET in mm/day.
# The last piece needed is the psychrometric constant (which varies a little with temperature).
# Dingman 3rd Edition Eq(6.21) describes this, and lists a common value as 0.066 kPa/K
# In[15]:


gamma = 0.066 #[kPa/C] = [kPa/K]
ET0 = (0.408*D*Rn+gamma*(900/(Tmean+273))*u2*(es-ea))/(D+gamma*(1+0.34*u2))
print('ET_0=',ET0,' mm/day')

# By this final potential ET equation, we see that increased radiation, wind speed, and 
# vapor pressure deficit all increase the potential for the atmosphere to evaporate moisture from
# the ground surface or from plants.