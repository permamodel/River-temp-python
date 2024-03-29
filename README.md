# Rivertemp

River water temperature model 

The river water temperature model is designed to be applied in Arctic rivers. Heat energy transfers considered include surface net solar radiation, net longwave radiation, latent heat due to evaporation and condensation, convective heat and the riverbed heat flux. 

You can run the model with the River_temp_portal.

Please note that the default riverbed deposit is sand and gravel, users can change the sedimentary material with alternative soil/substrate physical and thermal parameters. 

Inputs------------------------------------------------------------

His: incident solar (shortwave) radiation (w m-2)   --- Mfile 

Hdl: Downward thermal (Longwave) radiation (w m-2)  --- Mfile 

Ta: Air temperature (°C)                            --- Mfile 

Uw:  Wind speed (m/s)                               --- Mfile 

rh: Relative humidity (0-1)                         --- Mfile 

P: Surface pressure (pa)                            --- Mfile 

h:  River stage (m)                                 --- Sfile 

rivericeoff: First day of river induation (DOY)

rivericeon: Last day of river induation (DOY)

Output-----------------------------------------------------------

Tw: Water temperature (°C)


The underlying theory and assumptions of this model are published in JGR-ES: 
Changing Arctic River Dynamics Cause Localized Permafrost Thaw, Lei Zheng  Irina Overeem  Kang Wang  Gary D. Clow, 2019. Journal of Geophysical Research - Earth Surface. https://doi.org/10.1029/2019JF005060
