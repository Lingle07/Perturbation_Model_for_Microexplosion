# Perturbation Model for Microexplosion
Hello everyone. This is a code developed by Xinyi for a UGVR project on microexplosion. The basic idea is explained in the poster file.
## Basic wish
Zengbing Yang [2] has detailedly described in his PhD thesis about the perturbation model for predicting the onset of microexplosion.
However, though the model is exact and detailed, the great simulation cost makes it hard to give a fast prediction for the onset of microexplosion. The best way is to simplify it through assumption and approximation. This makes the prediction limited and we need to be careful.
## Models used
### Big idea - Dynamic model
Since the whold model is based on the work of Chia-fon Lee's group, we would like to make it clear how it work. Using the equations of momentum, temperature and diffusion, assuming the droplet and all its properties and process are all homogeneous, the equations can be simplified into 3, with boundary conditions. Still, many properties should be decided by the temperature, which make all equations strongly coupled.
We wish to reduce the calculation cost, so it's best to approximate a bit. Since the calculation is based on the dynamic process of gas bubble, we try to use only the dynamic discription of the bubble, and temperature profile will be approximated to obtain properties.
### Model of bubble growth - Mikic
The model of bubble growth has been thoroughly exammined before by a lot of research groups. Mathematicians [6] have studied some oscillating analytical solutions of the bubble, which is based on that there is no diffusion and that gas inside the bubble will not change its mass. At the same time, using Clausius-Clapeyron relation, which assumes that the gas is saturated, Mikic [5] developed a model for approximating its growth, which is especially helpful when the diffusion reaches any limitation.
For droplet growth, there should be more possible models which carefully incorporate the diffusion properties and the mixing of components. 
This model has omitted the boundary condition of the outer layer, which is easy to be added into the model since we can always get Rayleigh-Plesset equation with slight changes due to the boundary conditions. Then we will have more specific equations similar to Mikic's solution for this case.
Another important part which might need consideration is the evaporation of the inner layer, which is really important for temperature profile, but also might have some influence on pressure.
Also, there is still no good enough way for how properties of mixing fuels behave. That should be easier to solve, since that is chemistry things.
### Model for temperature profile
There is no proper idea of how to develop a temperature profile for all conditions. Basically, the progress is 1-dimensional. Using some approximation, there might be a proper evaluation for how gas bubble change.

# Acknowledgement
Great thanks to Pavan Govindaraju and Prof. Matthias Ihme in Stanford University for their instructions and help.

If you have any questions, please contact Xinyi by <hxy176@126.com>.

1. Zeng, Y., & Chia-fon, F. L. (2007). Modeling droplet breakup processes under micro-explosion conditions. *Proceedings of the Combustion Institute*, 31(2), 2185-2193.
2. Zeng, Y. (2000). Modeling of multicomponent fuel vaporization in internal combustion engines *(Doctoral dissertation, University of Illinois at Urbana-Champaign)*.
3. Law, C. K. (1977). A model for the combustion of oil/water emulsion droplets. *Combustion Science and Technology*, 17(1-2), 29-38.
4. Kadota, T., Tanaka, H., Segawa, D., Nakaya, S., & Yamasaki, H. (2007). Microexplosion of an emulsion droplet during Leidenfrost burning. *Proceedings of the Combustion institute*, 31(2), 2125-2131.
5. Mikic, B. B., Rohsenow, W. M., & Griffith, P. (1970). On bubble growth rates. *International Journal of Heat and Mass Transfer*, 13(4), 657-666.
6. Mehrem, A. M. Estudio teo rico de la dina mica de microburbujas bajo la accio n de un campo ultraso nico.
7. Lien, Y. C. (1969). Bubble growth rates at reduced pressure *(Doctoral dissertation, Massachusetts Institute of Technology)*.
