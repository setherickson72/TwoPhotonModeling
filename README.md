# TwoPhotonModeling

A group of codes for modeling two-photon transitions in a variety of atoms. No guarantee of accuracy, but hopefully this is helpful to someone

* RepRateCalculator: For calculating the apparent seperation in frequency of all two-photon transitions expected in your frequency window. You must manually enter the frequency of the transitions within the "fs" array, and the corresponding linestrength in the "as".

* CombExcitationCalc_ManyDifferentLengths: For calculating the effect of different bandwidths for direct comb excitation. Hidden behind the simulation is the assumption that bandwidth is significantly smaller than 1 THz (intermediate state detuning)

* ResidualDopplerBroadening: For a given spectral bandwidth FWHM, calculates the effective lineshape accounting for residual doppler broadening. Run this simulation several times to find what spectral bandwidths result in minimal residual doppler broadening.

* TransitTime: A very simple code for calculating the transit time broadening for a given temperature and beam size. Technically works for any atom, but the variables are preset for Rb87

* CombExcitationCalc_GDD: Similar to CombExcitationCalc_ManyDifferentLengths, but allows for simulation over a range of chirps.
