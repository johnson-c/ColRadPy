=============
Tutorial
=============
This tutorial will take the user through the basics of atomic data needed for CR modeling.
Give an overview of how to run ColRadPy.

Neither the local thermodynamic equilibrium (LTE) or conornal approximations are valid for fusion plasmas.
These plasma must be modeled with a collisional radiative (CR) model which includes effects of electron density.
Collisional Radiative theory has been updated to included general collisional radiative (GCR) theory, this is
outlined in (Summers 2006).

Atomic Data and the ADF04 File
==============================
Atomic data is generally made by calculating cross-sections (probablity) of transitions.
From these cross-sections rates and then made by convolving with a maxwellian temperature distriubtion of electron
energies with the cross section. The rates are what is used in the CR model.

The ADF04 file format from ADAS that is used to hold atomic data for input into the CR model. The ADF04 file is supported as
an input to colradpy.

Electron impact de/excitation
-----------------------------
Electron impact excitation and deexcitation is what causes transitions be the ground and excited states and within excited
states. There rates are stored as "effective collisional strengths" in the adf04 file. These can then be converted to the
excitiation and deexcitation rates that are needed for CR modeling.

Electron ionization
-------------------
Electrons can also ionize atoms to the next charge state. Ionization can be stored in the adf04 file as reduced ionization
rates. ColRadPy is also capable of making ECIP ionization rates. ECIP ionization is a very approximate way of making ionization
and should be avoid if possible. Ground and metastable states generally have calculated ionization rates while excited states
generally have rates that are supplimented with ECIP.

Radative Recombination
----------------------


Dielectronic Recombination
----------------------


Three Body Recombination
------------------------
Three body recombination is calculated from using detailed balance from the ionization


Running ColRadPy
================
ColRadPy requires multiple inputs from the user to run.

Atomic data input file - currently this is limitted to an ADF04 file but there is nothing special about it.

array of metastable levels - This is an array of levels that includes the ground and any levels that could be considered metastable.

temperature grid - This is an array of electron temperatures for the calculation in eV.

Density grid     - This is an array of electron densities for the calculation in cm-3.



In this example Be II (Be+) is used because it is a simple system that has a parent ion (the next charge state) that has a metastable.
This allows all of the different functionality to be shown and tested.


First the mode must be executed and temperature, density and metastable inputs defined. Ionization, recombination

.. code-block:: python
   :linenos:
      
   %run colradpy_class.py
   fil = 'cpb03_ls#be0.dat' #adf04 file
   temperature_arr = np.array([10,50,100]) #eV
   metastable_levels = np.array([0])   #metastable level, just ground chosen here
   density_arr =     np.array([1.e13,4.e14]) # cm-3
   be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
                 use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True)


'be' is now a colradpy class, there are many different calls that could be made from the class :ref:`documented here <colradpy_class>`

   


