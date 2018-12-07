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

:ref:`ColRadPy <colradpy>` requires multiple inputs from the user to run.
The basic quasistatic run will be discussed first. This should be sufficient for most users.
Once the code runs data is store in a dictionary accessed through the '.data' method.
The entries of this dictionary and descriptions are documented :ref:`here.<colradpy_dict>`

* Atomic data input file - currently this is limitted to an ADF04 file but there is nothing special about it.

* array of metastable levels - This is an array of levels that includes the ground and any levels that could be considered metastable.

* temperature grid - This is an array of electron temperatures for the calculation in eV.

* Density grid     - This is an array of electron densities for the calculation in cm-3.




In this example Be II (Be+) is used because it is a simple system that has a parent ion
(the next charge state) that has a metastable.
This allows all of the different functionality to be shown and tested.


#. The module  must be executed
#. Input file, temperature, density and metastable inputs defined.
#. Ionization, recombination options chosen default values are true


.. code-block:: python
   :linenos:
      
   %run colradpy_class.py
   fil = 'cpb03_ls#be0.dat' #adf04 file
   temperature_arr = np.array([10,50,100]) #eV
   metastable_levels = np.array([0])   #metastable level, just ground chosen here
   density_arr =     np.array([1.e13,4.e14]) # cm-3
   be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
                 use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True)


'be' is now a colradpy class, there are many different calls that could be made from the class :ref:`documented here <colradpy>`.
The general flow of the code goes as


#. create rates
#. populate matrix
#. solve matrix and create post processed data.

   
.. code-block:: python
   :linenos:

   if(be.data['user']['use_ionization']):
       be.make_ioniz_from_reduced_ionizrates()
   if(be.data['user']['suppliment_with_ecip']):
       be.make_ecip()
       be.suppliment_with_ecip()
   if(be.data['user']['use_recombination']):
       be.make_recombination_rates_from_file()
   if(be.data['user']['use_recombination_three_body']):
       be.make_three_body_recombination()
   be.make_electron_excitation_rates()
       
   
.. code-block:: python
   :linenos:

   be.populate_cr_matrix()


.. code-block:: python
   :linenos:

   be.solve_quasi_static()
   
   
This can be done from one call to colradpy. Which does the procedure above.

.. code-block:: python
   :linenos:

   be.solve_cr_qausistatic()


Data from the calculation is now avaible in the '.data' dictionary.
Various postpocessing can be done to now analysis the calcuation.



Post processing analysis
------------------------



Photon emissivity coefficients (PECs)
--------------------------------------

A theortical spectrum can be made from the PEC coefficients.
PEC coefficient are stored in array that has shape (#pecs,metastable,temperature,density).
The code below produces a PEC spectrum for on temperature and density.
The wavelength and pec arrays share the same length.

.. code-block:: python
   :linenos:
      
   import matplotlib.pyplot as plt
   plt.ion()
   met = 0 #metastable 0, this corresponds to the ground state
   te = 0 #first temperature in the grid
   ne = 0 #frist density in the grid
   plt.plot(be.data['processed']['wave_air'],be.data['processed']['PECS'][:,met,te,ne])
   
   
Often the index of a specific pec is wanted to find its temperature or density dependence.
This can be accomplished in two basic ways.

#. Upper and lower levels of the transitions are known
#. The wavelength of the transition is known

There is a map from transition numbers to pec index levels. .data['processed']['pec_levels'] has
the same order as .data['processed']['wave_air'] and .data['processed']['pecs'].


.. code-block:: python
   :linenos:

   print(np.shape(be.data['processed']['wave_air']),
         np.shape(be.data['processed']['pecs']),
	 np.shape(be.data['processed']['pec_levels']))
   #(320,) (320, 3, 1, 1) (320, 2)

   upper_ind = 10 #ninth excited state
   lower_ind = 0  #ground state

   pec_ind = np.where( (be.data['processed']['pec_levels'][:,0] == upper_ind) &\
                       (be.data['processed']['pec_levels'][:,1] == lower_ind))[0]

   plt.figure()
   #plot the temeprature dependence of the chosen pec at first density in the grid
   plt.title('Temperature dependence of line ' +\
              str(be.data['processed']['wave_air'][pec_ind]) +' nm')
   plt.plot(be.data['user']['temp_grid'],be.data['processed']['pecs'][pec_ind,met,:,ne])
   plt.xlabel('Temperature (eV)')
   plt.ylabel('PEC (ph cm-3 s-1)')

   plt.figure()
   #plot the density dependence of the chosen pec at first density in the grid
   plt.title('Density dependence of line ' +\
              str(be.data['processed']['wave_air'][pec_ind]) +' nm')
   plt.plot(be.data['user']['dens_grid'],be.data['processed']['pecs'][pec_ind,met,te,:])
   plt.xlabel('Density (cm-3)')
   plt.ylabel('PEC (ph cm-3 s-1)')   


If the wavelength of a line of interest is known, the index can be found by looking at the
wavelength array.
The indices of all pecs that fall within the upper and lower bound of the 'where' statement are
returned. PECs can generally be distinguished by the actual value, large lines that are of interest
have much large PEC values, this can allow 

      
.. code-block:: python
   :linenos:

   #want to find the index of Be I line at 351.55
   pec_ind = np.where( (be.data['processed']['wave_air'] <352) &\
                       (be.data['processed']['wave_air'] >351))
   print('Wavelength from file ' + str(be.data['processed']['wave_air'][pec_ind[0]]))
   print('PEC upper and lower levels '+ str(be.data['processed']['pec_levels'][pec_ind[0]]))
   

Generalized radiative coefficients (GCRs)
-----------------------------------------

The generalized collsional radiative coefficients are calculated by ColRadPy as well.
A description of these can be found in (Summers 2006), (Johnson 2019).
GCR coefficients are often use as inputs to plasma transport codes.
GCR coefficients are also use as inputs to ionization balance calculations which will be discussed
later. This allows for different ionization stages to be linked.




.. code-block:: python
   :linenos:
   #plotting the QCD
   plt.figure()
   plt.plot(be.data['user']['temp_grid'],
            be.data['processed']['qcd'][0,0,:,0],
	    label = 'metastable cross coupling coefficient 1->2')
	    
   plt.plot(be.data['user']['temp_grid'],
            be.data['processed']['qcd'][1,0,:,0],
	    label = 'metastable cross coupling coefficient 2->1')
   plt.legend()
   plt.xlabel('Temperature (eV)')
   plt.ylabel('QCD (cm-3 s-1)')
   

.. code-block:: python
   :linenos:
   #plotting the SCD
   plt.figure()
   plt.plot(be.data['user']['temp_grid'],
            be.data['processed']['scd'][0,0,:,0],
	    label = 'metastable cross coupling coefficient 1->1')
	    
   plt.plot(be.data['user']['temp_grid'],
            be.data['processed']['scd'][1,0,:,0],
	    label = 'metastable cross coupling coefficient 2->1')
   plt.legend()
   plt.xlabel('Temperature (eV)')
   plt.ylabel('SCD (ion cm-3 s-1)')
