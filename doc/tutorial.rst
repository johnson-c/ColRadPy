.. _tutorial:
=============
Tutorial
=============
This tutorial will take the user through the basics of atomic data needed for Collisional Radiative (CR) modeling.
The python scripts for these examples can be found in the 'examples' folder.

For a detailed introduction to CR theory see the ColRadPy paper included in the 'user_doc' folder.
Summers 2006 also provides a more detailed overview of CR theory.

Neither the local thermodynamic equilibrium (LTE) or conornal approximations are valid
for fusion plasmas.
These plasmas must be modeled with a collisional radiative (CR) model which includes effects
of electron density.
Collisional Radiative theory has been updated to included general collisional radiative (GCR) theory,
this is outlined in (Summers 2006).



The tutrial will go over the following major topics

Installing ColRadPy `Installing ColRadPy`_


The atomic data input for colradpy is overviewed in `Atomic Data and the ADF04 File`_

Running ColRadPy will be outlined in `Running ColRadPy`_

Built in post processes in ColRadPy is described in `Post processing analysis built in functions`_

Posted processing for advanced users will be outlined in `Post processing analysis advanced`_

The generalized collisional radiative coefficients will be discussed in `Generalized radiative coefficients (GCRs)`_

Functionality for advanced users is outlined in `Advanced functionality`_

Ionization balance is outlined in `Ionization Balance`_



Installing ColRadPy
===================
ColRadPy should work with any python 3.6 distribution.
However, ColRadPy is only tested with the anaconda distribution of python 3.6.
ColRadPy requires numpy, matplotlib, re, scipy, sys, collections packages.



Atomic Data and the ADF04 File
==============================
Atomic data is generally made by calculating cross-sections (probablity) of transitions.
From these cross-sections rates are then made by convolving with a maxwellian temperature distriubtion (or some of distribution) of electron
energies with the cross section.
The rates produced by this method are what is used in the CR model.

The ADF04 file format from `ADAS <http://www.adas.ac.uk/>`_ that is used to hold atomic rate data for input into the CR model.




Electron Impact De/Excitation
-----------------------------
Electron impact excitation and deexcitation is what causes transitions be the ground and excited states and within excited
states. There rates are stored as "effective collisional strengths" in the adf04 file.
These can then be converted to the excitiation and deexcitation rates that are needed for CR modeling.


Electron Impact Ionization
----------------------------
Electrons can also ionize atoms to the next charge state. Ionization can be stored in the adf04 file as reduced ionization
rates. ColRadPy is also capable of making ECIP ionization rates. ECIP ionization is a very approximate way of making ionization
and should be avoid if possible. Ground and metastable states generally have calculated ionization rates while excited states
generally have rates that are supplimented with ECIP.

Radative Recombination
-------------------------


Dielectronic Recombination
----------------------------

Dielectronic recombination is the dominant recombination process in many plasmas.
Dielectronic recombination is a two step process.
First, there is a resonance capture of a free electron.
In this processes a bound electron gains some energy and the free electron loses some energy.
Second, the atom radiates a phonton.


Three-Body Recombination
--------------------------
Three-body recombination is calculated from using detailed balance from the ionization.
In this process, free electrons and a ion recombine to produce a new ion and one free electron.


Solutions
-----------
Once rates from the basic atomic data have been input.
ColRadPy can then solve the set of rate equations known as the collisional radiative equations.
There are two ways that ColRadPy can solve the CR set of equations.
The general user will probably be interested in the first way of solving the CR equations, this is same way ADAS solves the equations.
The Quasistatic approximation solves the equations assuming that excited states do not have
time changing populations. It assumes that these populations and in their equilbrium values with the
ground and any metastable levels. This is the way that ADAS solves the CR set of equations.

ColRadPy is also capable of solving the system of equations time dependently using the matrix exponentialation method.
This solution makes no assumptions of equilbrium time scales for excited states.
This is an excact solution.
In equilbrium, the non-quasistatic and quastistic solutions will be the same.





Running ColRadPy
===================

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




There are are highlevel class definitions that a general user can use.
These high level definitions will be shown as well as the more basic definitions that
these high level definitions call.
A basic user needs to only care about the high level definitions

Below is an example of running ColRadPy.
The user chooses an input adf04 file, temperature and density grids as well as the number of metastables.
If you don't know how many metastables exist only choose level '0' (the ground) as being metastable.
Choices were also made to include ionization and recombination.
Note that metastable levels are only important if the plasma is not in equilbrium.
In equilbrium metastable fraction is set and will not vary for a given temperature and density.

The 'use_recombination' flag when true will use any recombination rates that are included in the adf04 file.

The 'use_recombination_three_body' flag when true will have ColRadPy make and include three body recombination rates.

Note inorder to have three body recombination there must be some ionization included in the calculation.

The 'use_ionization' flag when true will use any ionization rates that are included in the adf04 file.

The 'suppliment_with_ecip' when true will have ColRadPy make ECIP ionization rates and include these rates anywhere
that there are no ionization rates included in the adf04 file.


Note that ionization should always be included in the calculation even if it is just ECIP.
Ionization can significantly change both absolute values of PECs as well as relative values of PECs.


.. code-block:: python
   :linenos:

    import sys
    sys.path.append('../') #starting in 'examples' so need to go up one
    from colradpy_class import colradpy
    import numpy as np

    fil = 'cpb03_ls#be0.dat' #adf04 file
    temperature_arr = np.linspace(1,100,100) #eV
    metastable_levels = np.array([0])   #metastable level, just ground chosen here
    density_arr =     np.array([1.e13,4.e14]) # cm-3

    #calling the colradpy class with the various inputs
    be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True, 
		  use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True)

    be.solve_cr() #solve the CR equations with the quasistatic method
    

'be' is now a colradpy class that has been solved. There are various methods for getting the data out.
Data that required ColRadPy to solve the CR set of equations is now stored in the 'processed' sub dictionary.
There are many different calls that could be made from the class :ref:`documented here <colradpy>`.



.. hidden-code-block:: python
    :linenos:
    :label: --- Show/details of solve_cr()---

	    
    #The flow of the code is below
    # 1. create rates from adf04 or internal
    # 2. populate matrix
    # 3. solve matrix and create post processed data.
	    
    #the block of code below is what 'solve_cr()' is doing
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

    # This call takes in all of the rates that were provided by
    # the adf04 file as well as any rates that were bade in
    # ColRadPy and puts them into the matrix.
    
    be.populate_cr_matrix()

    # This call solves the CR matrix and creates the 'processed'
    # sub dictionary with all of the processed data.
    
    be.solve_quasi_static()


Data from the calculation is now avaible in the '.data' dictionary.
Various postpocessing can be done to now analysis the calcuation.



Post processing analysis built in functions
=============================================
There are some built in functions avaible for post processing analysis
these will of use for the general user.
The basic functions  require minimal knowelge of the underlying datastructure.
These basic functions will be overviewed first then a more complex analysis will be presented after.

The theorical spectrum can be plotted describe in `Plotting Theorical Spectrum (PEC sticks)`_

The line ratios versus temperature and density is describe in `Plotting PEC ratios`_



Plotting Theorical Spectrum (PEC sticks)
------------------------------------------
The theorical spectral spectrum from the adf04 file can be plotted with the below command.
The parameters are the lists or arrays of the index of the metastable,
temperature and density grids.
*WARNING* Note that wavelengths will not match NIST wavelengths
unless the adf04 energy levels have been shifted to the NIST values.
This generally hasn't been done in the past so there are many adf04 files that don't
use NIST energy values.

.. code-block:: python
   :linenos:

      
      be.plot_pec_sticks([0],[0],[0])

      

.. hidden-code-block:: python
    :linenos:
    :label: --- Show/details of plot_pec_sticks()---


            
      rc('axes', linewidth=2)
      rc('font', weight='semibold')

      if('processed' not in self.data.keys()):
          self.solve_cr()

      p_t = np.arange(0,len(self.data['user']['temp_grid']))
      p_n = np.arange(0,len(self.data['user']['dens_grid']))
      p_m = np.arange(0,len(self.data['atomic']['metas']))
      if(np.asarray(temp).size>0):
          p_t = p_t[temp]
      if(np.asarray(dens).size>0):
          p_n = p_n[dens]
      if(np.asarray(meta).size>0):
          p_m = p_m[meta]

      for i in p_n:
          for j in p_t:
              for k in p_m:
                  plt.figure()
                  scaling = int(np.floor(np.log2(np.max(self.data['processed']['pecs'][:,k,j,i]))/np.log2(10)))
                  plt.vlines(self.data['processed']['wave_air'],
                         np.zeros_like(self.data['processed']['wave_air']),
                                       self.data['processed']['pecs'][:,k,j,i]*10**np.abs(scaling))

                  plt.xlabel('Wavelength in air (nm)',weight='semibold')
                  plt.ylabel('PEC X 1E' +str(scaling) + ' (ph cm$^{-1}$ s$^{-1}$)',weight='semibold')
                  plt.title('Temperature ' + str(self.data['user']['temp_grid'][j]) + ' eV,  '+\
                            'Density ' + format(self.data['user']['dens_grid'][i],'.2e') + ' cm$^{-3}$, '+\
                            'Metastable ' + str(self.data['atomic']['metas'][k]),weight='semibold')
                  plt.xlim(0,1300)

      

Plotting PEC ratios
---------------------
Spectral line intenties are functions of both electron temperature and density as well as ion density.
Different spectral lines will have different functional forms on temperature and density.
It is therefore possible to find ratios of spectral lines that are depenended on either temperature or density as the ion density cancels out.
It is then possible to dianose electron temperature and density from line ratios where the charge state exists in the plasma.



Temperature ratios can be plotted with the below function.
The temperature ratio of the two pecs will be plotted with new figures made for each density and metastable requested.
The inputs are the indexes of pec1, pec2, array of densities and array of metastables

.. code-block:: python
   :linenos:

      
      be.plot_pec_ratio_temp(0,1,[0],[0])



.. hidden-code-block:: python
    :linenos:
    :label: --- Show/details of plot_pec_ratio_temp()---


        if('processed' not in self.data.keys()):
            self.solve_cr()

        
        dens = np.array(dens)
        p_n = np.arange(0,len(self.data['user']['dens_grid']),dtype=int)
        p_m = np.arange(0,len(self.data['atomic']['metas']),dtype=int)
        
        if(np.asarray(dens).size>0):
            p_n = p_n[dens]
        if(np.asarray(meta).size>0):
            p_m = p_m[meta]

        for k in p_m:
            plt.figure()                        
            for i in p_n:
                
                plt.plot(self.data['user']['temp_grid'],
                     self.data['processed']['pecs'][pec1,k,:,i]/ \
                         self.data['processed']['pecs'][pec2,k,:,i],
                         label='$n_e$ = ' + format(self.data['user']['dens_grid'][i],'.1e') + ' cm$^{-3}$')

                plt.xlabel('Temperature (eV)',weight='semibold')
                plt.ylabel('Ratio (-)',weight='semibold')
                plt.title('Ratio of PEC '+str(pec1)+', ' + format(self.data['processed']['wave_air'][pec1],'.2f') + ' nm'+\
                          ' to PEC '+str(pec2)+', ' + format(self.data['processed']['wave_air'][pec2],'.2f') + ' nm, '+\
                              'Metastable ' + str(self.data['atomic']['metas'][k]),weight='semibold')
                plt.legend(loc='best')


.. code-block:: python
   :linenos:

      
      be.plot_pec_ratio_dens([0],[0],[0])
      


.. hidden-code-block:: python
    :linenos:
    :label: --- Show/details of plot_pec_ratio_dens()---


        temp = np.array(temp)
        p_n = np.arange(0,len(self.data['user']['temp_grid']),dtype=int)
        p_m = np.arange(0,len(self.data['atomic']['metas']),dtype=int)
        
        if(np.asarray(temp).size>0):
            p_n = p_n[temp]
        if(np.asarray(meta).size>0):
            p_m = p_m[meta]

        for k in p_m:
            plt.figure()                        
            for i in p_n:
                
                plt.plot(self.data['user']['dens_grid'],
                     self.data['processed']['pecs'][pec1,k,i,:]/ \
                         self.data['processed']['pecs'][pec2,k,i,:],
                         label='$T_e$ = ' + format(self.data['user']['temp_grid'][i],'.1f') + ' eV')

                plt.xlabel('Density (cm$^{-3}$)',weight='semibold')
                plt.ylabel('Ratio (-)',weight='semibold')
                plt.title('Ratio of PEC '+str(pec1)+', ' + format(self.data['processed']['wave_air'][pec1],'.2f') + ' nm'+\
                          ' to PEC '+str(pec2)+', ' + format(self.data['processed']['wave_air'][pec2],'.2f') + ' nm, '+\
                              'Metastable ' + str(self.data['atomic']['metas'][k]),weight='semibold')
                if(scale=='log'):
                    plt.semilogx()
                plt.legend(loc='best')


            




Post processing analysis advanced
====================================

Photon emissivity coefficients (PECs)
----------------------------------------

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

   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.105,right=0.965)
   ax1.vlines(be.data['processed']['wave_air'],
	      np.zeros_like(be.data['processed']['wave_air']),
	      be.data['processed']['pecs'][:,met,te,ne])
   ax1.set_xlim(0,1000)
   ax1.set_title('PEC spectrum  T$_e$=' +str(be.data['user']['temp_grid'][te])+\
		 ' eV  ne=' + "%0.*e"%(2,be.data['user']['dens_grid'][ne]) + ' cm$^{-3}$',size=10)
   ax1.set_xlabel('Wavelength (nm)')
   ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')


.. figure:: be0_pec_0_1000.png
   :scale: 50 %
   :alt: Be I pecs 0-1000 nm


   
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

   upper_ind = 7 #ninth excited state
   lower_ind = 0  #ground state

   pec_ind = np.where( (be.data['processed']['pec_levels'][:,0] == upper_ind) &\
		       (be.data['processed']['pec_levels'][:,1] == lower_ind))[0]

   #plot the temeprature dependence of the chosen pec at first density in the grid
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.93,left=0.105,right=0.965)
   ax1.set_title('Temperature dependence of line ' +\
		 str(be.data['processed']['wave_air'][pec_ind]) +' nm',size=10)
   ax1.plot(be.data['user']['temp_grid'],be.data['processed']['pecs'][pec_ind[0],met,:,ne])
   ax1.set_xlabel('Temperature (eV)')
   ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')

   #plot the density dependence of the chosen pec at first density in the grid
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.93,left=0.125,right=0.965)
   ax1.set_title('Density dependence of line ' +\
		 str(be.data['processed']['wave_air'][pec_ind]) +' nm',size=10)
   ax1.plot(be.data['user']['dens_grid'],be.data['processed']['pecs'][pec_ind[0],met,te,:])
   ax1.set_xlabel('Density (cm$^{-3}$)')
   ax1.set_ylabel('PEC (ph cm$^{-3}$ s$^{-1}$)')


If the wavelength of a line of interest is known, the index can be found by looking at the
wavelength array.
The indices of all pecs that fall within the upper and lower bound of the 'where' statement are
returned. PECs can generally be distinguished by the actual value, large lines that are of interest
have much large PEC values, this can allow 


.. figure:: be0_pec_temp.png
   :scale: 50 %
   :alt: Be I temperature



.. figure:: be0_pec_dens.png
   :scale: 50 %
   :alt: Be I density



.. code-block:: python
   :linenos:

   #want to find the index of Be I line at 351.55
   pec_ind = np.where( (be.data['processed']['wave_air'] <352) &\
		       (be.data['processed']['wave_air'] >351))
   print('Wavelength from file ' + str(be.data['processed']['wave_air'][pec_ind[0]]))
   #Wavelength from file [351.55028742]
   print('PEC upper and lower levels '+ str(be.data['processed']['pec_levels'][pec_ind[0]]))
   #PEC upper and lower levels [[25  2]]
   

Generalized radiative coefficients (GCRs)
-----------------------------------------

The generalized collsional radiative coefficients are calculated by ColRadPy as well.
A description of these can be found in (Summers 2006), (Johnson 2019).
GCR coefficients are often used as inputs to plasma transport codes.
GCR coefficients are also used as inputs to ionization balance calculations which will be discussed
later. This allows for different ionization stages to be linked.

For example, the total ionization from one charge state to the other is defined as the SCD.
The total recombination from a charge state to the charge state of interest is defined as the ACD.
This gives the rate of population transfer from one ionization state to a lower ionization state.
The situation for systems with metastable states requires that the effective ionization and
recombination rates be metastable resolved.
In addition, it requires metastable cross coupling coefficients known as QCD and XCD coefficients.

Generally it is of interest to look at how the GCR coefficients change with some parameter such
as temperature. Plots are shown below of the different GCRs.


A physical description of the GCRs can be helpful in interpreting the meaning behind
them.


Metastable Cross Coupling Coefficient (QCD)
---------------------------------------------
The QCD coefficient represents the transfer of population from one metastable state to another within
the ionization state of interest and includes both direct population transfer between
metastable states as well as the transfer via an intermediate excited state.

GCR Ionization Coefficient (SCD)
-----------------------------------------
The total ionization from one charge state to the other is defined as the SCD.


GCR Recombination Coefficient (ACD)
----------------------------------------
The total recombination from a charge state to the charge state of interest is defined as the ACD.


Metastable Parent Cross Coupling Coefficient (XCD)
-------------------------------------------------------
The XCD coefficient represents the transfer of population between metastable states from
the ionization stage just above the stage of interest. Populations in the upper ionization
stage can recombine into the ionization state of interest from one metastable, redistribute
through all the states and then ionize back into a different metastable state of the upper
ionization state.


GCR Examples
---------------
For this example we will look at Be II this is soley because Be III has two metastable states.
This means that the XCD will have non-zero values. Remeber the call from before for Be I.

.. code-block:: python
   :linenos:

   import sys
   sys.path.append('../')
   from colradpy_class import colradpy
   import numpy as np

   fil = 'cpb03_ls#be1.dat' #adf04 file
   temperature_arr = np.linspace(1,100,20) #eV
   metastable_levels = np.array([0,1])   #ground and level 1 chosen to be metastable
   density_arr =     np.array([1.e13,8.e13,4.e14]) # cm-3
   beii = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
		 use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True)
   beii.solve_cr()

.. code-block:: python
   :linenos:
      
   #plotting the QCD
   import matplotlib.pyplot as plt
   plt.ion
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.125,right=0.965)
   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['qcd'][0,1,:,0]*1e5,
	    label = 'metastable cross coupling coefficient 1->2')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['qcd'][1,0,:,0]*1e5,
	    label = 'metastable cross coupling coefficient 2->1')
   ax1.legend()
   ax1.set_title('QCD plot')
   ax1.set_xlabel('Temperature (eV)')
   ax1.set_ylabel('QCD * 10$^5$ (cm$^{-3}$ s$^{-1}$)')


.. figure:: be1_qcd.png
   :scale: 50 %
   :alt: Be II QCD

	 
.. code-block:: python
   :linenos:
      
   #plotting the SCD
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.125,right=0.965)
   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['scd'][0,0,:,0],
	    label = 'metastable cross coupling coefficient 1->1+')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['scd'][0,1,:,0],
	    label = 'metastable cross coupling coefficient 1->2+')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['scd'][1,0,:,0],
	    label = 'metastable cross coupling coefficient 2->1+')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['scd'][1,1,:,0],
	    label = 'metastable cross coupling coefficient 2->2+')

   ax1.legend(fontsize='x-small',loc='best')
   ax1.set_title('SCD plot')
   ax1.set_xlabel('Temperature (eV)')
   ax1.set_ylabel('SCD (ion cm$^{-3}$ s$^{-1}$)')


.. figure:: be1_scd.png
   :scale: 50 %
   :alt: Be II SCD


.. code-block:: python
   :linenos:

   #plotting the ACD
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.075,right=0.965)
   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['acd'][0,0,:,0],
	    label = 'metastable cross coupling coefficient 1+->1')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['acd'][0,1,:,0],
	    label = 'metastable cross coupling coefficient 2+->1')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['acd'][1,0,:,0],
	    label = 'metastable cross coupling coefficient 1+->2')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['acd'][1,1,:,0],
	    label = 'metastable cross coupling coefficient 2+->2')

   ax1.legend(fontsize='x-small',loc='best')
   ax1.set_title('ACD plot')
   ax1.set_xlabel('Temperature (eV)')
   ax1.set_ylabel('ACD (rec cm$^{-3}$ s$^{-1}$)')



.. figure:: be1_acd.png
   :scale: 50 %
   :alt: Be II ACD

   
.. code-block:: python
   :linenos:

   #plotting the XCD
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.12,right=0.965)
   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['xcd'][0,1,:,0],
	    label = 'metastable cross coupling coefficient 1+->2+')

   ax1.plot(beii.data['user']['temp_grid'],
	    beii.data['processed']['scd'][1,0,:,0],
	    label = 'metastable cross coupling coefficient 2+->1+')
   ax1.legend(fontsize='x-small',loc='best')
   ax1.set_title('XCD plot')
   ax1.set_xlabel('Temperature (eV)')
   ax1.set_ylabel('XCD (cm$^{-3}$ s$^{-1}$)')



.. figure:: be1_xcd.png
   :scale: 50 %
   :alt: Be II XCD




	 
Determining Populating Mechanisms
---------------------------------
One feature unique to ColRadPy is the ability to determine the populating mechanism of levels.
This allows one to see which levels in a calculation are important to modeling the spectral lines of interest.
This allows those that generate the atomic data to know which transitions are required to accurately
model spectral lines. With this new analysis technique, it is possible to identify transitions that are
the most important and allow for complex systems such as high-Z near neutral systems to be simplified.



ColRadPy also allows the user to determine which intermediate levels populate a level of interest withThis is don if the summation is not carried out from the calculation of the QCD.
This allows one to see which levels in a calculation are important to modeling the spectral lines of interest. 




.. code-block:: python
   :linenos:
      
   #plotting the populating levels
   plt.figure()
   plt.figure();plt.plot(be.data['processed']['pop_lvl'][0,:,0,0,0]/\
                         np.sum(be.data['processed']['pop_lvl'][0,:,0,0,0]))

   plt.figure();plt.plot(be.data['processed']['pop_lvl'][0,:,0,10,0]/\
                         np.sum(be.data['processed']['pop_lvl'][0,:,0,10,0]))

   plt.figure();plt.plot(be.data['processed']['pop_lvl'][0,:,0,-1,0]/\
                         np.sum(be.data['processed']['pop_lvl'][0,:,0,-1,0]))
   
   plt.legend()
   plt.xlabel('Level number (#)')
   plt.ylabel('Populating fraction (-)')

   #plotting the populating fraction from the ground versus temperature
   plt.figure()
   plt.plot(be.data['user']['temp_grid'],
             be.data['processed']['pop_lvl'][10,0,0,:,0]/\
	     np.sum(be.data['processed']['pop_lvl'][10,:,0,:,0],axis=0))
	     
   plt.xlabel('Temperature (eV)')
   plt.ylabel('Populating fraction from ground (-)')



.. figure:: be0_pop_lvl.png
   :scale: 50 %
   :alt: Be I populating levels

   This shows that as temperature increase other excited levels contributed more and more
   to the first excited state

	 
.. figure:: be0_ground_contribution.png
   :scale: 50 %
   :alt: Be I ground contriubtion

   This shows that as the temperature increases the ground tributes less to the total population
   of level 1.
   

   

Advanced functionality
=======================

Time dependent CR modeling
-----------------------------


ColRadPy is also capable of solving the full collisional radiative matrix time-dependently.
This can be important for systems where there is significant population in
many excited states or where ultra fast timescales need to be considered.
Instead of the quasi-static approximation used in Equation 4 where excited states are assumed to
have no population change, the matrix is solved as a system of ordinary differential equations n (t) = An(t).
This method used to solve the system of equations was adapted from R. LeVeque.

Case in which with and without a source term can be considered in ColRadPy.
The case without a source term can used in a system like a linear machine with views that are
transverse to the direction of motion of the particles.

A source term can be used when the line of sight includes a source of particles.
The source term could also be used to model the pumping of specific levels with LIF.


.. code-block:: python
   :linenos:
      
   import sys
   sys.path.append('../')
   from colradpy_class import colradpy
   import numpy as np
   import matplotlib.pyplot as plt

   #Time dependent CR modeling
   td_t = np.geomspace(1.e-5,.1,1000)
   td_n0 = np.zeros(30)
   td_n0[0] = 1.

   fil = 'cpb03_ls#be0.dat' #adf04 file
   temperature_arr = np.array([10]) #eV
   metastable_levels = np.array([0])   #metastable level, just ground chosen here
   density_arr =     np.array([1.e9]) # cm-3
   be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
		 use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True,
		 td_t=td_t,td_n0=td_n0,td_source=td_s)
   be.solve_cr()
   be.solve_time_dependent()

   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.1,right=0.965)
   plt.plot(be.data['user']['td_t'],
	    be.data['processed']['td']['td_pop'][0,:,0,0],
	    label='Ground')
   plt.plot(be.data['user']['td_t'],
	    be.data['processed']['td']['td_pop'][1,:,0,0],
	    label='level 1')
   plt.plot(be.data['user']['td_t'],
	    be.data['processed']['td']['td_pop'][-1,:,0,0],
	    label='ion')
   ax1.legend(fontsize='x-small',loc='best')
   ax1.set_title('Time dependent solution of CR Be I no source term')
   ax1.set_xlabel('Time (s)')
   ax1.set_ylabel('Population (-)')



.. figure:: be0_time_dep_no_source.png
   :scale: 50 %
   :alt: Be I time dependence no source

   This time dependent collisional radiative model shows the time history for all Be I levels and
   the ground sate of Be II. This is the non-quasistatic solution, for a light system like Be the
   which only has one metastable the quasistatic approximation and non-quastatic solutions will
   give similiar results it is only for heavy species such as Mo and W where the quasistatic
   approximation starts to break down that this non-quasistatic solution is required.
   


.. code-block:: python
   :linenos:

   td_t = np.geomspace(1.e-5,1,1000)
   td_n0 = np.zeros(30)
   td_n0[0] = 1.
   td_s = np.zeros(30)
   td_s[0] = 1.
   fil = 'cpb03_ls#be0.dat' #adf04 file
   temperature_arr = np.array([10]) #eV
   metastable_levels = np.array([0])   #metastable level, just ground chosen here
   density_arr =     np.array([1.e8]) # cm-3
   be = colradpy(fil,metastable_levels,temperature_arr,density_arr,use_recombination=True,
		 use_recombination_three_body = True,use_ionization=True,suppliment_with_ecip=True,
		 td_t=td_t,td_n0=td_n0,td_source=td_s)

   be.solve_cr()
   be.solve_time_dependent()

   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.92,left=0.115,right=0.965)
   plt.plot(be.data['user']['td_t'],
	    be.data['processed']['td']['td_pop'][0,:,0,0],
	    label='Ground')
   plt.plot(be.data['user']['td_t'],
	    be.data['processed']['td']['td_pop'][1,:,0,0],
	    label='level 1')
   plt.plot(be.data['user']['td_t'],
	    be.data['processed']['td']['td_pop'][-1,:,0,0],
	    label='ion')
   ax1.legend(fontsize='x-small',loc='best')
   ax1.set_title('Time dependent solution of CR Be I with source term')
   ax1.set_xlabel('Time (s)')
   ax1.set_ylabel('Population (-)')
   

.. figure:: be0_time_dep_source.png
   :scale: 50 %
   :alt: Be I time dependence with source

   Time dependent solution with a constant source term of particles in the ground state.
   This could be used to model spectra where there is a constant erosion term from the
   wall. This could also be use to model level pumping in LIF systems.



Split LS resolved data to LSJ 
------------------------------

ColRadPy is able to split PECs from term resolved (LS) into level resolved (LSJ) values.
This currently does put PECs at the NIST wavelengths, a user must do this manually for now.
In the future this will be done automatically using the NIST database.



.. code-block:: python
   :linenos:


   import sys
   sys.path.append('../')
   from colradpy_class import *
   import numpy as np

   he = colradpy('./mom97_ls#he1.dat',[0],np.array([20]),np.array([1.e13]),use_recombination=False,
                  use_recombination_three_body = False,use_ionization=True)

   he.solve_cr()
   he.split_pec_multiplet()

   wave_8_3 = np.array([468.5376849,468.5757974,468.5704380])
   ind_8_3 = np.where( (he.data['processed']['pec_levels'][:,0] == 8) & \
                       (he.data['processed']['pec_levels'][:,1] == 3))[0]

   wave_6_5 = np.array([468.5407225,468.5568006])
   ind_6_5 = np.where( (he.data['processed']['pec_levels'][:,0] == 6) & \
                      (he.data['processed']['pec_levels'][:,1] == 5))[0]

   wave_7_3 = np.array([468.5524404,468.5905553])
   ind_7_3 = np.where( (he.data['processed']['pec_levels'][:,0] == 7) & \
                       (he.data['processed']['pec_levels'][:,1] == 3))[0]

   wave_9_4 = np.array([468.5703849, 468.5830890, 468.5804092])
   ind_9_4 = np.where( (he.data['processed']['pec_levels'][:,0] == 9) & \
                       (he.data['processed']['pec_levels'][:,1] == 4))[0]

   wave_6_4 = np.array([ 468.5917884, 468.5757080, 468.5884123])
   ind_6_4 = np.where( (he.data['processed']['pec_levels'][:,0] == 6) & \
                       (he.data['processed']['pec_levels'][:,1] == 4))[0]


   wave_468 = np.hstack((wave_8_3,wave_6_5,wave_7_3,wave_9_4,wave_6_4))
   pecs_468 = np.vstack((he.data['processed']['split']['pecs'][ind_8_3[0]],
			 he.data['processed']['split']['pecs'][ind_6_5[0]],
			 he.data['processed']['split']['pecs'][ind_7_3[0]],
			 he.data['processed']['split']['pecs'][ind_9_4[0]],
			 he.data['processed']['split']['pecs'][ind_6_4[0]]))[np.argsort(wave_468)]
   wave_468 = wave_468[np.argsort(wave_468)]



   plt.figure()
   plt.vlines(wave_468,np.zeros_like(wave_468),pecs_468[:,0,0,0])



Ionization Balance
====================
An ionization balance can be used to get the relative abundances of charge states in a given species.
The relative populations of charge states are solved using the GCR coefficient that are calculated
from the CR set of equations. A matrix similiar to the CR matrix is assembled using the GCR rate coefficients.
The QCD rates transfer population between metastable states in one ionization stage.
SCD is the ionization from one charge stage to the next.
ACD is the recombination from one charge stage to the previous stage and the XCD is population transfer between metastable
states through the next charge state.

ColRadPy is capable of preforming time dependent as well as time independent ionization balance calculations.
The values for time independent ionization balance are solved by looking at the sencond longested lived eigenvalue of the system.
The equations are then solved at eight times this value ensure that equilbrium of the system has been reached.



An example of the ionization balance code is run for Be from the example 'example/ion_bal.py'. First for a time dependent case and then for a time independent case.
In the plot of the time dependent abundances are shown as the solid lines and the time independent limits are shown as the dashed lines.



.. code-block:: python
   :linenos:

   import sys
   sys.path.append('../')
   from colradpy_class import colradpy
   import numpy as np
   import matplotlib.pyplot as plt
   from ionization_balance_class import ionization_balance

   #the adf04 files
   fils = np.array(['cpb03_ls#be0.dat','cpb03_ls#be1.dat','be2_adf04','be3_adf04'])
   temp = np.linspace(5,100,5) #temp grid
   dens = np.array([1.e11,1.e14]) #density grid
   metas = [np.array([0,1]),np.array([0]),np.array([0,1]),np.array([0])]#number of metastable
			      #this should match the metastables at the top of the adf04 file
			      #this information is used to calculate the QCD values
			      #without it only the SCD, ACD and XCD for a species will be calculated

   time = np.linspace(0,.01,1.e4)

   ion = ionization_balance(fils, metas, temp, dens, keep_species_data = True)
   ion.populate_ion_matrix()

   ion.solve_no_source(np.array([1,0,0,0,0,0,0]),time)

   ion.solve_time_independent()

   plt.ion
   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.99,left=0.11,right=0.99)
   ax1.plot(time*1e3,ion.data['pops'][0,:,1,1],label='be0, met0',color='b')
   ax1.hlines(ion.data['pops_ss'][0,0,1,1],0,10,color='b',linestyle=':')

   ax1.plot(time*1e3,ion.data['pops'][1,:,1,1],label='be0, met1',color='g')
   ax1.hlines(ion.data['pops_ss'][1,0,1,1],0,10,color='g',linestyle=':')

   ax1.plot(time*1e3,ion.data['pops'][2,:,1,1],label='be1, met0',color='r')
   ax1.hlines(ion.data['pops_ss'][2,0,1,1],0,10,color='r',linestyle=':')

   ax1.plot(time*1e3,ion.data['pops'][3,:,1,1],label='be2, met0',color='c')
   ax1.hlines(ion.data['pops_ss'][3,0,1,1],0,10,color='c',linestyle=':')

   ax1.plot(time*1e3,ion.data['pops'][4,:,1,1],label='be2, met1',color='m')
   ax1.hlines(ion.data['pops_ss'][4,0,1,1],0,10,color='m',linestyle=':')

   ax1.plot(time*1e3,ion.data['pops'][5,:,1,1],label='be3, met0',color='y')
   ax1.hlines(ion.data['pops_ss'][5,0,1,1],0,10,color='y',linestyle=':')

   ax1.plot(time*1e3,ion.data['pops'][6,:,1,1],label='be4',color='k')
   ax1.hlines(ion.data['pops_ss'][6,0,1,1],0,10,color='k',linestyle=':')

   ax1.legend(fontsize='x-small')

   ax1.set_xlabel('Time (ms)')
   ax1.set_ylabel('Fractional Abundance (-)')





.. figure:: be_ion_bal_time.png   
   :scale: 50 %
   :alt: Be time dependent ionization balance





The time independent solution as a function of electron temperature is shown below.



.. code-block:: python
   :linenos:

   temp = np.linspace(2,100,200) #temp grid
   ion = ionization_balance(fils, metas, temp, dens, keep_species_data = False)
   ion.populate_ion_matrix()
   ion.solve_time_independent()


   fig, ax1 = plt.subplots(1,1,figsize=(16/3.,9/3.),dpi=300)
   fig.subplots_adjust(bottom=0.15,top=0.99,left=0.11,right=0.99)

   ax1.plot(temp,ion.data['pops_ss'][0,0,:,1],label='be0, met0',color='b')

   ax1.plot(temp,ion.data['pops_ss'][1,0,:,1],label='be0, met1',color='g')

   ax1.plot(temp,ion.data['pops_ss'][2,0,:,1],label='be1, met0',color='r')

   ax1.plot(temp,ion.data['pops_ss'][3,0,:,1],label='be2, met0',color='c')

   ax1.plot(temp,ion.data['pops_ss'][4,0,:,1],label='be2, met1',color='m')

   ax1.plot(temp,ion.data['pops_ss'][5,0,:,1],label='be3, met0',color='y')

   ax1.plot(temp,ion.data['pops_ss'][6,0,:,1],label='be4',color='k')

   ax1.legend(fontsize='x-small')

   ax1.set_xlabel('Temperature (eV)')
   ax1.set_ylabel('Fractional Abundance (-)')



.. figure:: be_ion_bal_ind_time.png   
   :scale: 50 %
   :alt: Be time dependent ionization balance


	 

	    
   
Error bar analysis from atomic data
-----------------------------------

   
