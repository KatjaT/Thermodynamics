# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:33:02 2015

@author: katja

calculate thermodynamics for Katja

Reversibility index with adding known concentrations
Modified from Dans initial script 15/04/15

"""
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T

import csv
import numpy as np
import uncertainties.unumpy as unumpy  
import pandas as pd

def reaction2dG0(reaction_list):
    '''
        Calculates the dG0 of a list of a reaction.
        Uses the component-contribution package (Noor et al) to estimate
        the standard Gibbs Free Energy of reactions based on 
        component contribution  approach and measured values (NIST and Alberty)
        
        Arguments:
            List of reaction strings
        Returns:
            Array of dG0 values and standard deviation of estimates
    '''
    cc = ComponentContribution.init()
    
    Kmodel = KeggModel.from_formulas(reaction_list)

    Kmodel.add_thermo(cc)
    dG0_prime, dG0_std = Kmodel.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
    dG0_prime = np.array(map(lambda x: x[0,0], dG0_prime))
    
    dG0_prime = unumpy.uarray(dG0_prime, np.diag(dG0_std))          
    return dG0_prime

def reaction2Keq(reaction_list):
    '''
        Calculates the equilibrium constants of a reaction, using dG0.
        
        Arguments:
            List of cobra model reaction objects
        Returns:
            Array of K-equilibrium values
    '''
    dG0_prime = reaction2dG0(reaction_list)
    Keq = unumpy.exp( -dG0_prime / (R*default_T) )

    return Keq
        
def reaction2RI(reaction_list, metabolite_list, fixed_conc=0.1):
    '''
        Calculates the reversibility index (RI) of a reaction.
        The RI represent the change in concentrations of metabolites
        (from equal reaction reactants) that will make the reaction reversible.
        That is, the higher RI is, the more irreversible the reaction.
        A convenient threshold for reversibility is RI>=1000, that is a change of
        1000% in metabolite concentrations is required in order to flip the
        reaction direction. 
        
        ToDo
        Includes fixed meatbolite concentrations for the measured metabolites 
        
        Arguments:
            List of cobra model reaction objects
        Returns:
            Array of RI values
    '''
    keq = reaction2Keq(reaction_list)
    
    print metabolite_list
    sparse = map(lambda x: KeggReaction.parse_formula(x).sparse, reaction_list)

    # map concentrations to reactions
    N_u = np.array(map(lambda x: len(set(x.keys())-set(metabolite_list.keys())), sparse))

    concs = [{k:metabolite_list[k] if k in metabolite_list.keys() 
                                    else fixed_conc 
                                    for k in l.iterkeys()} 
                                    for l in sparse] #new array with concentrations instead of stoichiometric values
    
    Q =  np.zeros(len(concs))
    for i, (c, s) in enumerate(zip(concs, sparse)):
        Q[i] = np.prod(np.array([v**s[k] for k, v in c.iteritems()]))
    # ToDo: if all concentration are known, take only Keq/Q
    RI = ( keq/Q )**( 2.0/N_u )
    return RI,sparse, Q, N_u
    
    
        
if __name__ == "__main__":
    # Read reaction list
    reactions = csv.reader(open('CCMtbRxnsKEGG.txt', 'r'))
    names = []
    reaction_list = []
    for row in reactions:
        row = row[0].split("    ")
        names.append(row[0].replace("'", ''))
        reaction_list.append(row[1])
        
    # Read metabolite lists
    metabolite_list = pd.read_csv('metabolomics4.csv').set_index('Unnamed: 0')
    # for condition in ['Glc', 'Prop', 'Glu']:    
    metabolite_list = metabolite_list.loc['Glc'] * 1e-6 # conversion from uM to M
    # ToDo: remove NaNs
        
    # Calculate properties    
    dG0 = reaction2dG0(reaction_list)
    Keq = reaction2Keq(reaction_list)
    RI, sparse, Q,  N_u = reaction2RI(reaction_list,metabolite_list) 
    reversibility_index = dict(zip(names, RI))

    f = open('reversibility_index_concs.csv','w')
    w = csv.writer(f)
    for k,v in reversibility_index.iteritems():
        w.writerow([k, v]) #ksdjfk
    f.close()