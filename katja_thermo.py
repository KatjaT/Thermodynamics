# -*- coding: utf-8 -*-
"""
calculate thermodynamics for Katja
"""
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T

import csv
import numpy as np
import uncertainties.unumpy as unumpy  

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
        
def reaction2RI(reaction_list, fixed_conc=0.1):
    '''
        Calculates the reversibility index (RI) of a reaction.
        The RI represent the change in concentrations of metabolites
        (from equal reaction reactants) that will make the reaction reversible.
        That is, the higher RI is, the more irreversible the reaction.
        A convenient threshold for reversibility is RI>=1000, that is a change of
        1000% in metabolite concentrations is required in order to flip the
        reaction direction. 
        
        Arguments:
            List of cobra model reaction objects
        Returns:
            Array of RI values
    '''
    keq = reaction2Keq(reaction_list)
    
    sparse = map(lambda x: KeggReaction.parse_formula(x).sparse, reaction_list)
    
    N_P = np.zeros(len(sparse))
    N_S = np.zeros(len(sparse))    
    for i,s in enumerate(sparse):   
        N_P[i] = sum([v for v in s.itervalues() if v>0])
        N_S[i] = -sum([v for v in s.itervalues() if v<0])
    N = N_P + N_S
    Q_2prime = fixed_conc**(N_P-N_S)
    
    RI = ( keq*Q_2prime )**( 2.0/N )
    return RI
        
if __name__ == "__main__":

    reactions = csv.reader(open('CCMtbRxnsKEGG.txt', 'r'))
    names = []
    reaction_list = []
    for row in reactions:
        row = row[0].split("    ")
        names.append(row[0].replace("'", ''))
        reaction_list.append(row[1])
    
    dG0 = reaction2dG0(reaction_list)
    Keq = reaction2Keq(reaction_list)
    RI = reaction2RI(reaction_list)
    
    reversibility_index = dict(zip(names, RI))

    f = open('reversibility_index.csv','w')
    w = csv.writer(f)
    for k,v in reversibility_index.iteritems():
        w.writerow([k, v])
    f.close()