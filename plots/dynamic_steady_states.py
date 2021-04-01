import sys
sys.path.append('../')

from goldblatt import SteadyStateModel, DynamicModel
from utils.conversions import mol_to_O2_PAL
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# vary_beta = []
# for beta in [0.1, 0.01, 0.001, 0.0001]:
#     model = DynamicModel(pp=3.75e15, r=3e10, dt=0.01, beta=beta, total_time=1e5) 
#     model.find_steady_states('pp')
#     O2 = mol_to_O2_PAL(np.array(model.steady_states_pp['O']))
#     M = mol_to_O2_PAL(np.array(model.steady_states_pp['M']))
#     data_model = pd.DataFrame({'O2': O2, 'M': M, 'pp': model.pp_list, 'beta': "$%s$" % beta})
#     vary_beta.append(data_model)
# vary_beta = pd.concat(vary_beta)

# fig, axs = plt.subplots(2, 1, figsize=(12, 8))
# sns.scatterplot(data=vary_beta,x='pp', y='O2', hue='beta', ax=axs[0], legend='full') 
# axs[0].set_xlabel('Primary Productivity mol O2 / year')
# axs[0].set_ylabel('O2 PAL')
# axs[0].set_xscale('log')
# axs[0].set_yscale('log')

# sns.scatterplot(data=vary_beta,x='pp', y='M', hue='beta', ax=axs[1], legend='full') 
# axs[1].set_xlabel('Primary Productivity mol O2 / year')
# axs[1].set_ylabel('Methane O2 PAL eq')
# axs[1].set_xscale('log')
# axs[1].set_yscale('log')
# fig.savefig('DynamicVaryBeta.png', bbox_inches='tight')


#----------------------------------------------------------------------------
# vary_r = []
# for r in [1e10, 1e11, 1e12, 1e13]:
#     model = DynamicModel(pp=3.75e15, r=r, dt=0.01, beta=0.01, total_time=1e4) 
#     model.find_steady_states('pp')
#     O2 = mol_to_O2_PAL(np.array(model.steady_states_pp['O']))
#     M = mol_to_O2_PAL(np.array(model.steady_states_pp['M']))
#     data_model = pd.DataFrame({'O2': O2, 'M': M, 'pp': model.pp_list, 'r':  "$%s$" % r})
#     vary_r.append(data_model)
# vary_r = pd.concat(vary_r)

# fig, axs = plt.subplots(2, 1, figsize=(12, 8))
# sns.scatterplot(data=vary_r,x='pp', y='O2', hue='r', ax=axs[0], legend='full') 
# axs[0].set_xlabel('Primary Productivity mol O2 / year')
# axs[0].set_ylabel('O2 PAL')
# axs[0].set_xscale('log')
# axs[0].set_yscale('log')

# sns.scatterplot(data=vary_r,x='pp', y='M', hue='r', ax=axs[1], legend='full') 
# axs[1].set_xlabel('Primary Productivity mol O2 / year')
# axs[1].set_ylabel('Methane O2 PAL eq')
# axs[1].set_xscale('log')
# axs[1].set_yscale('log')
# fig.savefig('DynamicVaryR.png', bbox_inches='tight')

#----------------------------------------------------------------------------
vary_t = []
for t in [1e2, 1e3, 1e4, 1e5]:
    model = DynamicModel(pp=3.375e15, r=3e10, dt=0.01, beta=0.01, total_time=t) 
    model.find_steady_states('pp')
    O2 = mol_to_O2_PAL(np.array(model.steady_states_pp['O']))
    M = mol_to_O2_PAL(np.array(model.steady_states_pp['M']))
    data_model = pd.DataFrame({'O2': O2, 'M': M, 'pp': model.pp_list, 't':  "$%s$" % t})
    vary_t.append(data_model)
vary_t = pd.concat(vary_t)

fig, axs = plt.subplots(2, 1, figsize=(12, 8))
sns.scatterplot(data=vary_t,x='pp', y='O2', hue='t', ax=axs[0], legend='full') 
axs[0].set_xlabel('Primary Productivity mol O2 / year')
axs[0].set_ylabel('O2 PAL')
axs[0].set_xscale('log')
axs[0].set_yscale('log')

sns.scatterplot(data=vary_t,x='pp', y='M', hue='t', ax=axs[1], legend='full') 
axs[1].set_xlabel('Primary Productivity mol O2 / year')
axs[1].set_ylabel('Methane O2 PAL eq')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
fig.savefig('DynamicVaryTime.png', bbox_inches='tight')