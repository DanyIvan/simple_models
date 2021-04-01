import sys
sys.path.append('../')

from goldblatt import SteadyStateModel, DynamicModel
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

model1 = SteadyStateModel()
model1.find_steady_states()

# model2 = SteadyStateModel(r=7.5e10)
model2 = SteadyStateModel(pp=1e15, r=7.5e10)
model2.find_steady_states()

fig, axs = plt.subplots(2, 1, figsize=(12, 10))
plt.sca(axs[0])
plt.scatter(model1.primary_prod_list, model1.steady_states['primary_prod']\
    ['oxygen'], label='r= 3e10')
plt.scatter(model2.primary_prod_list, model2.steady_states['primary_prod']\
    ['oxygen'], label='r = 7.5e10')
plt.xlabel('Primary Productivity mol O2 / year')
plt.ylabel('O2 PAL')
plt.legend()
plt.xscale('log')
plt.yscale('log')


plt.sca(axs[1])
plt.scatter(model1.reductant_input_list,
    model1.steady_states['reductant_input']['oxygen'], label='N = 3.75e15')
plt.scatter(model2.reductant_input_list,
    model2.steady_states['reductant_input']['oxygen'], label='N = 1.8e15')
plt.xlabel('Reductant input mol O2')
plt.ylabel('O2 PAL')
plt.legend()
plt.xscale('log')
plt.yscale('log')
fig.savefig('GoldblattFigure.png', bbox_inches='tight')


fig1, axs = plt.subplots(2, 1, figsize=(12, 10))
plt.sca(axs[0])
plt.scatter(model1.primary_prod_list, model1.steady_states['primary_prod']\
    ['methane'], label='r= 3e10')
plt.scatter(model2.primary_prod_list, model2.steady_states['primary_prod']\
    ['methane'], label='r = 7.5e10')
plt.xlabel('Primary Productivity mol O2 / year')
plt.ylabel('O2 PAL')
plt.legend()
plt.xscale('log')
plt.yscale('log')


plt.sca(axs[1])
plt.scatter(model1.reductant_input_list,
    model1.steady_states['reductant_input']['methane'], label='N = 3.75e15')
plt.scatter(model2.reductant_input_list,
    model2.steady_states['reductant_input']['methane'], label='N = 1.8e15')
plt.xlabel('Reductant input mol O2')
plt.ylabel('O2 PAL')
plt.legend()
plt.xscale('log')
plt.yscale('log')
fig1.savefig('GoldblattFigureMethane.png', bbox_inches='tight')