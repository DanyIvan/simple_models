import sys
sys.path.append('../')

from goldblatt import SteadyStateModel, DynamicModel
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

gregory_data_file= '../../atmos-python/experiments_data/gregory_case1_20200322_prliminar_data.csv'
gregory_data = pd.read_csv(gregory_data_file)
gregory_data.O2_flow = gregory_data.O2_flow * 267
gregory_data.O2_fixedmr = gregory_data.O2_fixedmr / 0.21
gregory_data = gregory_data[['O2_flow', 'O2_fixedmr', 'experiment', 'converged']]
gregory_data.experiment = 'Gregory ' + gregory_data.experiment.str[28:33] 
data1 = gregory_data[gregory_data.converged == 1]
data2 = gregory_data[gregory_data.converged == 0]
data2.experiment = 'Not converged'

goldblatt_model = SteadyStateModel()
goldblatt_model.find_steady_states()
goldblatt_O2 = goldblatt_model.steady_states['primary_prod']['oxygen']
goldblatt_data = pd.DataFrame({'O2_flow': goldblatt_model.primary_prod_list,
     'O2_fixedmr': goldblatt_O2, 'experiment': 'Goldblatt r=3e10'})

goldblatt_model1 = SteadyStateModel(r=7.5e10)
goldblatt_model1.find_steady_states()
goldblatt_O21 = goldblatt_model1.steady_states['primary_prod']['oxygen']
goldblatt_data1 = pd.DataFrame({'O2_flow': goldblatt_model1.primary_prod_list,
     'O2_fixedmr': goldblatt_O21, 'experiment': 'Goldblatt r=7.5e10'})

goldblatt_model2 = SteadyStateModel(r=5e11)
goldblatt_model2.find_steady_states()
goldblatt_O22 = goldblatt_model2.steady_states['primary_prod']['oxygen']
goldblatt_data2 = pd.DataFrame({'O2_flow': goldblatt_model2.primary_prod_list,
     'O2_fixedmr': goldblatt_O22, 'experiment': 'Goldblatt r=5e11'})

goldblatt_model3 = SteadyStateModel(r=1e9)
goldblatt_model3.find_steady_states()
goldblatt_O23 = goldblatt_model3.steady_states['primary_prod']['oxygen']
goldblatt_data3 = pd.DataFrame({'O2_flow': goldblatt_model3.primary_prod_list,
     'O2_fixedmr': goldblatt_O23, 'experiment': 'Goldblatt r=1e9'})


data_gregory = pd.concat([data1, data2])
data_goldblat = pd.concat([goldblatt_data3, goldblatt_data, goldblatt_data1, goldblatt_data2])
fig = plt.figure(figsize=(9, 6)) 

sns.scatterplot(data=data_goldblat, x='O2_flow', y='O2_fixedmr', style='experiment',
    legend='full',alpha=0.85, palette=['w'])
sns.scatterplot(data=data_gregory, x='O2_flow', y='O2_fixedmr', style='experiment',
    legend='full',alpha=0.85, palette=['b'])

plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-10, 1e2])
plt.xlim([1e12, 1e17])
plt.ylabel(r'$O_2$ PAL', fontsize=15)
plt.xlabel(r'mol  $O_2$ / year', fontsize=15)


fig.savefig('GoldblattVSGregory.png', bbox_inches='tight')