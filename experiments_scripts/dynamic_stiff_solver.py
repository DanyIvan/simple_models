import sys
sys.path.append('..')
from goldblatt import DynamicSteadyStates
from utils.objects import save_object
from datetime import datetime
import numpy as np


folder = '../experiments_data/'



betas = [1e-2, 1e-3, 1e-4]
pp_lists = [np.logspace(13, 16, 50),np.logspace(16, 13, 50)]

for pp_list in pp_lists:
    for beta in betas:
        print('beta=' + format(beta, 'E'))
        states = DynamicSteadyStates(r=3e10, beta=beta, dt=1, total_time=3e5,
            save_every=100)
        states.find_steady_states(pp_list=pp_list)
        date = datetime.today().strftime('%Y%m%d_%H%M%S')
        NAME = 'dynamic_stiff_solver_' + date + '.obj'
        save_object(states, folder + NAME)

