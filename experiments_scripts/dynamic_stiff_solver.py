import sys
sys.path.append('..')
from goldblatt import DynamicSteadyStates
from utils.objects import save_object
from datetime import datetime
import numpy as np


folder = '../experiments_data/'



betas = [1e-2, 1e-3, 1e-4]
rs = [3e10, 1e11, 1e12]

for r in rs:
    for beta in betas:
        print('r='+ format(r, 'E'))
        print('beta=' + format(beta, 'E'))
        states = DynamicSteadyStates(r=r, beta=beta, dt=1, total_time=1e6,
            save_every=1000)
        states.find_steady_states(pp_list=np.logspace(13, 16, 50))
        date = datetime.today().strftime('%Y%m%d_%H%M%S')
        NAME = 'dynamic_stiff_solver_' + date + '.obj'
        save_object(states, folder + NAME + '.obj')

