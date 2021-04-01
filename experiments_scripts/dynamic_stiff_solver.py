import sys
sys.path.append('..')
from goldblatt import DynamicSteadyStates
from utils.objects import save_object
from datetime import datetime

date = datetime.today().strftime('%Y%m%d_%H%M%S')
NAME = 'dynamic_stiff_solver_' + date
folder = '../experiments_data/'

save_everys = [100, 100, 1000]
times = [1e4, 1e5, 1e6]
betas = [1e-1, 1e-2, 1e-3, 1e-4]
rs = [1e10, 1e11, 1e12]
for r in rs:
    for beta in betas:
        for idx, time in enumerate(times):
            print('r='+ format(r, 'E'))
            print('beta=' + format(beta, 'E'))
            print('time=' + format(time, 'E'))
            states = DynamicSteadyStates(r=r, beta=beta, dt=0.1, total_time=time,
                save_every=save_everys[idx])
            states.find_steady_states()
            save_object(states, folder + NAME + '.obj')

