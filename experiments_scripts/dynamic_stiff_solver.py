import sys
sys.path.append('..')
from goldblatt import DynamicSteadyStates
from utils.objects import save_object
from datetime import datetime

date = datetime.today().strftime('%Y%m%d_%H%M%S')
NAME = 'dynamic_stiff_solver_' + date
folder = '../experiments_data/'


states = DynamicSteadyStates(r=3e10, beta=0.01, dt=0.1, total_time=1e5,
    save_every=100)
states.find_steady_states()
save_object(states, folder + NAME + '.obj')

