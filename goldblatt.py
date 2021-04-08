import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from scipy.integrate import ode
import matplotlib.pyplot as plt
#from scipy.integrate import RK45
from utils.conversions import mol_to_O2_PAL
# from numba import njit, generated_jit
# from numba.typed import List


class Base():
    '''This class stores properties common to the SteadyStateModel and the 
    DynamicModel classes'''
    def __init__(self, *kwargs):
        # hydrogen escape coefficient
        self.s = 3.7e-5
        # rate of exposure of buried organic carbon from the crust 
        self.w = 6 * 10**(-9)
        #  inhibition ofaerobic respiration at the Pasteur point [mol]
        self.d_gamma = 1.36 * 10**19
        # fraction of OC decomposed by heterotrophic aerobic respiration
        self.gamma = lambda O : O / (O + self.d_gamma)
        # coefficient of methane used by methanotrophs[mol]
        self.d_delta = 2.73 * 10**17
        # fraction of methane produced used by methanotrophs
        self.delta = lambda O : O / (O + self.d_delta)
        # fraction of the organic carbonavailable to decomposers 
        self.omega = lambda O : (1 - self.gamma(O)) * (1 - self.delta(O))
        # coefficients of methane oxidation parametrization
        self.p = [0.0030084, -0.1655405, 3.2305351, -25.8343054,
            71.5397861]
        # methane oxidation parametrization
        self.Phi = lambda x : 10**(self.p[0] * np.log10(x)**4 + self.p[1] *\
            np.log10(x)**3 + self.p[2] * np.log10(x)**2 + self.p[3] *\
            np.log10(x) + self.p[4] )

class SteadyStateModel(Base):
    def __init__(self, **kwargs):
        super().__init__()
        self.oxygen_eq = lambda N, r: lambda O: self.omega(O) * N + ((self.omega(O) -\
            2) * r) - self.Phi(O) * ((r / self.s)**0.7)
        self.beta = kwargs.get('beta') or 0.01
        self.r = kwargs.get('r') or 3e10
        self.pp = kwargs.get('pp') or 3.75e15
        self.M = self.r / self.s
        self.C = self.beta * (self.pp +self.r) / self.w
        self.steady_states = { 
            'primary_prod': {
                'oxygen': [], 'methane': [], 'organic_carbon': []
                },
            'reductant_input': {
                'oxygen': [], 'methane': [], 'organic_carbon': []
                }
        }
        self.reductant_input_list = None 
        self.primary_prod_list = None

    def atmosphere_flux(self, primary_prod): 
        flux = [] 
        for pp in primary_prod: 
            flux.append(pp * (self.omega(pp) + self.beta*(1 - self.omega(pp)))) 
        return flux 
            
    
    def find_steady_states(self, N_points=1000, equation=1, atmosphere_flux=False):
        # alternative oxygen equation
        oxygen_eq = lambda M, C: lambda pp, r: lambda O:  self.omega(O) * pp -\
                 (1 - self.omega(O)) * r +\
                 self.beta * (1 - self.omega(O)) * (pp + r) -\
                 (1 - self.omega(O)) * self.w * C -\
                 self.s * M -\
                 self.Phi(O) * M**0.7
        # solve for primary productivity
        # solve upwards
        # ----------------------------------------------------------
        x0 = 1e11
        primary_prod_upwards  = np.logspace(11, 18, 1000)
        if atmosphere_flux:
            primary_prod_upwards = self.atmosphere_flux(primary_prod_upwards)
            print(primary_prod_upwards[-1])
        # get point of bistability
        # x0_list = []
        # for idx, primary_prod in enumerate(primary_prod_upwards):
        #     M = self.r / self.s
        #     C = self.beta * (primary_prod + self.r) / self.w
        #     eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
        #     solver = lambda N, r: fsolve(eq(N, r), x0)
        #     x0 = solver(primary_prod, self.r)[0]
        #     x0_list.append(x0)
        #     if idx > 2:
        #         if x0 == x0_list[idx -2]:
        #             x0_max = x0
        #             print(x0_max)
        #             break
        # first solve
        x0 = 1e11
        x0_list = []
        for idx, primary_prod in enumerate(primary_prod_upwards):
            M = self.r / self.s
            C = self.beta * (primary_prod + self.r) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0)
            x0 = solver(primary_prod, self.r)[0]
            x0_list.append(x0)
            # help the solver in the bistability jump
            x0 = x0 if x0 < 739943473002414.6 else 1e18
            # x0 = x0 if x0 < x0_max else 1e18
        # second solve
        oxygen_upwards, methane_upwards, OC_upwards = [], [], []
        for idx, primary_prod in enumerate(primary_prod_upwards):
            M = self.r / self.s
            C = self.beta * (primary_prod + self.r) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0_list[idx])
            oxygen_upwards.append(mol_to_O2_PAL(solver(primary_prod, self.r)[0]))
            methane_upwards.append(mol_to_O2_PAL(M))   
            OC_upwards.append(mol_to_O2_PAL(C))  
        
        # solve backwards
        # ----------------------------------------------------------
        x0 = 1e19
        primary_prod_backwards = np.logspace(18, 11, 1000)
        if atmosphere_flux:
            primary_prod_backwards = self.atmosphere_flux(primary_prod_backwards)
        # get point of bistability
        # x0_list = []
        # for idx, primary_prod in enumerate(primary_prod_backwards):
        #     M = self.r / self.s
        #     C = self.beta * (primary_prod + self.r) / self.w
        #     eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
        #     solver = lambda N, r: fsolve(eq(N, r), x0)
        #     x0_list.append(solver(primary_prod, self.r)[0])
        #     x0 = solver(primary_prod, self.r)[0]
        #     # help the solver in the bistability jump
        # x0_min = min(x0_list)

        x0 = 1e19
        x0_list = []
        for primary_prod in primary_prod_backwards:
            M = self.r / self.s
            C = self.beta * (primary_prod + self.r) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0)
            x0_list.append(solver(primary_prod, self.r)[0])
            x0 = solver(primary_prod, self.r)[0]
            # help the solver in the bistability jump
            # x0 = x0 if x0 > x0_min else 1e11
            x0 = x0 if x0 > 1.8e17 else 1e11
        # second solve
        oxygen_backwards, methane_backwards, OC_backwards = [], [], []
        for idx, primary_prod in enumerate(primary_prod_backwards):
            M = self.r / self.s
            C = self.beta * (primary_prod + self.r) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0_list[idx])
            oxygen_backwards.append(mol_to_O2_PAL(
                solver(primary_prod, self.r)[0]))
            methane_backwards.append(mol_to_O2_PAL(M))   
            OC_backwards.append(mol_to_O2_PAL(C))

        #save steady states
        self.steady_states['primary_prod']['oxygen'] = oxygen_upwards +\
             oxygen_backwards
        self.steady_states['primary_prod']['methane'] = methane_upwards +\
             methane_backwards  
        self.steady_states['primary_prod']['organic_carbon'] = OC_upwards +\
             OC_backwards  
        self.primary_prod_list = list(primary_prod_upwards) + list(primary_prod_backwards)  

        # solve for reductant input
        # solve upwards
        # ----------------------------------------------------------
        x0 = 1e19
        reduct_upwards  = np.logspace(8, 16, N_points)
         # first solve
        x0_list = []
        for reduct in reduct_upwards:
            M = reduct / self.s
            C = self.beta * (self.pp + reduct) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0)
            x0_list.append(solver(self.pp, reduct)[0])
            x0 = solver(self.pp, reduct)[0]
            # help the solver in the bistability jump
            x0 = x0 if x0 > 1.7e17 else 1e9
        # second solve
        oxygen_upwards, methane_upwards, OC_upwards = [], [], []
        for idx, reduct in enumerate(reduct_upwards):
            M = reduct / self.s
            C = self.beta * (self.pp + reduct) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0_list[idx])
            oxygen_upwards.append(mol_to_O2_PAL(solver(self.pp, reduct)[0]))
            methane_upwards.append(mol_to_O2_PAL(M))   
            OC_upwards.append(mol_to_O2_PAL(C))           
        # solve backwards
        # ----------------------------------------------------------
        x0 = 1e9
        reduct_backwards = np.logspace(16, 8, 1000)
        # first solve
        x0_list = []
        for reduct in reduct_backwards:
            M = reduct / self.s
            C = self.beta * (self.pp + reduct) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0)
            x0_list.append(solver(self.pp, reduct)[0])
            x0 = solver(self.pp, reduct)[0]
            self.steady_states['reductant_input']['oxygen'].append(
                 mol_to_O2_PAL(solver(self.pp, reduct)))
            # help the solver in the bistability jump
            x0 = x0 if x0 < 7.2e14 else 1e18
        # # second solve
        oxygen_backwards, methane_backwards, OC_backwards = [], [], []
        for idx, reduct in enumerate(reduct_backwards):
            M = reduct / self.s
            C = self.beta * (self.pp + reduct) / self.w
            eq = oxygen_eq(M, C) if equation == 1 else self.oxygen_eq
            solver = lambda N, r: fsolve(eq(N, r), x0_list[idx])
            oxygen_backwards.append(mol_to_O2_PAL(solver(self.pp, reduct)[0]))
            methane_backwards.append(mol_to_O2_PAL(M))   
            OC_backwards.append(mol_to_O2_PAL(C))
            
        #save steady states
        self.steady_states['reductant_input']['oxygen'] = oxygen_upwards +\
             oxygen_backwards
        self.steady_states['reductant_input']['methane'] = methane_upwards +\
             methane_backwards  
        self.steady_states['reductant_input']['organic_carbon'] = OC_upwards +\
             OC_backwards  
        self.reductant_input_list = list(reduct_upwards) + list(reduct_backwards)                 

    def plot(self):
        fig, axs = plt.subplots(2, 1, figsize=(4, 10))
        plt.sca(axs[0])
        plt.scatter(self.primary_prod_list, self.steady_states['primary_prod']\
            ['oxygen'], label='O')
        # plt.scatter(self.primary_prod_list, self.steady_states['primary_prod']\
        #     ['methane'], label='M')
        # plt.scatter(self.primary_prod_list, self.steady_states['primary_prod']\
        #     ['organic_carbon'], label='C')
        plt.xlabel('Primary Productivity mol O2')
        plt.ylabel('mol O2')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')


        plt.sca(axs[1])
        plt.scatter(self.reductant_input_list,
            self.steady_states['reductant_input']['oxygen'], label='O')
        # plt.scatter(self.reductant_input_list,
        #     self.steady_states['reductant_input']['methane'], label='M')
        # plt.scatter(self.reductant_input_list,
        #     self.steady_states['reductant_input']['organic_carbon'], label='C')
        plt.xlabel('reductant_input mol O2')
        plt.ylabel('mol O2')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

        
class DynamicModel(Base):
    def __init__(self, pp, r, **kwargs):
        super().__init__()
        #primary productivity
        self.pp = pp
        # reductabt input
        self.r = r
        # total time to run the model
        self.total_time = int(kwargs.get('total_time') or 1e7)
        # time step
        self.dt = kwargs.get('dt') or 1
        self.save_every = kwargs.get('save_every') or 100
        # number of steps
        self.N = int(self.total_time / self.dt / self.save_every)
        # fraction,of the available organic carbon in the surface ocean
        # buried in sediments,
        self.beta = kwargs.get('beta') or 0.4
        # methaane surface reservoir
        self.M = []
        # initial value of M
        self.M0 = kwargs.get('M0') or self.r / self.s
        self.M.append(self.M0)
        # oxygen surface reservoir
        self.O = []
        # initial value of O
        self.O0 = kwargs.get('O0') or 2 * self.M0
        self.O.append(self.O0)
        # organic matter in crust reservoir
        self.C = []
        # initial value of C
        self.C0 = kwargs.get('C0') or self.beta * (self.r + self.pp) / self.w
        self.C.append(self.C0)
        self.time = []
        self.time.append(0.0)


    @staticmethod
    # @njit
    def _run(total_time, dt, save_every,  M, O, C, time, beta, pp, r, w, s):
        d_gamma = 1.36 * 10**19
        gamma = lambda O : O / (O + d_gamma)
        d_delta = 2.73 * 10**17
        delta = lambda O : O / (O + d_delta)
        omega = lambda O : (1 - gamma(O)) * (1 - delta(O))
        p = [0.0030084, -0.1655405, 3.2305351, -25.8343054,
            71.5397861]
        phi = lambda x : 10**(p[0] * np.log10(x)**4 + p[1] *\
            np.log10(x)**3 + p[2] * np.log10(x)**2 + p[3] *\
            np.log10(x) + p[4])

        M_last, O_last, C_last = M[0], O[0], C[0]
        for i in range(total_time // dt -1):
            # print('Step: ' + str(i))
            # print('---------------------')
            # print('M_last:' + str(M_last))
            # print('O_last:' +  str(O_last))
            # print('C_last:' +  str(C_last))
            # print('omega(O_last):' +  str(omega(O_last)))
            # print('phi(O_last):' +  str(phi(O_last)))

            M_last = M_last + dt / 2 * (
                omega(O_last) * (1 - beta) * (pp + r) +\
                omega(O_last) * w * C_last -\
                2 * s * M_last -\
                phi(O_last) * M_last**0.7
            )
            O_last = O_last + dt * (
                omega(O_last) * pp -\
                (1 - omega(O_last)) * r +\
                beta * (1 - omega(O_last)) * (pp + r) -\
                (1 - omega(O_last)) * w * C_last -\
                s * M_last -\
                phi(O_last) * M_last**0.7
            )
            C_last = C_last + dt * (
                beta * (pp + r) - w * C_last
            )

            # k1_M = dt / 2 * (
            #     omega(O_last) * (1 - beta) * (pp + r) +\
            #     omega(O_last) * w * C_last -\
            #     2 * s * M_last -\
            #     phi(O_last) * M_last**0.7
            # )
            # k2_M = dt / 2 * (
            #     omega(O_last) * (1 - beta) * (pp + r) +\
            #     omega(O_last) * w * C_last -\
            #     2 * s * (M_last + k1_M / 2) -\
            #     phi(O_last) * (M_last + k1_M / 2)**0.7
            # )
            # k3_M = dt / 2 * (
            #     omega(O_last) * (1 - beta) * (pp + r) +\
            #     omega(O_last) * w * C_last -\
            #     2 * s * (M_last + k2_M / 2) -\
            #     phi(O_last) * (M_last + k2_M / 2)**0.7
            # )
            # k4_M = dt / 2 * (
            #     omega(O_last) * (1 - beta) * (pp + r) +\
            #     omega(O_last) * w * C_last -\
            #     2 * s * (M_last + k3_M) -\
            #     phi(O_last) * (M_last + k3_M)**0.7
            # )
            # M_next = M_last + (k1_M / 2 + k2_M + k3_M + k4_M / 2) / 3

            # k1_O = dt * (omega(O_last) * pp -\
            #     (1 - omega(O_last)) * r +\
            #     beta * (1 - omega(O_last)) * (pp + r) -\
            #     (1 - omega(O_last)) * w * C_last -\
            #     s * M_last -\
            #     phi(O_last) * M_last**0.7
            # )
            # k2_O = dt * (omega(O_last + k1_O / 2) * pp -\
            #     (1 - omega(O_last + k1_O / 2)) * r +\
            #     beta * (1 - omega(O_last + k1_O / 2)) * (pp + r) -\
            #     (1 - omega(O_last + k1_O / 2)) * w * C_last -\
            #     s * M_last -\
            #     phi(O_last + k1_O / 2) * M_last**0.7
            # )
            # k3_O = dt * (omega(O_last + k2_O / 2) * pp -\
            #     (1 - omega(O_last + k2_O / 2)) * r +\
            #     beta * (1 - omega(O_last + k2_O / 2)) * (pp + r) -\
            #     (1 - omega(O_last + k2_O / 2)) * w * C_last -\
            #     s * M_last -\
            #     phi(O_last + k2_O / 2) * M_last**0.7
            # )
            # k4_O = dt * (omega(O_last + k3_O) * pp -\
            #     (1 - omega(O_last + k3_O)) * r +\
            #     beta * (1 - omega(O_last + k3_O)) * (pp + r) -\
            #     (1 - omega(O_last + k3_O)) * w * C_last -\
            #     s * M_last -\
            #     phi(O_last + k3_O) * M_last**0.7
            # )
            # O_next = O_last + (k1_O / 2 + k2_O + k3_O + k4_O / 2) / 3

            # k1_C = dt * (beta * (pp + r) - w * C_last)
            # k2_C = dt * (beta * (pp + r) - w * (C_last + k1_C / 2))
            # k3_C = dt * (beta * (pp + r) - w * (C_last + k2_C / 2))
            # k4_C = dt * (beta * (pp + r) - w * (C_last + k3_C))
            # C_next = C_last + (k1_C / 2 + k2_C + k3_C + k4_C / 2) / 3

            #M_last, O_last, C_last = M_next, O_next, C_next
            if i%save_every == 0:
                M.append(M_last)
                O.append(O_last)
                C.append(C_last)
                time.append(dt * i)
        return M, O, C, time


    def run_rk5(self):
        M_prime = lambda M, O, C: \
                self.omega(O) * (1 - self.beta) * (self.pp + self.r) +\
                self.omega(O) * self.w * C -\
                2 * self.s * M -\
                self.Phi(O) * M**0.7
        O_prime = lambda M, O, C: self.omega(O) * self.pp -\
            (1 - self.omega(O)) * self.r +\
            self.beta * (1 - self.omega(O)) * (self.pp + self.r) -\
            (1 - self.omega(O)) * self.w * C -\
            self.s * M -\
            self.Phi(O) * M**0.7
        C_prime = lambda C: self.beta * (self.pp + self.r) - self.w * C
        y0 = np.array([self.M[0], self.O[0], self.C[0]])
        func = lambda t, y: np.array([
                                M_prime(y[0], y[1], y[2]),
                                O_prime(y[0], y[1], y[2]),
                                C_prime(y[2])                                
                                ])
        solver = RK45(func, 0, y0, self.total_time + 1, first_step=1)
        while solver.t < self.total_time:
            solver.step()
            print(solver.y)
            print(solver.t)
            self.M.append(solver.y[0])
            self.O.append(solver.y[1])
            self.C.append(solver.y[2])
            self.time.append(solver.t)

    def run_stiff(self):
        M_prime = lambda M, O, C: \
                self.omega(O) * (1 - self.beta) * (self.pp + self.r) +\
                self.omega(O) * self.w * C -\
                2 * self.s * M -\
                self.Phi(O) * M**0.7
        O_prime = lambda M, O, C: self.omega(O) * self.pp -\
            (1 - self.omega(O)) * self.r +\
            self.beta * (1 - self.omega(O)) * (self.pp + self.r) -\
            (1 - self.omega(O)) * self.w * C -\
            self.s * M -\
            self.Phi(O) * M**0.7
        C_prime = lambda C: self.beta * (self.pp + self.r) - self.w * C
        y0 = np.array([self.M[0], self.O[0], self.C[0]])
        func = lambda t, y: np.array([
                                M_prime(y[0], y[1], y[2]),
                                O_prime(y[0], y[1], y[2]),
                                C_prime(y[2])                                
                                ])
        solver = ode(func).set_integrator('zvode', method='bdf', first_step=1).\
            set_initial_value(y0, 0)
        while solver.successful() and  solver.t < self.total_time:
            sol = solver.integrate(solver.t + self.dt)
            if round(solver.t, 1) % self.save_every == 0:
                #print(solver.t)
                self.M.append(np.real(solver.y[0]))
                self.O.append(np.real(solver.y[1]))
                self.C.append(np.real(solver.y[2]))
                self.time.append(round(solver.t, 1))

    def run(self):
        M, O, C, time = List(), List(), List(), List()
        M.append(self.M[0])
        O.append(self.O[0])
        C.append(self.C[0])
        time.append(0.0)
        self.M, self.O, self.C, self.time = self._run(self.total_time, self.dt,
            self.save_every, M, O, C, time, self.beta, self.pp, self.r, self.w,
            self.s)
        
    def plot(self, n=10):
        plt.plot(self.time[0:-1:n], (np.array(self.O[0:-1:n])), label='O')
        plt.plot(self.time[0:-1:n], (np.array(self.M[0:-1:n])), label='M')
        plt.plot(self.time[0:-1:n], (np.array(self.C[0:-1:n])), label='C')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('years')
        plt.ylabel('mol O2')
        plt.show()
    
class DynamicSteadyStates:
    def __init__(self, r, beta, dt, total_time, save_every=100):
        self.r = r 
        self.beta = beta
        self.dt = dt
        self.total_time = total_time
        self.save_every = save_every
        self.pp_list = None
        self.data = None

    def find_steady_states(self, pp_list=np.logspace(11, 17, 100)):
        self.pp_list = pp_list
        data = []
        new_O0, new_M0 = None, None
        for pp in pp_list:
            print(format(pp, 'E'))
            model = DynamicModel(pp=pp, r=self.r, dt=self.dt, beta=self.beta,
                total_time=self.total_time, save_every=self.save_every)
            if new_O0 and new_M0:
                print('new initial conds')
                model.O0 = new_O0
                model.M0 = new_M0
            model.run_stiff()
            new_O0 = model.O[-1]
            new_M0 = model.M[-1]
            data_model = pd.DataFrame({'O': model.O, 'M': model.M, 'C': model.C,
                 'time': model.time, 'pp': pp})
            data.append(data_model)
        self.data = pd.concat(data)
            





