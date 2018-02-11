from neurodesign import geneticalgorithm, generate, msequence
from random import randint
import os
RandSeed=randint(1000, 9999)
#path = '%s/'%(os.path.abspath('..'))
path='/DATA/schonberglab/Classic_CAT/CreateOnsets/probe2/'


EXP = geneticalgorithm.experiment( 
    TR = 1.0, 
    P = [0.5, 0.5], 
    C = [[1.0, 0.0], [0.0, 1.0], [0.5, -0.5]], 
    rho = 0.3, 
    n_stimuli = 2, 
    n_trials = 72, 
    duration = 360.0, 
    resolution = 0.1,
    stim_duration = 1.5, 
    t_pre = 0.0, 
    t_post = 0.0, 
    maxrep = 6, 
    hardprob = True, 
    confoundorder = 3, 
    ITImodel = 'exponential', 
    ITImin = 1.0, 
    ITImean = 3.5, 
    ITImax = 12.0, 
    restnum = 0, 
    restdur = 0.0) 


POP = geneticalgorithm.population( 
    experiment = EXP, 
    G = 20, 
    R = [0, 1, 0], 
    q = 0.01, 
    weights = [0.0, 0.5, 0.25, 0.25], 
    I = 4, 
    preruncycles = 100,
    cycles = 100,
    convergence = 100,
    seed = RandSeed,
    folder = path)


POP.naturalselection()
POP.download()
