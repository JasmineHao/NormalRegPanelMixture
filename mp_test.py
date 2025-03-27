
# %%
#import Pool
from multiprocessing import Pool
#Define a worker — a function which will be executed in parallel
def worker(x):
 return x*x
#Assuming you want to use 3 processors
num_processors = 3
#Create a pool of processors
p=Pool(processes = num_processors)
#get them to work in parallel
output = p.map(worker,[i for i in range(0,3)])
print(output)

def worker(x):
    return x*x

from multiprocessing import Pool
import workers
if __name__ ==  '__main__': 
 num_processors = 3
 p=Pool(processes = num_processors)
 output = p.map(workers.worker,[i for i in range(0,3)])
 print(output)
# %%

from numba import jit
import math

# This is the function decorator syntax and is equivalent to `hypot = jit(hypot)`.
# The Numba compiler is just a function you can call whenever you want!
@jit
def hypot(x, y):
    # Implementation from https://en.wikipedia.org/wiki/Hypot
    x = abs(x);
    y = abs(y);
    t = min(x, y);
    x = max(x, y);
    t = t / x;
    return x * math.sqrt(1+t*t)
# %%
