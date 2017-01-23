
import numpy as np
import numpy.linalg as la

from cgs import *
from pyminiconf import *

params = pyminiconf_read(open("stary.txt","r"))

print params['temp_eff']
