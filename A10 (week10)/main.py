from sineModel_MultiRes import sineModel_MultiRes
import sys
sys.path.append('../../software/models/')
import utilFunctions as UF
from scipy.signal import get_window

fs, x = UF.wavread('input_orchestra.wav')

w3 = get_window('blackman',1023)
w2 = get_window('blackman',2047)
w1 = get_window('blackman',4096)

N3 = 1024
N2 = 2048
N1 = 4096

B1 = 300
B2 = 1000
B3 = 22050

t = -100

y = sineModel_MultiRes(x, fs, w1, w2, w3, N1, N2, N3, t, B1, B2, B3)

rs = x - y

UF.wavwrite(y, fs, 'output_orchestra.wav')
