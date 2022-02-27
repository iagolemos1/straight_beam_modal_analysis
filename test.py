import modal_straigth_beam as msb
import numpy as np

poutre_model = msb.straight_beam(beam_length = 0.7)
freqs1 = poutre_model.modal_traction(n = [1,2,3], limite_0 = 'free', limite_1 = 'fix')
freqs2 = poutre_model.modal_traction(n = [1,2,3], limite_0 = 'fix', limite_1 = 'free')