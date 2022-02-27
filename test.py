import modal_straigth_beam as msb
import numpy as np

poutre_model = msb.straight_beam(beam_length = 0.7)
freqs = poutre_model.modal_traction(limite_1 = 'fix')
freqs2 = poutre_model.modal_torsion(limite_1 = 'fix')
