from grasp import *

charges = np.arange(6,13)
masses_number = [12,14,16,19,20,23]
masses_amu  = [12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, 22.98976928, 24.3050]
nuclei = [Rnucleus(Z = charge, A = number, neutralMass = amu, I = 1, NDM = 1, NQM = 1) for charge,number, amu in zip(charges,masses_number,masses_amu)]

csf = Rcsfgenerate(core = 'None', ordering = 'Default', csflist = ['1s(2,i)2s(1,i)'], activeset = [2,2],jlower = 1, jhigher = 3, exc = 0) + Rcsfgenerate(core = 'None', ordering = 'Default', csflist = ['1s(2,i)2p(1,i)'], activeset = [2,2], jlower = 1, jhigher = 3, exc = 0, write_csf = 'DF.c')


