import pymc3 as pm
import theano.tensor as tt

with pm.Model() as gp_fit:
    rho = pm.Gamma('rho', 1, 1)
    eta = pm.Gamma('eta', 1, 1)
    K = rho * pm.gp.cov.Matern32(1, rho)
with gp_fit:
    M = pm.gp.mean.Zero()
    det = pm.HalfCauchy('det', 2.5)