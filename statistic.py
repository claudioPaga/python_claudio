import numpy as np
import scipy.misc


class Stat:
    def __init__(self, name):
        self.name = name

    def resid(self, data, model, err=None):
        """
        For the least-squares minimisers that require residuals.
        """
        pass
        return None

    def __call__(self, data, model, err=None):
        """
        Evaluate the statistic
        """
        pass
        return None


class ChiSqr(Stat):
    """
    Classical Chi-square statistic, with weights provided by the
    data variance. That is,

    chisqr = Sum_i ( (data(i) - model(i) / err(i))**2 )
    """
    def __init__(self, name="ChiSqr"):
        Stat.__init__(self, name)

    # resid is provided for levmar and cmpfit
    def resid(self, data, model, err):
        return (data - model) / err

    def __call__(self, data, model, err):
        chisqr = self.resid(data, model, err)**2.0
        return chisqr.sum()


class ChiSqrP(Stat):
    """
    Pearson's Chi-square statistic, with weights provided by the model.
    That is,

    chisqr_p = Sum_i ( (data(i) - model(i) )**2 / model(i) )
    """
    def __init__(self, name="ChiSqrP"):
        Stat.__init__(self, name)

    # resid is provided for levmar and cmpfit
    def resid(self, data, model, err=None):
        msqrt = np.sqrt(model)
        chi = (data - model) / msqrt
        return chi

    def __call__(self, data, model, err=None):
        chisqr = ((data - model)**2.0 / model)
        return chisqr.sum()


class Cash(Stat):
    """
    Cash ML statistic for Poisson data.

    C = 2 Sum_i ( model(i) - data(i) log (model(i)) + log(data(i)!) )

    NB - the log(data(i)!) term is usually dropped as it is constant
    for the given dataset.

    """

    def __init__(self, name="Cash", data=None, logzero=-32.0):
        Stat.__init__(self, name)

        self.logzero = logzero

        # log (data(i)!)
        self.logfactdata = 0.0
        if data is not None:
            self.logfactdata = np.log(scipy.misc.factorial(data))

    # resid is provided for levmar and cmpfit
    def resid(self, data, model, err=None):

        # From Craig :
        # Assume you are going to do something like TOTCASH = SUM(CASH_i).
        # If you define
        # RESID_i = SIGN(DATA_i - MODEL_i) * SQRT(CASH_i)
        # then
        # TOTCASH = SUM(RESID_i^2)
        # which looks like a sum if squares.

        lmodel = np.zeros(model.shape)

        # trap -ve values before logs are taken
        # lmodel[model <= 0.0] = -32.0
        lmodel[model <= 0.0] = self.logzero

        lmodel[model > 0.0] = np.log(model[model > 0.0])

        cash_i = 2.0 * (model - data * lmodel + self.logfactdata)

        rtn = np.sign(data - model) * np.sqrt(cash_i)

        return rtn

    def __call__(self, data, model, err=None):

        lmodel = np.zeros(model.shape)

        # trap -ve values before logs are taken
        # lmodel[model <= 0.0] = -32.0
        lmodel[model <= 0.0] = self.logzero

        lmodel[model > 0.0] = np.log(model[model > 0.0])

        cash_i = model - data * lmodel + self.logfactdata

        return 2.0 * cash_i.sum()


class CStat(Stat):
    """
    cstat - xspec's version of the Cash ML statistic for Poisson data.

    C = 2 Sum_i ( model(i) - data(i) + data(i) (log(data(i)) - log(model(i))) )

    The advantage of this statistic over Cash is that it tends to ChiSqr when
    the number of counts per bin is large, so a goodness of fit can be inferred
    in the usual way by comparing the returned statistic with the number of
    degrees of freedom in the fit.
    """

    def __init__(self, name="CStat", logzero=-32.0):
        Stat.__init__(self, name)

        self.logzero = logzero

    # resid is provided for levmar and cmpfit
    def resid(self, data, model, err=None):

        lmodel = np.zeros(model.shape)
        ldata = np.zeros(data.shape)

        # trap -ve values before logs are taken
        # lmodel[model <= 0.0] = -32.0
        lmodel[model <= 0.0] = self.logzero
        lmodel[model > 0.0] = np.log(model[model > 0.0])

        # ldata[data <= 0.0] = -32.0
        ldata[data <= 0.0] = self.logzero
        ldata[data > 0.0] = np.log(data[data > 0.0])

        # according to the xspec manual, this definition approaches
        # chisq in the limit of a large number of counts
        cstat_i = model + data * ((ldata - lmodel) - 1.0)

        ohoh = np.argwhere(cstat_i < 0.0)

        if len(ohoh) > 0:
            # floating point round off seems to cause a slightly
            # negative cstat_i value. E.g. when data=1202.0 and
            # model=1202.00001.
            # cstat_i should (probably) never be negative provided
            # data >=0 and model >= 0
            print("WARNING: Cstat_i < 0!")
            print("i:", ohoh)
            np.set_printoptions(precision=32, floatmode='fixed')
            print("Cstat_i:", cstat_i[ohoh])
            print("model[i]:", model[ohoh])
            print("lmodel[i]:", lmodel[ohoh])
            print("data[i]:", data[ohoh])
            print("ldata[i]:", ldata[ohoh])
            np.set_printoptions(precision=8, floatmode='maxprec_equal')

            cstat_i[ohoh] = 0.0

        cstat_i *= 2.0

        rtn = np.sign(data - model) * np.sqrt(cstat_i)

        return rtn

    def __call__(self, data, model, err=None):

        lmodel = np.zeros(model.shape)
        ldata = np.zeros(data.shape)

        # trap -ve values before logs are taken
        # xspec uses -32.0 here - why not 0.0 ?
        # lmodel[model <= 0.0] = -32.0
        lmodel[model <= 0.0] = self.logzero
        lmodel[model > 0.0] = np.log(model[model > 0.0])

        # ldata[data <= 0.0] = -32.0
        ldata[data <= 0.0] = self.logzero
        ldata[data > 0.0] = np.log(data[data > 0.0])

        # according to the xspec manual, this definition approaches
        # chisq in the limit of a large number of counts
        cstat_i = model + data * ((ldata - lmodel) - 1.0)

        return 2.0 * cstat_i.sum()


class Whittle(Stat):
    """
    Whittle statistic for fitting a single powerspectrum :

    W = 2 Sum_j ( I_j / S_j + log(S_j) )

    where I_j = powerspectrum (i.e. data being fit)
          S_j = model of the powerspectrum being fit

    See S. Vaughan 2010 MNRAS 402 307
    """

    def __init__(self, name="Whittle"):
        Stat.__init__(self, name)

    def __call__(self, data, model, err=None):

        # trap zeros before reciprocal is taken
        invmodel = np.zeros(model.shape)

        invmodel[model > 0.0] = 1.0 / model[model > 0.0]

        lmodel = np.zeros(model.shape)

        # trap -ve values before logs are taken
        lmodel[model <= 0.0] = -32.0

        lmodel[model > 0.0] = np.log(model[model > 0.0])

        fitstat = data * invmodel + lmodel

        return 2.0 * fitstat.sum()


class WhittleM(Stat):
    """
    Whittle-like statistic for fitting the average of M powerspectra.

    Derived from a corrected version of eqn 8 of Appourchaux
    2003 AA RN 412 903, which defines the probability density function
    for the average of M periodograms to be:

    p(I_j) = M**(M - m) (1/Gamma(M)) (I_j**(M-1) / S_j**M) exp(-M I_j / S_j)

    The new corrected derivation suggests that the first term should be
    M**M and not M**(M-1), as in Appourchaux 2003. This new derivation
    has been confirmed by simulations. The value of when m can be set to
    0 or 1 to switch between derivations (m=0 for new derivation, m=1 for
    Appourchaux derivation). This only effects the constant c(M) below,
    which acts as a global offset to the log-likelihood.


    This gives a log-likelihood of :

    W = 2M Sum_j ( I_j / S_j + log(S_j) + ((1 / M) - 1) log(I_j) + c(M) )

    where I_j = average of M powerspectra (i.e. data being fit)
          S_j = model of the average powerspectrum being fit
          c(M) = a constant for a fixed M
               = -(1/M) log(M**(M-m) / Gamma(M))
               = (1/M) [log(Gamma(M) - (M-m)log(M)]

    When M is large (>50), we use an approximation for log(Gamma(M))
    (based on Stirling approx for the factorial) to avoid numeric
    overflows in Gamma(M)
    """

    def __init__(self, name="WhittleM", M=1, m=0, c=None):
        """
        M = number of PSDs in average.

        m = 0 to compute c(M) using new derivation of p(I_j)
            1 ro compute c(M) using  Appourchaux derivation of p(I_j)

        c = override for c calculation. Used if WhittleM is constructed
        with c set to a value other than None.
        """
        Stat.__init__(self, name)

        def log_Gamma(M):
            import math
            if (M < 51):
                lg = np.log(math.gamma(M))
            else:
                # log(Gamma(M) using Stirling's approx for the factorial
                # in Gamma
                lg = M * np.log(M) - M - 0.5 * np.log(M / (2. * math.pi)) + \
                     1.0/(12.0 * M)

            return lg

        M = int(M)
        m = int(m)

        self.M = M
        self.m = m

        if c is not None:
            self.c = float(c)

        else:
            if M == 1:
                self.c = 0.0

            else:
                log_Gamma_M = log_Gamma(M)

                self.c = (log_Gamma_M - (M - m) * np.log(M)) / M

    def __call__(self, data, model, err=None):

        # trap zeros before reciprocal is taken
        invmodel = np.zeros(model.shape)

        invmodel[model > 0.0] = 1.0 / model[model > 0.0]

        lmodel = np.zeros(model.shape)

        # trap -ve values before logs are taken
        lmodel[model <= 0.0] = -32.0

        lmodel[model > 0.0] = np.log(model[model > 0.0])

        ldata = np.zeros(data.shape)

        # trap -ve values before logs are taken
        ldata[data <= 0.0] = -32.0

        ldata[data > 0.0] = np.log(data[data > 0.0])

        fitstat = data * invmodel + lmodel + ((1.0 / self.M) - 1) * ldata + self.c

        return 2.0 * self.M * fitstat.sum()


class WhittleM_BV(Stat):
    """
    Whittle statistic for fitting the average of M powerspectra :

    W = nu Sum_j ( I_j / S_j + log(S_j) + ((2 / nu) - 1) log(I_j) + c(nu) )

    where I_j = average of M powerspectra (i.e. data being fit)
          S_j = model of the average powerspectrum being fit
          nu  = 2 M
          c(nu) = a constant for a fixed nu (and is dropped from the
                  fit statistic below)

    See appendix A of D. Barret and S. Vaughan 2012 ApJ 746 131
    """

    def __init__(self, name="WhittleM", M=1):
        Stat.__init__(self, name)
        # nu = 2M dof
        self.nu = 2 * M

    def __call__(self, data, model, err=None):

        # trap zeros before reciprocal is taken
        invmodel = np.zeros(model.shape)

        invmodel[model > 0.0] = 1.0 / model[model > 0.0]

        lmodel = np.zeros(model.shape)

        # trap -ve values before logs are taken
        lmodel[model <= 0.0] = -32.0

        lmodel[model > 0.0] = np.log(model[model > 0.0])

        ldata = np.zeros(data.shape)

        # trap -ve values before logs are taken
        ldata[data <= 0.0] = -32.0

        ldata[data > 0.0] = np.log(data[data > 0.0])

        fitstat = data * invmodel + lmodel + (2.0 / self.nu - 1) * ldata

        return self.nu * fitstat.sum()
