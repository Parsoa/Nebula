import math

class NormalDistribution():
    def __init__(self, mean, std):
        self.mean = mean
        self.std = std
        self.var = std ** 2

    def pmf(self, x):
        return (1 / math.sqrt(2 * math.pi * self.var)) * math.exp(-1 * ((x - self.mean) ** 2) / (2 * self.var))

    def log_pmf(self, x):
        return math.log(1 / math.sqrt(2 * math.pi * self.var)) - (((x - self.mean) ** 2) / (2 * self.var))

class ErrorDistribution():
    def __init__(self, p):
        self.p = p

    def pmf(self, x):
        return self.p ** x

    def log_pmf(self, x):
        return x * math.log(self.p)
