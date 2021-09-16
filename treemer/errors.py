from Bio.Phylo.BaseTree import Clade


class TipNotMatchedError(Exception):
    def __init__(self, tip: Clade):
        self.tip = tip

    def __str__(self):
        return "{} is not found in alignment".format(self.tip.name)


class SiteNotValidError(ValueError):
    def __init__(self, site: int):
        self.site = site

    def __str__(self):
        return "{} is either not a integer or not within the range of alignment".format(self.site)


class MethodNotFoundError(Exception):
    def __init__(self, method: str, available: list):
        self.method = method
        self.available = available

    def __str__(self):
        return "{} is not one of the available methods {}".format(self.method, self.available)


class ThresholdError(ValueError):
    def __init__(self, threshold, message: str, bound):
        self.threshold = threshold
        self.message = message
        self.bound = bound

    def __str__(self):
        return "Threshold ({}) cannot be {} bound ({})".format(
            self.threshold, self.message, self.bound
        )


class FastMLError(Exception):
    def __init__(self, message: str) -> None:
        self.message = message

    def __str__(self):
        return self.message
