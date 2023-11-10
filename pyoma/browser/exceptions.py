from __future__ import division, print_function, unicode_literals


class PyOmaException(Exception):
    pass


class IdException(PyOmaException):
    pass


class InvalidTaxonId(IdException):
    pass


class DBVersionError(PyOmaException):
    pass


class DBConsistencyError(PyOmaException):
    pass


class DBOutdatedError(PyOmaException):
    pass


class InvalidId(IdException):
    pass


class InvalidOmaId(InvalidId):
    pass


class UnknownIdType(IdException):
    pass


class UnknownSpecies(IdException):
    pass


class OutdatedHogId(InvalidId):
    def __init__(self, hog_id):
        super().__init__("Outdated HOG ID: {}".format(hog_id), hog_id)
        self.outdated_hog_id = hog_id


class Singleton(Exception):
    def __init__(self, entry, msg=None):
        super(Singleton, self).__init__(msg)
        self.entry = entry


class NoReprEntry(Exception):
    pass


class AmbiguousID(Exception):
    def __init__(self, message, candidates):
        super(AmbiguousID, self).__init__(message, candidates)
        self.candidates = candidates


class TooUnspecificQuery(PyOmaException):
    def __init__(self, query, hits):
        super().__init__("Too unspecific query '{}' results in {} hits".format(query, hits))
        self.query = query
        self.hits = hits
