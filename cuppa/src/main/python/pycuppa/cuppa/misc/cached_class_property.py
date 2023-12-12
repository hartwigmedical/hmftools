from functools import cache
def cached_class_property(f):
    return classmethod(property(cache(f)))