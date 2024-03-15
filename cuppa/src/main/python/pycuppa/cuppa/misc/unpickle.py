from __future__ import annotations

from typing import Any

from cuppa.logger import LoggerMixin


class MissingAttrHandler(LoggerMixin):
    def __init__(self, obj: object):
        self.obj = obj

    def check_and_set_missing_attr(self, name: str, value: Any) -> None:
        if hasattr(self.obj, name):
            return

        self.logger.warning(
            f"In `{self.obj.__class__.__name__}` object, the missing attr `{name}` was assigned the value `{str(value)}`. "
            f"`{name}` is likely a new attr that did not exist when this object was pickled"
        )

        setattr(self.obj, name, value)
