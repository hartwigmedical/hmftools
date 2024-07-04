import logging
import sys
from typing import Optional

logging.addLevelName(logging.CRITICAL,  "[CRIT ]")
logging.addLevelName(logging.ERROR,     "[ERROR]")
logging.addLevelName(logging.WARNING,   "[WARN ]")
logging.addLevelName(logging.INFO,      "[INFO ]")
logging.addLevelName(logging.DEBUG,     "[DEBUG]")
logging.addLevelName(logging.NOTSET,    "[NOTSET]")


def initialize_logging(
    filename: Optional[str] = None,
    level: int = logging.DEBUG,
    format: str = "%(asctime)s [%(threadName)s] %(levelname)s %(name)s.%(funcName)s | %(message)s",
    datefmt: str = "%H:%M:%S",
    capture_warnings: bool = True
):
    ## Remove all handlers associated with the root logger object to allow changing the log path at runtime.
    ## Solution from this answer:
    ## https://stackoverflow.com/questions/12158048/changing-loggings-basicconfig-which-is-already-set
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)

    ## Set up handlers
    handlers = []

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(level)
    handlers.append(stream_handler)

    if filename is not None:
        file_handler = logging.FileHandler(filename=filename, mode='w')
        file_handler.setLevel(level)
        handlers.append(file_handler)

    ## Apply config
    logging.basicConfig(level=level, format=format, datefmt=datefmt, handlers=handlers)
    logging.captureWarnings(capture=capture_warnings)


class LoggerMixin:
    @property
    def logger(self):
        cls_name = self.__class__.__name__
        return logging.getLogger(f"{cls_name}")

    @staticmethod
    def get_class_logger(obj):
        cls_name = obj.__name__
        return logging.getLogger(f"{cls_name}")
