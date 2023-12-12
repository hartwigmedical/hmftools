import logging
import re
import sys

## Disable matplotlib logging
logging.getLogger('matplotlib').disabled = True
logging.getLogger('matplotlib.font_manager').disabled = True
logging.getLogger('matplotlib.pyplot').disabled = True
logging.getLogger('matplotlib.backends.backend_pdf').disabled = True


def reset_logging_basic_config(
    filename: str | None = None,
    level: int = logging.DEBUG,
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
    logging.basicConfig(
        level=level,
        format='%(asctime)s %(threadName)s %(levelname)s | %(name)s.%(funcName)s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )

    logging.captureWarnings(capture=capture_warnings)

## Set default basic config upon module load
reset_logging_basic_config()


class LoggerMixin:
    @property
    def logger(self):
        module_name = re.sub("^.+[.]", "", self.__module__)
        cls_name = self.__class__.__name__
        return logging.getLogger(f"{module_name}.{cls_name}")

    @staticmethod
    def get_class_logger(obj):
        module_name = re.sub("^.+[.]", "", obj.__module__)
        cls_name = obj.__name__
        return logging.getLogger(f"{module_name}.{cls_name}")
