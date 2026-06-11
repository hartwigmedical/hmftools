import logging
import sys

import torch

LOGGER = logging.getLogger(__name__)

logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s %(levelname)5s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.DEBUG)

LOGGER.setLevel(logging.DEBUG)

logging.getLogger('matplotlib').setLevel(logging.WARNING)


DEVICE = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
