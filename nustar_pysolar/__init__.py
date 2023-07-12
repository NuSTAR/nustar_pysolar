# from pkg_resources import get_distribution, DistributionNotFound
# try:
#     __version__ = get_distribution(__name__).version
# except DistributionNotFound:
#     pass  # package is not installed
    
__version__='0.5'


from . import convert
from . import filter
from . import io
from . import map
from . import planning
from . import tracking
from . import utils
__all__ = ['convert', 'filter','io', 'map', 'planning', 'tracking', 'utils']
