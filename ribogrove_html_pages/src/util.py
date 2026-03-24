
import pandas as pd


def is_validlike_species_name(input_arg):
    if isinstance(input_arg, str):
        return _seems_like_valid(input_arg)
    elif isinstance(input_arg, pd.Series):
       return pd.Series([
            _seems_like_valid(x) for x in input_arg
       ])
    else:
        raise TypeError('Incompatible input_arg type: "{}"'.format(type(input_arg)))
   # end if
# end def


def _seems_like_valid(species_str):
    return not species_str.endswith('sp.') \
        and not ' sp. ' in species_str
# end def
