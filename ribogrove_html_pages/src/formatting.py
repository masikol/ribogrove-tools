
def format_int_number(number, thousand_separator):

    formatted_number = '{:_}' \
        .format(int(number)) \
        .replace('_', thousand_separator)

    return formatted_number
# end def


def format_float_number(number, thousand_separator, decimal_separator, digits):

    number = round(number, digits)

    format_string = '{:_.%df}' % digits

    formatted_number = format_string \
        .format(number) \
        .replace('_', thousand_separator) \
        .replace('.', decimal_separator)

    return formatted_number
# end def
