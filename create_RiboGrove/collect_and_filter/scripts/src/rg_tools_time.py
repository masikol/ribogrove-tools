
import sys
import time


def print_time():
    sys.stdout.write(get_time())
    sys.stdout.flush()
# end def

def get_time():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
# end def
