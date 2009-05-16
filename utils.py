# this file is part of gosam (generator of simple atomistic models)
# Licence: GNU General Public License version 2

"misc utilities"

import sys

def get_command_line():
    "return command used to call the program (from sys.argv)"
    def quote(s):
        need_quote =  "|&;<>()$`\\' \t\n*?[#~=%"
        t = s.replace('"', '\\"')
        if set(need_quote) & set(s):
            return '"%s"' % t
        else:
            return t
    return " ".join(quote(i) for i in sys.argv)


