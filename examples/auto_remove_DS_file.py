# remove DS file generated automatically by system

import os
p = os.path.abspath('')
os.system("cd " + p)
os.system("find . -name '.DS_Store' -type f -delete")


