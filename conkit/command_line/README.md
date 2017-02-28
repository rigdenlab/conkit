# Script instructions

If you want to add a new command-line script, please follow these instructions.

1. Start all scripts with the word "conkit"
2. Avoid using non-standard characters to separate words, please only use "_" or "."

If you follow the above steps, then re-running the "python setup.py install" command will automatically place a new script in your bin directory. Note, "_" and "." characters will be converted to a dash, i.e. the "conkit_convert.py" file will be installed as "conkit-convert" script.

