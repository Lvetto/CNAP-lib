from CLI.utils import *
from CLI.commands import commands
import readline # supposedly fixes arrow behaviour in input()

"""
This program is meant as a Command Line Interface to use some of the functionalities from the CNAP library and access a database of known CNAPs from common structures

"""

if __name__ == "__main__":
    main_loop(commands, compile_names(commands))

