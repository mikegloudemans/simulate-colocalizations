# Author: Mike Gloudemans
# Date created: 7/17/2018
# Load and parse config file.
#

import json

# Function load_config
#
# Input: Filename
# Returns: dictionary containing config settings.
def load_config(filename):

    with open(filename) as data_file:    
            config = json.load(data_file)

    return config
