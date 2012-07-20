#!/bin/env python

# call this function with sys.argv to parse the config file arguments
def getConfigArguments(filename,format):
    # check format string
    if format not in ['nokeys', 'short', 'full']:
        print "getConfigArguments(): format = {nokeys, short, full}"
        exit(1)

    fp = open(filename)
    lines = fp.readlines()
    fp.close()
    
    # iterate thru lines
    arguments = ""
    for line in lines:
        # strip comment lines
        if line[0] != '#' and line[0] != '\n':
            chunks = line.split() # split into columns
            if len(chunks) < 2:
                print "Line does not contain 2 columns:"
                print line
                exit(1)
            keyword = chunks[0]
            value = chunks[1]
            if format == 'full':
                arguments += "--" + keyword + "=" + value + " "
            if format == 'short':
                arguments += "-" + keyword[0] + " " + value + " "
            if format == 'nokeys':
                arguments += value + " "
    return arguments
