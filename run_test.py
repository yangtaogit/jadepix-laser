#!/usr/bin/env python

import sys,os
import re
import logging
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')
from decode import Decode

def main():


    input = '../data/WeakFe_CHIPA1_180329094037.df'
    output = 'test.root'

    decode = Decode(input,output)
    decode.run()


if __name__ == '__main__':
    main()