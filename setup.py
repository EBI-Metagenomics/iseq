#!/usr/bin/env python

import sys
from setuptools import setup

url = {
    "linux": "https://iseq-py.s3.eu-west-2.amazonaws.com/binhouse/hmmsearch_manylinux2010_x86_64",
    "darwin": "https://iseq-py.s3.eu-west-2.amazonaws.com/binhouse/hmmsearch_macosx_10_9_x86_64",
}

if __name__ == "__main__":
    setup()
