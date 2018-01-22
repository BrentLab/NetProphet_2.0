#!/usr/bin/python
import numpy as np
import json

with open('config.json') as data:
	resources = json.load(data)

for k in resources.keys():
	print(k)