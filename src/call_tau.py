'''
Module to calculate and combine tau values for specified time interval
-Inputs:
    * dictionary (output from mfa_transform.py) with
        nx3 array of background-subtracted B-field in MFA
        nx1 array of time in dt sec
    * user input for time interval (or a default configuation)
-Outputs:
    * individual returns from functions, particularly
        mx1 tau
        mx1 D_LL
        mx1 PSD
        mx1 time
        mx1 frequencies
-Used in:
    * main.py
-Functions:
    * define_window -- creates a window to split into smaller intervals
    * concat_tau -- gathers all B-field and time arrays into a dictionary
    * call_tau -- calculate tau for each time interval by looping over data set
'''

import get_tau

def define_window(ttimespan):
  window = 0
  return window

def concat_tau(Bfields, t_arrays):
  tau_dict = 0
  return tau_dict

def call_tau():
  concat_tau = []
  for time_interval in ?:
    one_tau = get_tau.get_tau()
    concat_tau.append(one_tau)
  return concat_tau
