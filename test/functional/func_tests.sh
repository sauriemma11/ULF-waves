test -e ssshtest || wget -qhttps://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

#run <test name> <program> <argument 1> <argument 2> <...>
# assert_exit_code 0

run test_main python main.py
assert_exit_code 0