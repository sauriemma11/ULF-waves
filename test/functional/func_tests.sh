test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

MAIN_DIR="../../src"
DATA_DIR="../../data/"
FILE="dn_magn-l2-hires_g16_d20230227_v1-0-1.nc"

# ../data/dn_magn-l2-hires_g16_d20230227_v1-0-1.nc

# # successful run -- SAVING PICKLE FILE WHEN TESTING DOESN'T WORK?
# run test_success python $MAIN_DIR/main.py --filename $DATA_DIR$FILE
# assert_exit_code 0


# TO DO: CHECK THAT OUTPUTTING PLOT CORRECTLY


#### testing `filename` input ####
# file does not exist
run test_no_file python $MAIN_DIR/main.py --filename $DATA_DIR'DOES_NOT_EXIST.nc'
assert_exit_code 1

# file is wrong type
run test_file_wrong_type python $MAIN_DIR/main.py --filename $DATA_DIR'wrong_file.txt'
assert_exit_code 2


#### testing `timespan` input ####
# `timespan` isn't a factor of 24
run test_timespan_not_factor python $MAIN_DIR/main.py --filename $DATA_DIR$FILE --timespan 17
assert_exit_code 11

#### testing `fband` input ####
# `fband` has MORE than two elements
run test_fband_wrong_size python $MAIN_DIR/main.py --filename $DATA_DIR$FILE --fband 0.001,0.1,0.5
assert_exit_code 31

# `fband` has LESS than two elements
run test_fband_wrong_size python $MAIN_DIR/main.py --filename $DATA_DIR$FILE --fband 0.001
assert_exit_code 31

# # `fband` wrong type
# run test_fband_wrong_type python $MAIN_DIR/main.py --filename $DATA_DIR$FILE --fband 1,2
# assert_exit_code 32
