test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

MAIN_DIR="../../src"  # all paths relative to MAIN_DIR
DATA_DIR="../test/unit/test_data"
OUTPUTS_PATH="../docs"

## downloading test data ##
FILE="dn_magn-l2-hires_g16_d20230227_v1-0-1.nc"  # test dataset
GOOGLE_ID="161_mW7XwKO-Ta1amOsM1VaQjVTs19FXC"

# cd to where main file is
cd $MAIN_DIR

# if does not already exist, download it -- using one example file in google drive
if [[ ! -f "$DATA_DIR/$FILE" ]]; then

    cd $DATA_DIR  # cd to DATA_DIR so data is saved in correct location

    # method for downloading big file from google drive: https://chadrick-kwag.net/wget-google-drive-large-files-bypassing-virus-check/#google_vignette
    wget --save-cookies cookies.txt 'https://docs.google.com/uc?export=download&id='$GOOGLE_ID -O- \
     | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt
    wget --load-cookies cookies.txt -O $FILE \
        'https://docs.google.com/uc?export=download&id='$GOOGLE_ID'&confirm='$(<confirm.txt)
    rm confirm.txt cookies.txt

    # if cd'ed to DATA_DIR, go back to MAIN_DIR
    cd ../../../src
fi

#### successful run ####
run test_success python main.py --filename $DATA_DIR/$FILE
assert_exit_code 0
assert_equal $OUTPUTS_PATH/tau_dict.pkl $( ls $OUTPUTS_PATH/tau_dict.pkl )  # output file created
assert_equal $OUTPUTS_PATH/output_plot.png $( ls $OUTPUTS_PATH/output_plot.png )  # output plot created
rm $OUTPUTS_PATH/tau_dict.pkl $OUTPUTS_PATH/output_plot.png  # remove test files

#### testing `filename` input ####
# file does not exist
run test_no_file python main.py --filename $DATA_DIR/DOES_NOT_EXIST.nc
assert_exit_code 1

# file is wrong type
run test_file_wrong_type python main.py --filename $DATA_DIR/b_av_test_set.txt
assert_exit_code 2

#### testing `timespan` input ####
# `timespan` isn't a factor of 24
run test_timespan_not_factor python main.py --filename $DATA_DIR/$FILE --timespan 17
assert_exit_code 11

#### testing `fband` input ####
# `fband` has MORE than two elements
run test_fband_wrong_size python main.py --filename $DATA_DIR/$FILE --fband 0.001,0.1,0.5
assert_exit_code 31

# `fband` has LESS than two elements
run test_fband_wrong_size python main.py --filename $DATA_DIR/$FILE --fband 0.001
assert_exit_code 31

# # `fband` wrong type
# run test_fband_wrong_type python $MAIN_DIR/main.py --filename $DATA_DIR/$FILE --fband 1,2
# assert_exit_code 32


# REMOVE TEST FILE
rm $DATA_DIR/$FILE