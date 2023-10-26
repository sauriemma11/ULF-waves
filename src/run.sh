set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

DIR="../data/"
FILE="dn_magn-l2-hires_g16_d20230227_v1-0-1.nc"
GOOGLE_ID="161_mW7XwKO-Ta1amOsM1VaQjVTs19FXC"


if [ -d "$DIR" ]; then
    cd $DIR  # if data directory exists, go there
else
    mkdir -p $DIR && cd $DIR  # if data directory does NOT exist, make it and go there
fi

if [[ ! -f "$FILE" ]]; then # if does not already exist, download it -- using one example file in google drive

    # method for downloading big file from google drive: https://chadrick-kwag.net/wget-google-drive-large-files-bypassing-virus-check/#google_vignette
    wget --save-cookies cookies.txt 'https://docs.google.com/uc?export=download&id='$GOOGLE_ID -O- \
     | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt
    wget --load-cookies cookies.txt -O $FILE \
        'https://docs.google.com/uc?export=download&id='$GOOGLE_ID'&confirm='$(<confirm.txt)
    rm confirm.txt cookies.txt

    # wget -O $FILE --no-check-certificate "https://drive.google.com/uc?export=download&id="$GOOGLE_ID  # THIS DOESN'T DOWNLOAD .NC FILE CORRECTLY
fi

cd ../src
python main.py --filename $DIR$FILE
