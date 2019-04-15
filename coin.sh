tar -xzvf coin.tar.gz
#mkdir coin
#cp coinbrew ./coin
cd coin
#chmod +x coinbrew
#./coinbrew fetch --no-prompt --main-proj=Clp
./coinbrew build --no-prompt --skip-update --main-proj=Clp --test
./coinbrew install --no-prompt --skip-update --main-proj=Clp

