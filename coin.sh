mkdir coin
cd coin
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod +x coinbrew
./coinbrew fetch --no-prompt --main-proj=Clp
./coinbrew build --no-prompt --main-proj=Clp --test
./coinbrew install --no-prompt --main-proj=Clp

