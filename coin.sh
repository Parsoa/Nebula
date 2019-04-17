mkdir coin
cp coinbrew ./coin
cd coin
chmod +x coinbrew
./coinbrew fetch --no-prompt --main-proj=Clp
./coinbrew build --no-prompt --main-proj=Clp --test
./coinbrew install --no-prompt --main-proj=Clp

