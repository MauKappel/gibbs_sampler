#!/bin/zsh

wget -r -np -R "index.html*" https://services.healthtech.dtu.dk/services/NetMHCII-2.3/suppl/
mv services.healthtech.dtu.dk/services/NetMHCII-2.3/suppl .
rm -r services.healthtech.dtu.dk