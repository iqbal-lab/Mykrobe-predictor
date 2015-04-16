#!/bin/bash
#export CODESIGN_ALLOCATE="/Applications/Xcode.app/Contents/Developer/usr/bin/codesign_allocate"
#Run the following to get a list of certs
#
#security find-identity
app="$1"
#home
identity="Developer ID Application: Simon Heys Limited"

echo "### signing frameworks"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/crash_inspector"
#codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/node-webkit Framework.framework/node-webkit Framework.tmp"
#codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/node-webkit Framework.framework/node-webkit Framework.TOC"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/nwjs Framework.framework/"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/nwjs Helper EH.app/"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/nwjs Helper NP.app/"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Frameworks/nwjs Helper.app/"

echo "### signing binary"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Resources/app.nw/bin/predictor-s-aureus/osx/Mykrobe.predictor.staph"
codesign --force --verify --verbose --sign "$identity" "$app/Contents/Resources/app.nw/bin/predictor-tb/osx/Mykrobe.predictor.tb"

echo "### signing app"
codesign --force --verify --verbose --sign "$identity" "$app"

echo "### verifying signature"
codesign -vvv -d "$app"
sudo spctl -a -vvvv "$app"