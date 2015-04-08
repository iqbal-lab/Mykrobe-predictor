Mykrobe app
===========

The app uses [node-webkit](https://github.com/rogerwang/node-webkit) to provide a cross-platform front-end UI for the Mykrobe binary.

### How do I get set up? ###

First you will need to install [Node.js](http://nodejs.org/), I suggest using homebrew via nvm for this. 

~~~~
brew install nvm
nvm install 0.12.0
nvm use 0.12.0
~~~~

Then run the following terminal commands to install node followed by the dependencies used by the app. From the root folder of the source code:

~~~~
npm install -g nw
npm install -g grunt-cli
npm install -g sass
npm install
bower install
~~~~

Note; `npm` should not require `sudo` to run. If it does (as happened to me) then it will prevent a lot of modules from working. You should [clear out and reinstall](http://stackoverflow.com/a/11178106/857998).

#### Running the app ####

The following terminal command runs the app in the `app` folder. You can terminate it with `ctrl-C`.

~~~~
nw app
~~~~

#### Configuration ####

The app uses [grunt](http://gruntjs.com/) to configure and build the app;

Command 				| Description 																								
----------------------- | --------------------------------------------------------------------------------------------------------- 
`grunt sass`  			| Convert the [sass](http://sass-lang.com/) into the `css` used by the app. 								
`grunt watch` 			| Watch for changes to the sass files and recompile them automatically. Terminate it with `ctrl-C`. 		
`grunt minify` 			| Minify assets and place in the `dist` folder for final packaging. Can be tested using `nw dist`. 	
`grunt mac-icons`		| Create Mac icon files from `icon/mac/icon.pdf`.															
`grunt clean`			| Erase temporary files.																					

### Targets ###

~~~~
grunt set-target
~~~~

Sets the current app target. This reads the targets set in the targets.json file, then updates the `package.json` file accordingly. 

#### Appearance ####

This target value is read at runtime by `MykrobeTarget.js` which sets appropriate global flags. These can be used to tailor the appearance and functionality of individual targets.

#### Binaries ####

TODO: put binaries in target-specific directories, and only copy the specific ones during dist build.
The app distribution build will only contain the binaries for the specified target.

### Deployment ###

~~~~
npm install grunt-contrib-sass --save-dev
grunt dist
~~~~

Compiles, minifies, packages and compiles the binaries. 

#### Mac Deployment ####

A compressed disk image containing the signed app is created in the `dist` folder. The signing identity is set in the `resources/mac/sign.sh` script. The signing step ocassionally fails, just run the task again. The disk image creation procress uses AppleScript so you may a window flash up in the Finder. Once created the disk image is ready to upload.

#### Windows Deployment ####

Binary and associated files are created in the `build/releases/<target>/win/` folder. 

Change the exe icon first using an exe icon resource editor

Then all the files can be bundled into a single executable using the free [Enigma Virtual Box](http://enigmaprotector.com/assets/files/enigmavb.exe)

### Who do I talk to? ###

[si@simonheys.com](mailto:si@simonheys.com)