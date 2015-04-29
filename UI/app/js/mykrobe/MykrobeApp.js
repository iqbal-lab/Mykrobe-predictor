var MykrobeApp = Class.extend({
	init: function() {
		var that = this,
			gui = require('nw.gui'),
			name = require('./package.json').name,
			toolbar = require('./package.json').window.toolbar,
			mainWindow,
			platform = require("os").platform,
        	version = require('./package.json').displayVersion;

		that.mainWindow = gui.Window.get();
		if ( platform() == 'darwin' ) {
			var nativeMenuBar = new gui.Menu({ 
				type: "menubar" 
			});
			try {
				nativeMenuBar.createMacBuiltin(name);
				that.mainWindow.menu = nativeMenuBar;
			} 
			catch (e) {
			}
		}

		// close about window when app quits
	    that.mainWindow.on('close',function() {
	    	if ( that.aboutWindow ) {
	            that.aboutWindow.close();
	        }
	        // kill processing
	        that.model.cancelLoadFile();
	        this.close(true);
        });

        $('body').addClass(MykrobeTarget.targetName);
        
	    that.model = new Model();

	    $('#version-footer .version-number').html('Version '+version);
	    $('#version-footer .version-footer-about-link').click(function(e){
	        e.preventDefault();
	        that.about();
	    });

	    

	    that.fileDragAndDrop = new FileDragAndDrop($('#file-drag-and-drop'),that.model);
	    that.backgroundScreen = new BackgroundScreen($('#background-screen'));		
	    that.dropFileScreen = new DropFileScreen($('#drop-file-screen'),that.model);
	    that.processingScreen = new ProcessingScreen($('#processing-screen'),that.model);
	    that.resultScreen = new ResultScreen($('#result-screen'));

	    that.currentScreen = undefined;
	    that.screens = [
	    	that.dropFileScreen,
	    	that.processingScreen,
	    	that.resultScreen
	    ];

	    $(that.model).on( "Model:error", function( event_, error_ ) {
	        alert('Error: '+error_.description);
	        // if 
	        if ( !that.resultScreen.isShowing()) {
	        	that.showScreen(that.dropFileScreen);
	        }
	    });

	    $(that.model).on( "Model:willLoad", function( event_ ) {
	        // display progress screen
	        console.log('willLoad');
	        that.showScreen(that.processingScreen);
	    });

	    $(that.model).on( "Model:didLoad", function( event_ ) {
	        // display results screen
	        console.log('didLoad');	      
	        that.resultScreen.setViewModel(that.model);  
	        try {
	        	// that.resultScreen.setViewModel(that.model.json);
	        	that.showScreen(that.resultScreen);
	        }
	        catch (err) {
	        	if ( !that.resultScreen.isShowing()) {
		        	that.showScreen(that.dropFileScreen);
		        }
		        setTimeout(function(){
		        	alert("Error reading file");
		        },0);
	        }
	    });

	    $(that.processingScreen).on( "ProcessingScreen:cancel", function() {
	        console.log('ProcessingScreen:cancel:');
	        that.model.cancelLoadFile();
	        that.showScreen(that.dropFileScreen);
	    });

	    $(that.fileDragAndDrop).on( "FileDragAndDrop:didDropFileWithPath", function( event_, path_ ) {
	        console.log('FileDragAndDrop:didDropFileWithPath:'+path_);
	        that.model.loadFileWithPath(path_);
	    });

	    $(that.dropFileScreen).on( "DropFileScreen:didOpenFileWithPath", function( event_, path_ ) {
	        console.log('DropFileScreen:didOpenFileWithPath:'+path_);
	        that.model.loadFileWithPath(path_);
	    });

	    $(that.resultScreen).on( "ResultScreen:new",function() {
	    	that.showScreen(that.dropFileScreen);
	    });

	    $(that.resultScreen).on( "ResultScreen:save",function(event_, path_) {
	        that.model.saveFileWithPath(path_);
	        // that.dropFileScreen.setShowing(true);
	        // that.resultScreen.setShowing(false);
	    });

	    that.showScreen(that.dropFileScreen);


	    // Listen to `open` event, called when file is opened (dragged onto icon) while app is already open
		gui.App.on('open', function(path_) {
			that.model.loadFileWithPath(path_);
		});

	    // open file on launch
	    console.log('gui.App.argv:'+gui.App.argv);
	    if ( gui.App.argv.length ) {
	    	path_ = gui.App.argv;
	    	if ( that.model.canLoadFileWithPath(path_)) {
	    		that.model.loadFileWithPath(path_);
	    	}
	    }

	    if ( toolbar ) {
	    // testing
	    // that.model.loadFileWithPath('/Users/simon/C00007086_R00000022.bam');
	    // that.model.loadFileWithPath('./model/output.json');
	    // that.error('Test');
	    	if ( kTargetSpeciesTB === MykrobeTarget.species ) {
			    that.model.loadFileWithPath('../spec/fixtures/tb/MTBC_NTM_mixed.json');
		   	}
		   	else {
		   		that.model.loadFileWithPath('../spec/fixtures/staph/mixed_staph.json');
		   	}
	    }
	},

	showScreen:function(screen_) {
		var that = this,
			direction,
			screenIndex,
			currentScreenIndex,
			i,
			screen;
		if ( screen_ === that.currentScreen ) {
			return that;
		}
		if ( that.currentScreen ) {
			screenIndex = that.screens.indexOf(screen_);		
			currentScreenIndex = that.screens.indexOf(that.currentScreen);
			direction = screenIndex > currentScreenIndex ? BaseScreenDirectionForward : BaseScreenDirectionBackward;
			that.currentScreen = screen_;
			for ( i = 0; i < that.screens.length; i++ ) {
				screen = that.screens[i];
				screen.setShowingAnimated(screen === screen_,direction);
			}
			that.backgroundScreen.setShowingAnimated(that.resultScreen !== screen_, BaseScreenDirectionNone);
		}
		else {
			// initial
			that.currentScreen = screen_;
			for ( i = 0; i < that.screens.length; i++ ) {
				screen = that.screens[i];
				screen.setShowing(screen === screen_,direction);
			}
			that.backgroundScreen.setShowing(that.resultScreen !== screen_, BaseScreenDirectionNone);
		}
		return that;
	},

	error:function(error_) {
		var that = this;
		$('#error-modal').modal();
		return that;
	},

	about:function() {
	    var that = this,
	    	name = require('./package.json').name;

	    if ( !that.aboutWindow ) {
	        that.aboutWindow = gui.Window.open('about.html', {
	            'title': 'About '+name,
	            'position': 'center',
	            'width': 600,
	            'height': 400,
	            'toolbar': false,
	            'resizable':false,
	            'focus':true
	        });  

	        that.aboutWindow.on('close',function() {
	            this.close(true);
	            that.aboutWindow = null;
	        });
	    }
	    else {
	        that.aboutWindow.focus();
	    }
	}
});