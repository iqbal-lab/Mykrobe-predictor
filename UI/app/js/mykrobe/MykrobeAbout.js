var MykrobeAbout = Class.extend({
	init: function() {
		var that = this,
			gui = require('nw.gui'),
			name = require('./package.json').name,
        	version = require('./package.json').displayVersion;

        $('body').addClass(MykrobeTarget.targetName);

	    $('.version-name').html(name);
	    $('.version-number').html('Version '+version);
		
		$('a[data-open-external]').click(function(e){
			e.preventDefault();
			gui.Shell.openExternal($(this).attr('href'));
		});
	}
});