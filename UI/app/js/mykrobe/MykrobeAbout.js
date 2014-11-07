var MykrobeAbout = Class.extend({
	init: function() {
		var that = this,
			name = require('./package.json').name,
        	version = require('./package.json').displayVersion;

        $('body').addClass(MykrobeTarget.targetName);

	    $('.version-name').html(name);
	    $('.version-number').html('Version '+version);
	
	}
});