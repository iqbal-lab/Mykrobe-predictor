var BaseMessageScreen = BaseScreen.extend({
	init: function(container_, model_) {
		var that = this;

		that._super( container_ );
		that.model = model_;
		that.messageContainer = $('.message',that.container);

		return that;
	},

	update:function() {
		var that = this,
	        windowWidth,
	        windowHeight,
	        w,
	        h;
	    windowWidth = $(window).innerWidth();
	    windowHeight = $(window).innerHeight();

	    w = $(that.messageContainer).outerWidth();
	    h = $(that.messageContainer).outerHeight();

	    $(that.messageContainer).css({
	    	'left': (0.5*(windowWidth-w)) + 'px',
	    	'top': (0.5*(windowHeight-h)) + 'px'
	    });

	    that._super();
		return that;
	}
});