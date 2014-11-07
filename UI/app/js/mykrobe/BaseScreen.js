var BaseScreenDirectionNone = 0,
	BaseScreenDirectionForward = 1,
	BaseScreenDirectionBackward = 2;

var BaseScreen = Class.extend({
	init: function(container_) {
		var that = this;
		that.container = $(container_);
		$(window).bind('orientationchange.BaseScreen resize.BaseScreen', function(e) {
            that.update();
        });
       	that._updateDeferred();
       	that.showing = undefined;
	},

	isShowing:function() {
		return this.showing;
	},

	willShow:function(animated_) {
		var that = this;
		that._updateDeferred();
		return this;
	},

	didShow:function(animated_) {
		return this;
	},

	willHide:function(animated_) {
		return this;
	},

	didHide:function(animated_) {
		return this;
	},

	setShowing:function(showing_) {
		var that = this;
		if ( that.showing === showing_ ) {
			return that;
		}
		if ( showing_ ) {
			that.willShow(false);
			$(that.container).show();
			that.showing = showing_;
			that.didShow(false);
		}
		else {
			that.willHide(false);
			$(that.container).hide();
			that.showing = showing_;
			that.didHide(false);
		}
		return that;
	},

	setShowingAnimated:function(showing_,direction_) {
		var that = this,
			l,
			w;
		if ( that.showing === showing_ ) {
			return that;
		}
		w = that.container.innerWidth();
		if ( showing_ ) {
			that.willShow(true);
			that.container.stop(true);
			if ( direction_ === BaseScreenDirectionNone ) {
				$(that.container).fadeIn(600, function() {
					that.didShow(true);
				});
			}
			else {
				l = ( direction_ === BaseScreenDirectionBackward ) ? -w : w;
				that.container.show();
				that.container.css({
					'left':l+'px',
					'opacity':1
				});
				that.container.animate({
					'left':'0px'
				},600, 'easeInOutQuad', function() {
					that.didShow(true);
				});
			}
			// $(that.container).fadeIn(150);
			that.showing = showing_;
		}
		else {
			that.willHide(true);
			that.container.stop(true);
			if ( direction_ === BaseScreenDirectionNone ) {
				$(that.container).fadeOut(600, function() {
					that.didHide(true);
				});
			}
			else {				
				that.container.show();
				that.container.css({
					'left':'0px'
				});
				l = ( direction_ === BaseScreenDirectionBackward ) ? w : -w;
				that.container.animate({
					'left':l+'px'
				},600, 'easeInOutQuad', function() {
					that.container.hide();
					that.didHide(true);
				});
			}
			that.showing = showing_;		
		}
		return that;
	},

	_updateDeferred:function() {
		var that = this;
		setTimeout(function() {
			that.update();
		},0);
		return that;
	},

	update:function() {
		var that = this,
	        windowWidth,
	        windowHeight;

    	windowWidth = $(window).innerWidth();
    	windowHeight = $(window).innerHeight();

	    that.container.css({
            'width': windowWidth + 'px',
            'height': windowHeight + 'px'
        });

		return that;
	}
});