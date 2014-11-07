var kProcessingScreenProgressStrokeWidth = 36;

var ProcessingScreen = BaseMessageScreen.extend({
	init: function(container_, model_) {
		var that = this,
			diameter,
			button;

		that._super( container_, model_ );

		that.indicatorDots = $('.message-dots',that.container);

		that.progressContainer = $('.progress-indicator',that.container);
		that.paper = Raphael(that.progressContainer.get(0), $(that.progressContainer).innerWidth(), $(that.progressContainer).innerHeight());
		diameter = $(that.progressContainer).innerWidth();

		that.paper.customAttributes.arc = function(xloc, yloc, value, total, R) {
            var alpha = 360 / total * value,
                a = (90 - alpha) * Math.PI / 180,
                x = xloc + R * Math.cos(a),
                y = yloc - R * Math.sin(a),
                path;
            if (total == value) {
                path = [
                    ["M", xloc, yloc - R],
                    ["A", R, R, 0, 1, 1, xloc - 0.01, yloc - R]
                ];
            } else {
                path = [
                    ["M", xloc, yloc - R],
                    ["A", R, R, 0, +(alpha > 180), 1, x, y]
                ];
            }
            return {
                path: path
            };
        };

        // $('.message',that.container).css({
        // 	'background':''
        // });

        that.paper.path().attr({
            "fill": '#fff',
            'stroke':'#fff',
            "stroke-width": kProcessingScreenProgressStrokeWidth,
            arc: [diameter/2, diameter/2, 100, 100, diameter/2-kProcessingScreenProgressStrokeWidth/2]
        });

        that.indicatorArc = that.paper.path().attr({
            "stroke": '#ceccc6',
            "stroke-width": kProcessingScreenProgressStrokeWidth,
            arc: [diameter/2, diameter/2, 0, 100, diameter/2-kProcessingScreenProgressStrokeWidth/2]
        });

		that.indicatorText = that.paper.text(diameter/2, diameter/2-6, "X").attr({
			"font-size": 96,
			"font-family": 'Bryant2-Bold',
			"fill": '#ceccc6'
		});

		that.setProgress(0);
		// that.setProgressAnimated(33);

		that.functionButtonsContainer = $('.processing-screen-message-functions nav',that.container);

		button = $('<li>Cancel</li>');
		button.on('click',function(e) {
			e.preventDefault();
			$(that).trigger("ProcessingScreen:cancel", false);
		});
		$('ul',that.functionButtonsContainer).append(button);

		$(that.model).on( "Model:willLoad", function( event_ ) {
			that.processingStartDate = new Date();
			that.processingPreviousProgressDate = that.processingStartDate;
			that.processingAverageTimePerPercentage = 0;
			that.processingLastTimePerPercentage = 0;
			that.SMOOTHING_FACTOR = 0.02;
			that.stopProgressAnimation();
			that.setProgress(0);
			that.processingLastPercent = 0;
		});

		$(that.model).on( "Model:progress", function( event_, progress_ ) {
	        var percent = 100 * progress_.progress / progress_.total;
	        var now = new Date();
	        var timeSinceLastProgress = now.getTime() - that.processingPreviousProgressDate.getTime();
	        that.processingLastTimePerPercentage = timeSinceLastProgress / (percent - that.processingLastPercent);
	        if ( 0 === that.processingAverageTimePerPercentage ) {
	        	that.processingAverageTimePerPercentage = that.processingLastTimePerPercentage;
	        }
	        that.processingAverageTimePerPercentage = that.SMOOTHING_FACTOR * that.processingLastTimePerPercentage + (1-that.SMOOTHING_FACTOR) * that.processingAverageTimePerPercentage;
	        that.setProgressAnimated(percent,(now.getTime() - that.processingPreviousProgressDate) / 1000);
	        that.processingPreviousProgressDate = now;
	        that.processingLastPercent = percent;

	        // console.log('timeSinceLastProgress:'+timeSinceLastProgress+' percent:'+percent+' that.processingLastPercent:'+that.processingLastPercent);
	    });

		return that;
	},

	stopProgressAnimation:function() {
		var that = this;
		if ( that.progressTween ) {
			that.progressTween.kill();
			that.progressTween = null;
		}
		return that;
	},

	setMessageTitle:function(title_) {
		var that = this;
		$('.message-title',that.container).html(title_);
		return that;
	},

	// http://stackoverflow.com/questions/2779600/how-to-estimate-download-time-remaining-accurately
	// averageSpeed = SMOOTHING_FACTOR * lastSpeed + (1-SMOOTHING_FACTOR) * averageSpeed;
	// 0.005 provides a pretty good smoothing value for an average download speed.

	setProgress:function(progress_) {
		var that = this,
			diameter = $(that.progressContainer).innerWidth(),
			stroke = kProcessingScreenProgressStrokeWidth,
			p;

		that.progress = progress_;

        that.indicatorArc.attr({
            arc: [diameter/2, diameter/2, that.progress, 100, diameter/2-stroke/2]
        });

        // milliseconds
        if ( that.progress > 0 ) {
        	if ( 100 === that.progress ) {
        		that.setMessageTitle('Check species and<br>scan for resistance');
        		that.indicatorText.attr({text:''});
        		that.indicatorDots.show();
        	}
        	else {
        		that.indicatorDots.hide();
        		that.setMessageTitle('Constructing genome');
        		/*
		        var now = new Date();
		        var timeElapsed = now.getTime() - that.processingStartDate.getTime();
		        var timeRemaining = (100 - that.progress) * timeElapsed / that.progress;
				var timeRemainingAlt = (100 - that.progress) * that.processingAverageTimePerPercentage;

			    var seconds = Math.floor(timeRemaining / 1000) % 60 ;
				var minutes = Math.floor((timeRemaining / (1000*60)) % 60);
		        //that.setMessageTitle('Analysing '+Math.round(percent)+'% '+humanize.relativeTime(timeToFinish));
			    var secondsAlt = Math.floor(timeRemainingAlt / 1000) % 60 ;
				var minutesAlt = Math.floor((timeRemainingAlt / (1000*60)) % 60);

				seconds = secondsAlt;
				minutes = minutesAlt;

		        // that.setMessageTitle(minutes+':'+('0'+seconds).substr(-2)+' remaining / '+minutesAlt+':'+('0'+secondsAlt).substr(-2)+' remaining');

		        if ( minutes > 0 ) {
		        	var description = minutes + ' minutes';
		        	if ( 1 === minutes ) {
		        		description = 'one minute';
		        	}
		        	that.setMessageTitle('About '+description+' to go');
		        }
		        else {
		        	if ( seconds > 30 ) {
		        		that.setMessageTitle('Less than a minute to go');
		        	}
		        	else if ( seconds > 1 ) {
			        	that.setMessageTitle(seconds+' seconds to go');
		        	}
		        	else {
		        		that.setMessageTitle('1 second to go');
		        	}
		        }
		        */
		        that.indicatorText.attr({text: Math.round(that.progress)+'%'});
		    }
	    }
	    else {
	    	that.setMessageTitle('Analysing');
	    	that.indicatorText.attr({text:''});
	    	that.indicatorDots.show();
	    }
		return that;
	},

	setProgressAnimated:function(progress_,duration_) {
		var that = this;
		duration_ = duration_ ? duration_ : 0.3;
		that.stopProgressAnimation();
        that.progressTween = TweenLite.to(that, duration_, {
            setProgress:progress_,
            ease:Linear.easeNone
            // ease:Quad.easeOut
        });
        return that;
	},

	getProgress:function() {
		var that = this;
		return that.progress;
	}
});