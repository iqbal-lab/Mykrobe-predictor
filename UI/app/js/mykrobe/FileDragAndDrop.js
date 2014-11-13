var FileDragAndDrop = Class.extend({
	init: function(container_, model_) {
		var that = this;

		that.model = model_;
		that.container = container_;

		// HTML5 events fire on each dom element
		// so workarounds for the container showing / hiding
		// which cause additional events to fire

		// dragover event, only fire if it's an initial drag over when the container isn't showing

	    $(window).on('dragover', function(e) {
	    	var path, isContainer;
	    	e.preventDefault();
	    	isContainer = ( e.target === that.container.get(0) );
	    	if ( !isContainer ) {
		    	var path;
		    	path = that._firstFilePathForDragEvent(e.originalEvent);
		    	$('body').addClass('hover');
		    	if ( !that.model.canLoadFileWithPath(path)) {
		    		$('body').addClass('hover-invalid');
		    	}
	    	}
	        return false;
	    });

		// dragover event, only fire if it's the container

	    $(window).on('dragleave', function(e) {
	    	var isContainer;
	    	e.preventDefault();
	    	isContainer = ( e.target === that.container.get(0) );
	    	if ( isContainer ) {
	        	$('body').removeClass('hover');
	        	$('body').removeClass('hover-invalid');
	    	}
	        return false;
	    });

	    $(window).on('drop', function(e) {
	    	e.preventDefault();
	    	var path;
	    	path = that._firstFilePathForDragEvent(e.originalEvent);
	        $('body').removeClass('hover');
	        $('body').removeClass('hover-invalid');
	        $(that).trigger("FileDragAndDrop:didDropFileWithPath", path);
	        return false;
	    });
		
		return that;
	},

	_firstFilePathForDragEvent:function(e) {
		return e.dataTransfer.files[0].path;	
	}
});