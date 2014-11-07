var DropFileScreen = BaseMessageScreen.extend({
	init: function(container_,model_) {
		var that = this,
			chooser;

		that._super( container_ );

		that.functionButtonsContainer = $('.drop-file-screen-screen-message-functions nav',that.container);

		
		chooser = $('<input type="file" accept=".json,.bam,.gz,.fastq" style="display:none" />');
		chooser.on('change',function(e) {
			var path = $(this).val();
			console.log(path);
			$(that).trigger("DropFileScreen:didOpenFileWithPath", path);
			// reset to change event fires next time even if we re-open the same file
			$(this).val(null);
		});
		$(that.container).append(chooser);

		button = $('<li>Browse&hellip;</li>');
		button.on('click',function(e) {
			e.preventDefault();
			chooser.trigger('click'); 
		});
		$('ul',that.functionButtonsContainer).append(button);


		// $('.message',that.container).on('click',function(e) {
		// 	e.preventDefault();
		// 	chooser.trigger('click'); 
		// });	
	

		return that;
	}
});

