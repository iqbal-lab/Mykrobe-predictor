var ResultScreenCell = BaseScreen.extend({
	init: function( title_, container_, footerClasses_) {
		var that = this;
		that._super( container_ );
		that.title = title_;
		that.footerClasses = footerClasses_;
		that.container.mCustomScrollbar({
            verticalScroll:true,
            scrollInertia:0,
            scrollButtons:{
                enable:false
            }
        });
		return that;
	},
	
	update:function() {
		var that = this;
		that.container.mCustomScrollbar("update");		
		return that;
	},

	setViewModel:function(viewModel_) {
		var that = this;
		that.viewModel = viewModel_;
		return that;
	}
});