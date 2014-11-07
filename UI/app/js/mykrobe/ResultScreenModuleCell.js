var ResultScreenModuleCell = ResultScreenCell.extend({
	init: function(title_, container_, footerClasses_) {
		var that = this;
		that._super( title_, container_, footerClasses_ );
		that.modules = $('.modules',that.container);
		that.masonry = new Masonry( that.modules.get(0), {
			itemSelector: '.module',
			containerStyle: null,
			gutter:20,
			transitionDuration: 0
		});
		return that;
	},
	
    minimumNumberOfColumns:function() {
        return 3;
    },

	update:function() {
		var that = this,
			availableWidth,
			moduleWidth,
			w,
			h;
		
		MIN_MODULE_WIDTH = ((1000-48)-40) / that.minimumNumberOfColumns();
		MODULE_GUTTER = 20;

		availableWidth = $(that.modules).innerWidth() - 48;
		columns = Math.floor((availableWidth+MODULE_GUTTER) / (MIN_MODULE_WIDTH+MODULE_GUTTER));
		columns = Math.min(columns, $('.module',that.modules).length);
		moduleWidth = (availableWidth-(columns-1)*MODULE_GUTTER)/columns;

		// moduleWidth = (availableWidth-4*24)/3;
		$('.module',that.container).css({
			'width':moduleWidth+'px',
			'height':''
		});
		that.masonry.layout();

		// setTimeout(function(){
			maxHeight = 0;
			// console.log('calculating…');
			$('.module.equalise-height',that.container).each(function() {
				oh = $(this).outerHeight();
				// console.log('outerHeight:'+oh);
				maxHeight = Math.max(maxHeight,oh);
			});
			if ( 0 !== maxHeight ) {
				// console.log('maxHeight:'+maxHeight);
				$('.module.equalise-height',that.container).css({
					'height':maxHeight+'px'
				});
				that.masonry.layout();
			}

			that.container.mCustomScrollbar("update");
		// },0)
		
		return that;
	},

	setViewModel:function(viewModel_) {
		var that = this;
			
		that._super(viewModel_);
		$('.module',that.modules).remove();
		
		that.populateModules();

		$("[data-toggle=tooltip]").tooltip();
		that.masonry.reloadItems();
		that._updateDeferred();
		return that;
	},

	populateModules:function() {
		var that = this,
			o,
			iconTitle,
			html;
		for ( i in that.viewModel ) {
			o = that.viewModel[i];
			iconTitle = o['icon-title'] ? o['icon-title'] : '';
			additionalClasses = o['equalise-height'] ? 'equalise-height' : '';

			// console.log(o);
			html = '';

			html += '\
				<div class="module threecol module-'+o['class']+' '+additionalClasses+'">\
                    <header>\
                        <i class="'+o['class']+'">'+iconTitle+'</i><h2>'+o['title']+'</h2>\
                    </header>\
                    <article>\
                        <ul>\
            ';
            for ( j in o['contents']) {
            	c = o['contents'][j];
            	cssClass = c['class'];
            	if ( "string" === typeof(c)) {
	            	html += '\
	            		<li><span class="'+cssClass+'">'+c+'</span></li>\
	            	';
	            }
	            else {
	            	html += '\
		            	<li>\
		            ';
	            	if ( cssClass === o['class'] ) {
	            		// don't show the same icon as the container
	            		html += '\
		            		<span class="'+cssClass+'">'+c['title']+'</span>\
		            	';
	            	}
	            	else {
						html += '\
	            			<span data-toggle="tooltip" title="'+c['tooltip']+'" data-placement="right">\
	            				<span class="'+cssClass+'">'+c['title']+'</span> <i class="'+c['class']+'"></i>\
	            			</span>\
		            	';
	            	}
	            	html += '\
	            		</li>\
	            	';
	            }
            }

            html += '\
                        </ul>\
                    </article>\
                </div>\
			';
			that.modules.append($(html));

		}
		return that;
	}
});