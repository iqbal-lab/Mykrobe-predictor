var ResultScreen = BaseScreen.extend({
	init: function(container_) {
		var that = this,
			chooser,
			button;

		that._super( container_ );

		if ( kTargetSpeciesTB === MykrobeTarget.species ) {
			$('.result-screen-article',that.container).html('\
				<div class="screen screen-view screen-view-all"><div class="modules"></div></div>\
                <div class="screen screen-view screen-view-drugs"><div class="modules"></div></div>\
                <div class="screen screen-view screen-view-species"><div class="modules"></div></div>\
                <div class="screen screen-view screen-view-evidence"><div class="modules"></div></div>\
            ');
			// no inducible for TB
			that.viewAllScreenCell = new ResultScreenModuleCell('View all',$('.screen-view-all',that.container),['susceptible','resistant','inconclusive']);
			that.drugsScreenCell = new ResultScreenDrugsCell('Drugs',$('.screen-view-drugs',that.container),['susceptible','resistant','inconclusive']);
			that.speciesScreenCell = new ResultScreenSpeciesCell('Species',$('.screen-view-species',that.container),['none']);
			that.evidenceScreenCell = new ResultScreenEvidenceCell('Evidence',$('.screen-view-evidence',that.container),['none']);

			that.screens = [
				that.viewAllScreenCell,
				that.drugsScreenCell,
				that.speciesScreenCell,
				that.evidenceScreenCell
			];
		}
		else {
			$('.result-screen-article',that.container).html('\
				<div class="screen screen-view screen-view-all"><div class="modules"></div></div>\
                <div class="screen screen-view screen-view-antibiotic-class"><div class="modules"></div></div>\
                <div class="screen screen-view screen-view-virulence"><div class="modules"></div></div>\
                <div class="screen screen-view screen-view-evidence"><div class="modules"></div></div>\
            ');
			that.viewAllScreenCell = new ResultScreenModuleCell('View all',$('.screen-view-all',that.container),['susceptible','resistant','inducible','inconclusive']);
			that.antibioticClassScreenCell = new ResultScreenModuleCell('Antibiotic class',$('.screen-view-antibiotic-class',that.container),['susceptible','resistant','inducible','inconclusive']);
			that.virulenceScreenCell = new ResultScreenModuleCell('Virulence',$('.screen-view-virulence',that.container),['positive','negative']);
			that.evidenceScreenCell = new ResultScreenEvidenceCell('Evidence',$('.screen-view-evidence',that.container),['none']);

			that.screens = [
				that.viewAllScreenCell,
				that.antibioticClassScreenCell,
				that.virulenceScreenCell,
				that.evidenceScreenCell
			];
		}

		that.screenButtons = [];

		that.screenButtonsContainer = $('nav.result-screen-header-sections',that.container);

		for ( i = 0; i < that.screens.length; i++ ) {
			button = $('<li>'+that.screens[i].title+'</li>');
			$(button).attr('data-index',i);
			$('ul',that.screenButtonsContainer).append(button);
			that.screenButtons.push(button);
			button.on('click',function(e) {
				e.preventDefault();
				that.displayCellAtIndex(parseInt($(this).data('index')));
			});
		}

		chooser = $('<input type="file" nwsaveas="mykrobe.json" style="display:none" />');
		chooser.on('change',function(e) {
			var path = $(this).val();
			console.log(path);
			$(that).trigger("ResultScreen:save", path);
		});
		$('nav.result-screen-header-functions',that.container).append(chooser);

		button = $('<li>Save</li>');
		button.on('click',function(e) {
			e.preventDefault();
			chooser.trigger('click'); 
		});
		$('nav.result-screen-header-functions ul',that.container).append(button);

		button = $('<li>New</li>');
		button.on('click',function(e) {
			e.preventDefault();
			$(that).trigger("ResultScreen:new", false);
		});
		$('nav.result-screen-header-functions ul',that.container).append(button);

		return that;
	},

	willShow:function(animated_) {
		var that = this;
		that._super(animated_);
		that.displayCellAtIndex(0);
		return that;
	},

	displayCellAtIndex:function(index_) {
		var that = this,
			i,
			button,
			footerClasses,
			className,
			windowWidth;
		for ( i = 0; i < that.screenButtons.length; i++ ) {
			button = that.screenButtons[i];
			if ( index_ === i ) {
				button.addClass('selected');
			}
			else {
				button.removeClass('selected');
			}
		}
		that.currentIndex = index_;
		// that.update();
		windowWidth = $(window).innerWidth();
		$('article.result-screen-article',that.container).stop(true);
		$('article.result-screen-article',that.container).animate({
			'left':(-that.currentIndex*windowWidth) + 'px'
		},300, 'easeInOutQuad', function() {
		});
		footerClasses = that.screens[index_].footerClasses;
		$('.result-screen-footer li').each(function() {
			console.log($(this).attr('class'));
			className = $(this).attr('class');
			if ( footerClasses.indexOf(className) !== -1 ) {
				$(this).show();
			}
			else {
				$(this).hide();
			}
			// for ( i = 0; i < footerClasses.length; i++ ) {
			// 	className = footerClasses[i];
			// 	if ( $(this).hasClass(className) ) {
			// 		$(this).show();
			// 	}
			// 	else {
			// 		$(this).hide();
			// 	}
			// }
		});
		return that;
	},

	update:function() {
		var that = this,
	    	combinedWidth,
	        windowWidth,
	        screenCell,
	        headerHeight,
	        footerHeight,
	        i,
	        w,
	        windowHeight;

	    that._super();

	    windowWidth = $(window).innerWidth();
	    windowHeight = $(window).innerHeight();
	    combinedWidth = 0;
	    $('article.result-screen-article',that.container).css({
	    	'left': (-that.currentIndex*windowWidth) + 'px',
	        'width': (i * that.screens.length) + 'px'
	    });
	    w = $(that.screenButtonsContainer).outerWidth();
	    $(that.screenButtonsContainer).css({
	    	'left': (0.5*(windowWidth-w)) + 'px'
	    });
	    headerHeight = $('header.result-screen-header',that.container).outerHeight();
	    footerHeight = $('footer.result-screen-footer',that.container).outerHeight();
	   	for ( i = 0; i < that.screens.length; i++ ) {
	    	screenCell = that.screens[i];
	    	screenCell.container.css({
	    		'top': headerHeight + 'px',
	            'left': (i * windowWidth) + 'px',
	            'width': windowWidth + 'px',
	            'height': (windowHeight - headerHeight - footerHeight) + 'px'
	        });
	        screenCell.update();
	    }

		return that;
	},

	setViewModel:function(viewModel_) {
		var that = this;
		that.model = viewModel_;
		if ( kTargetSpeciesTB === MykrobeTarget.species ) {
			that.viewAllScreenCell.setViewModel(that.model.viewAllModel);
			that.drugsScreenCell.setResistanceModel(that.model.drugsResistanceModel);
			that.drugsScreenCell.setViewModel(that.model.drugsModel);
			that.speciesScreenCell.setViewModel(that.model.speciesModel);
			// that.antibioticClassScreenCell.setViewModel(that.model.antibioticClassModel);
			// that.infectionSiteScreenCell.setViewModel(that.model.infectionSiteModel);
			// that.virulenceScreenCell.setViewModel(that.model.virulenceMarkersModel);
			that.evidenceScreenCell.setViewModel(that.model.evidenceModel);
		}
		else {
			that.viewAllScreenCell.setViewModel(that.model.viewAllModel);
			that.antibioticClassScreenCell.setViewModel(that.model.antibioticClassModel);
			that.virulenceScreenCell.setViewModel(that.model.virulenceMarkersModel);
			that.evidenceScreenCell.setViewModel(that.model.evidenceModel);
		}
		return that;
	}
});