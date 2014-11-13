var gui = require('nw.gui'),
	toolbar = require('./package.json').window.toolbar;

if ( toolbar ) {
	gui.Window.get().showDevTools();
}

$(document).ready(function() {
	var that = this;
	if ( $('body').hasClass('app') ) {
		that.app = new MykrobeApp();
	}    
	else if ( $('body').hasClass('about') ) {
		that.about = new MykrobeAbout();
	}
});

$.modal.defaults = {
	overlay: "#000",        // Overlay color
	opacity: 0.75,          // Overlay opacity
	zIndex: 1,              // Overlay z-index.
	escapeClose: true,      // Allows the user to close the modal by pressing `ESC`
	clickClose: true,       // Allows the user to close the modal by clicking the overlay
	closeText: 'Close',     // Text content for the close <a> tag.
	closeClass: '',         // Add additional class(es) to the close <a> tag.
	showClose: false,        // Shows a (X) icon/link in the top-right corner
	modalClass: "modal",    // CSS class added to the element being displayed in the modal.
	spinnerHtml: null,      // HTML appended to the default spinner during AJAX requests.
	showSpinner: true,      // Enable/disable the default spinner during AJAX requests.
	fadeDuration: null,     // Number of milliseconds the fade transition takes (null means no transition)
	fadeDelay: 1.0          // Point during the overlay's fade-in that the modal begins to fade in (.5 = 50%, 1.5 = 150%, etc.)
};