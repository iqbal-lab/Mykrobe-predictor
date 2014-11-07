var ResultScreenSpeciesCell = ResultScreenModuleCell.extend({

    populateModules:function() {
        var that = this,
            html;

        html = '\
            <div class="module threecol module-species equalise-height">\
                <header>\
                    <i class="species">S</i><h2>Species</h2>\
                </header>\
                <article>\
                    <ul>\
                        <li>'+that.viewModel+'</li>\
                    </ul>\
                </article>\
            </div>\
        ';
        that.modules.append($(html));
        
        return that;
    }
});
