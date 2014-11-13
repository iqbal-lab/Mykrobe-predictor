var ResultScreenEvidenceCell = ResultScreenModuleCell.extend({
    minimumNumberOfColumns:function() {
        return 2;
    },

    populateModules:function() {
        var that = this,
            mutation,
            html,
            i,
            j,
            k,
            key,
            iconTitle,
            gene;

        i = 1;

        for ( key in that.viewModel ) {
            evidence = that.viewModel[key];
            title = key;
            iconTitle = title.substr(0,1).toUpperCase();
            html = '\
                <div class="module threecol module-evidence module-generic-'+i+'">\
                    <header>\
                        <i class="generic-'+i+'">'+iconTitle+'</i><h2>'+title+'</h2>\
                    </header>\
            ';
            html += '\
                    <article>\
            ';
            for ( j = 0; j < evidence.length; j++ ) {
                html += '\
                        <ul>\
                ';
                for ( k = 0; k < evidence[j].length; k++ ) {
                    html += '\
                            <li>'+evidence[j][k]+'</li>\
                    ';
                }
                html += '\
                        </ul>\
                ';
            }
            html += '\
                    </article>\
            ';
            that.modules.append($(html));
            i++;
        }

        return that;
    }
});
