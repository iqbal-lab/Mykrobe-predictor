var ResultScreenDrugsCell = ResultScreenModuleCell.extend({

    setResistanceModel:function(resistanceModel_){
        var that = this;
        that.resistanceModel = resistanceModel_;
        return that;
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

        that._super();
        
        if ( that.resistanceModel.mdr || that.resistanceModel.xdr ) {
            i = (1+that.viewModel.length);
            html = '\
                <div class="module threecol module-evidence module-drug-resistance equalise-height">\
                    <header>\
                        <i class="drug-resistance">R</i><h2>Resistance</h2>\
                    </header>\
                    <article>\
                        <ul>\
            ';
            if ( that.resistanceModel.xdr ) {
                // xdr + mdr, only display xdr
                html += '\
                            <li>Extensively Drug Resistant (XDR)</li>\
                ';
            }
            else {
                // mdr only
                html += '\
                            <li>Multi-Drug Resistant (MDR)</li>\
                ';
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
