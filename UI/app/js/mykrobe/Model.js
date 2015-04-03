var Model = Class.extend({
    init: function() {
    	var that = this;
    	return that;
    },

    canLoadFileWithPath:function(path_) {
        var that = this,
            extension;
        extension = that._extensionForPath(path_);
        return ( '.json' === extension || '.bam' === extension || '.gz' === extension || '.fastq' === extension );
    },

    loadFileWithPath:function(path_) {
        var that = this,
            extension;
        that.cancelLoadFile();
        extension = that._extensionForPath(path_);
        if ( !that.canLoadFileWithPath(path_)) {
            $(that).trigger("Model:error", {
                'path':path_, 
                'description':"Cannot load file with extension '"+extension+"'"
            });
            return that;
        }
        if ( '.json' === extension ) {
            that._loadJsonWithPath(path_);
        }
        else {
            that._loadFileWithPath(path_);
        }
        return that;
    },

    _extensionForPath:function(path_) {
        var that = this,
            path,
            extension;

        extension = path_.substr(path_.lastIndexOf('.'));
        // path = require('path');
        // extension = path.extname(path_).toLowerCase();
        return extension.toLowerCase();
    },

    _loadJsonWithPath:function(path_) {
        var that = this,
            isValid;
        console.log('_loadJsonWithPath:'+path_);
        $.ajax({
            type: "GET",
            url: path_,
            success:function(dataString) {
                that._loadString(dataString,path_);
            },
            error:function(jqXHR, textStatus, errorThrown) {
                $(that).trigger("Model:error", {
                    'path':path_, 
                    'description':"Parsing error: "+textStatus.reason
                });
            }
        });
        return that;
    },

    _pathToBin:function() {
        var that = this,
            path = require("path"),
            platform = require("os").platform(),
            platformFolder = '',
            chmodSync = require("fs").chmodSync,
            pathToBin = '';

        // will be 'darwin', 'win64' or 'win64'
        // FIXME: not tested for Linux
        platformFolder = platform;

        switch (platform) {
            case "darwin":
                // use 'osx' folder for Mac
                platformFolder = 'osx';
                break;
        };

        if ( kTargetSpeciesTB === MykrobeTarget.species ) {
            pathToBin = path.join('bin',MykrobeTarget.targetName,platformFolder,'Mykrobe.predictor.tb');
        }
        else {
            pathToBin = path.join('bin',MykrobeTarget.targetName,platformFolder,'Mykrobe.predictor.staph');
        }

        chmodSync(pathToBin, 0755);
        return pathToBin;
    },

    cancelLoadFile:function() {
        var that = this;
        if ( that.child ) {
            //that.child.stdout.on("data", that._stdoutData);
            // that.child.stdout.removeListener('data',that._stdoutData);
            // that.child.stdout.unpipe();
            // that.child.stdio.unpipe();
            that.child.kill();
            that.child = null;
        }
        return that;
    },

    _loadString:function(string_,path_) {
        var that = this,
            obj, 
            first, 
            last, 
            species,
            isValid,
            extracted;

        try {
            // extract just the portion in curly braces {}
            first = string_.indexOf('{');
            last = string_.lastIndexOf('}');
            extracted = string_.substr(first,1+last-first);

            // replace escaped tabs, quotes, newlines
            extracted = extracted.replace(/\\n/g,'\n').replace(/\\t/g,'\t').replace(/\\"/g,'"');
            console.log(extracted);
            that.json = JSON.parse(extracted);

            // update all our model objects;
            that._createModels();

            // valid if there is one species, 'S.aureus'
            // TODO: this should be much more flexible
            isValid = false;
            requiredSpecies = [];
            if ( kTargetSpeciesTB === MykrobeTarget.species ) {
                requiredSpecies = 'MTBC';
            }
            else {
                requiredSpecies = 'Staphylococcus aureus';
            }
            // isValid = requiredSpecies.indexOf(that.species[0]) !== -1;
            isValid = _.includes(that.phyloGroup,requiredSpecies);
            if ( isValid ) {
                $(that).trigger("Model:didLoad", false);
            }
            else {
                /*
                "This sample seems to be {WHATEVER NAME IS IN JSON}, not S. aureus, and therefore the predictor does not give susceptibility predictions"

                In an ideal world, if it was Penicillin or Methicillin resistant, I would say
                "..is methicillin/penicillin resistant S. lugdunensis"
                */

                sampleName = that.phyloGroup.join(' / ');
                if ( that.resistant.length ) {
                    resistantName = that.resistant.join(' / ');
                    sampleName = resistantName +' resistant ' + sampleName;
                }

                $(that).trigger("Model:error", {
                    'path':path_, 
                    'description':'This sample seems to be '+sampleName+', not '+requiredSpecies+', and therefore the predictor does not give susceptibility predictions' 
                });
            }

        } 
        catch (e) {
            console.error("Parsing error:", e); 
            console.log(string_);
            $(that).trigger("Model:error", {
                'path':path_, 
                'string':string_,
                'description':"Parsing error: "+e
            });
        }
        return that;
    },

    _loadFileWithPath:function(path_) {
        var that = this,
            pathToBin = that._pathToBin(),
            dirToBin,
            path = require('path'),
            didReceiveError,
            args,
            digitGroups,
            progress,
            total,
            spawn = require("child_process").spawn;

        that.cancelLoadFile();
        console.log('_loadFileWithPath:'+path_);
        $(that).trigger("Model:willLoad", false);
        dirToBin = path.join(process.cwd(),'bin',MykrobeTarget.targetName);
        that.didReceiveError = false;
        args = ['--file', path_, '--install_dir', dirToBin, '--format', 'JSON', '--progress'];
        that.child = spawn(pathToBin, args);
        that.jsonBuffer = '';
        that.isBufferingJson = false;
        that.processExited = false;
        that.child.stdout.on("data", function(data) {
            if ( that.didReceiveError ) {
                return;
            }
            dataString = data.toString('utf8');
            console.log(dataString);
            if ( dataString.indexOf('Progress') === 0 ) {
                // console.log('progress');
                // we get a string like "Progress 1000000/1660554"
                // extract groups of digits
                digitGroups = dataString.match(/\d+/g);
                if ( digitGroups.length > 1 ) {
                    progress = parseInt(digitGroups[0]);
                    total = parseInt(digitGroups[1]);
                    // console.log('progress:'+progress);
                    // console.log('total:'+total);
                    $(that).trigger("Model:progress", {
                        'progress':progress,
                        'total':total
                    });
                }
            }
            
            if ( that.isBufferingJson ) {
                that.jsonBuffer += dataString;
            }
            else if ( dataString.indexOf('{') !== -1 ) {
                // start collecting as soon as we see { in the buffer
                that.isBufferingJson = true;
                that.jsonBuffer = dataString;
            }

            // sometimes receive json after process has exited
            if ( that.isBufferingJson && that.processExited ) {
                if ( that.jsonBuffer.length ) {
                    setTimeout(function() {
                        that._loadString(that.jsonBuffer,path_);
                    },0);
                }
            }
        });
        that.child.stderr.on("data", function(data) {
            that.didReceiveError = true;
            console.log("ERROR: " + data);
            // deferring seems to allow the spawn to exit cleanly
            setTimeout(function() {
                $(that).trigger("Model:error", {
                    'path':path_, 
                    'description':"Processing failed with error: "+data
                });
            },0);
            // that.cancelLoadFile();
        });
        that.child.on("exit", function(code) {
            console.log("Processing exited with code: " + code);
            // that.child = null;
            // deferring seems to allow the spawn to exit cleanly
            if ( 0 === code ) {
                if ( that.jsonBuffer.length ) {
                    setTimeout(function() {
                        that._loadString(that.jsonBuffer,path_);
                    },0);
                }
            }
            that.processExited = true;
        });
        return that;
    },

    saveFileWithPath:function(path_) {
        var that = this,
            fs = require('fs');

        fs.writeFile(path_, JSON.stringify(that.json, null, 4), function(err) {
            if(err) {
                console.log(err);
            } 
            else {
                console.log("JSON saved to " + path_);
            }
        });
    },

    _sortObject:function(o) {
        var sorted = {},
        key, a = [];

        for (key in o) {
            if (o.hasOwnProperty(key)) {
                a.push(key);
            }
        }

        a.sort();

        for (key = 0; key < a.length; key++) {
            sorted[a[key]] = o[a[key]];
        }
        return sorted;
    },

    _createModels:function() {
        var that = this,
            i,
            o,
            susceptibilityModel,
            virulenceModel,
            calledVariants,
            calledGenes,
            key,
            value,
            templateModel,
            isInducible;

        that.susceptible = [];
        that.resistant = [];
        that.inconclusive = []
        that.positive = [];
        that.negative = []
        that.inducible = [];

        susceptibilityModel = that.json['susceptibility'];


        calledVariants = that.json['called_variants'];
        calledGenes = that.json['called_genes'];

        that.evidenceModel = {};

        for ( key in calledVariants ) {
            mutation = calledVariants[key];
            title = mutation['induced_resistance'];
            genes = key.split('_');
            // group by title
            o = [];
            if ( that.evidenceModel[title] ) {
                o = that.evidenceModel[title];
            }
            else {
                // initialise
                that.evidenceModel[title] = o;
            }
            o.push([
                'Resistance mutation found: '+genes[1]+' in gene '+genes[0],
                'Resistant allele seen '+(mutation['R_median_cov'])+' times',
                'Susceptible allele seen '+(mutation['S_median_cov'])+' times'
            ]);
        }

        for ( key in calledGenes) {
            spot = calledGenes[key];
            title = spot['induced_resistance'];
            // group by title
            o = [];
            if ( that.evidenceModel[title] ) {
                o = that.evidenceModel[title];
            }
            else {
                // initialise
                that.evidenceModel[title] = o;
            }
            o.push([
                key+' gene found',
                'Percent recovered: '+spot['per_cov']+'%',
                'Median coverage: '+spot['median_cov']
            ]);
        }

        that.evidenceModel = that._sortObject(that.evidenceModel);

        // ignore the values
        that.phyloGroup = _.keys(that.json.phylogenetics.phylo_group);

        // build array of included species
        that.species = _.keys(that.json.phylogenetics.species);
                // if ( kTargetSpeciesTB === MykrobeTarget.species ) {
            // sourceSpecies = that.json.phylogenetics.species;
        // }
        // else {
        //     sourceSpecies = that.json.species;
        // }
        // for ( key in sourceSpecies ) {
        //     value = sourceSpecies[key].toLowerCase();
        //     if ( 'major' === value ) {
        //         that.species.push(key);
        //     }
        // }

        that.lineage = [];
        if ( kTargetSpeciesTB === MykrobeTarget.species ) {
            that.lineage = _.keys(that.json.phylogenetics.lineage);
            // sourceLineage = that.json.phylogenetics.lineage;
            // for ( key in sourceLineage ) {
                // value = sourceLineage[key].toLowerCase();
                // if ( 'major' === value ) {
                // that.lineage.push(key);
                // }
            // }
        }

        for ( key in susceptibilityModel ) {
            value = susceptibilityModel[key].substr(0,1).toUpperCase();
            isInducible = susceptibilityModel[key].toUpperCase().indexOf('INDUCIBLE') !== -1;
            if ( 'S' === value ) {
                that.susceptible.push(key);
            }
            else if ( 'R' === value ) {
                that.resistant.push(key);
            }           
            else if ( 'N' === value ) {
                that.inconclusive.push(key);
            }
            if ( isInducible ) {
                that.inducible.push(key);
            }
        }

        if ( 'virulence_toxins' in that.json) {
            virulenceModel = that.json['virulence_toxins'];
            for ( key in virulenceModel ) {
                value = virulenceModel[key].toUpperCase();
                if ( 'POSITIVE' === value ) {
                    that.positive.push(key);
                }       
                else if ( 'NEGATIVE' === value ) {
                    that.negative.push(key);
                }
            }
        }

        templateModel = [
            {
                "class":"susceptible",
                "equalise-height":true,
                "title":"Susceptible",
                "contents":that.susceptible
            },
            {
                "class":"resistant",
                "equalise-height":true,
                "title":"Resistant",
                "contents":that.resistant
            },
            {
                "class":"inconclusive",
                "equalise-height":true,
                "title":"Inconclusive",
                "contents":that.inconclusive
            }
        ];
        if ( !that.inconclusive.length ) {
            templateModel.pop();
        }
        that.viewAllModel = that._viewModelWithTemplate(templateModel);

        // that.viewAllScreenCell.setViewModel(viewAllModel);

        //
        // antibiotic class
        //

        templateModel = [
            {
                "title":"Aminoglycosides",
                "class":"aminoglycosides",
                "equalise-height":true,
                "icon-title":"A",
                "contents":[
                    "Gentamicin"
                ]
            },
            {
                "title":"Beta Lactams",
                "class":"beta-lactams",
                "equalise-height":true,
                "icon-title":"BL",
                "contents":[
                    "Penicillin",
                    "Methicillin"
                ]
            },
            {
                "title":"Macrolides /\rLincosamides",
                "class":"macrolides-lincosamides",
                "equalise-height":true,
                "icon-title":"M/L",
                "contents":[
                    "Clindamycin",
                    "Erythromycin"
                ]
            },
            {
                "title":"Glycopeptides",
                "class":"glycopeptides",
                "equalise-height":true,
                "icon-title":"G",
                "contents":[
                    "Vancomycin"
                ]
            },
            {
                "title":"Tetracyclines",
                "class":"tetracyclines",
                "equalise-height":true,
                "icon-title":"T",
                "contents":[
                    "Tetracycline"
                ]
            },
            {
                "title":"Quinolone",
                "class":"quinolone",                
                "equalise-height":true,
                "icon-title":"Q",
                "contents":[
                    "Ciprofloxacin"
                ]
            },
            {
                "title":"Miscellaneous",
                "class":"miscellaneous",                
                "equalise-height":true,
                "icon-title":"#",
                "contents":[
                    "Rifampicin",
                    "FusidicAcid",
                    "Trimethoprim",
                    "Mupirocin"
                ]
            }
        ];

        that.antibioticClassModel = that._viewModelWithTemplate(templateModel);
        // that.antibioticClassScreenCell.setViewModel(antibioticClassModel);

        //
        // infection site
        //

        templateModel = [
            {
                "title":"Bloodstream",
                "class":"bloodstream",
                "icon-title":"B",
                "contents":[
                    "Penicillin",
                    "Methicillin",
                    "Vancomycin"
                ]
            },
            {
                "title":"Endocarditis",
                "class":"endocarditis",
                "icon-title":"E",
                "contents":[
                    "Penicillin",
                    "Methicillin",
                    "Vancomycin",
                    "Gentamicin",
                    "Rifampicin"
                ]
            },
            {
                "title":"Skin / Soft tissue",
                "class":"skin-soft-tissue",
                "icon-title":"S/T",
                "contents":[
                    "Penicillin",
                    "Methicillin",
                    "Vancomycin",
                    "Clindamycin",
                    "Tetracycline",
                    "Erythromycin"
                ]
            },
            {
                "title":"Joint / Bone",
                "class":"joint-bone",
                "icon-title":"J/B",
                "contents":[
                    "Penicillin",
                    "Methicillin",
                    "Vancomycin",
                    "Tetracycline", 
                    "Erythromycin",
                    "Clindamycin",
                    "Ciprofloxacin",
                    "Rifampicin",
                    "FusidicAcid", 
                    "Trimethoprim"
                ]
            },
            {
                "title":"CNS",
                "class":"cns",
                "icon-title":"C",
                "contents":[
                    "Penicillin",
                    "Methicillin",
                    "Vancomycin"
                ]
            },
            {
                "title":"Lung",
                "class":"lung",
                "icon-title":"L",
                "contents":[
                    "Penicillin",
                    "Methicillin",
                    "Vancomycin",
                    "Tetracycline",
                    "Clindamycin"
                ]
            },
            {
                "title":"MRSA Decolonisation",
                "class":"mrsa-decolonisation",
                "icon-title":"M",
                "contents":[
                    "FusidicAcid",
                    "Mupirocin"
                ]
            }
        ];

        that.infectionSiteModel = that._viewModelWithTemplate(templateModel);

        //
        // virulence
        //

        templateModel = [
            {
                "title":"Virulence markers",
                "class":"virulence-markers",
                "icon-title":"VM",
                "contents":[
                    "PVL",
                    "Toxic shock protein"
                ]
            }
        ];

        that.virulenceMarkersModel = that._viewModelWithTemplate(templateModel);

        //
        // drugs (tb)
        //

        templateModel = [
            {
                "title":"First Line Drugs",
                "class":"first-line-drugs",
                "equalise-height":true,
                "icon-title":"1",
                "contents":[
                    "Isoniazid",
                    "Rifampicin",
                    "Ethambutol",
                    "Pyrazinamide"
                ]
            },
            {
                "title":"Second Line Drugs",
                "class":"second-line-drugs",
                "equalise-height":true,
                "icon-title":"2",
                "contents":[
                    "Quinolones",
                    "Streptomycin",
                    "Amikacin",
                    "Capreomycin",
                    "Kanamycin"
                ]
            }
        ];

        that.drugsModel = that._viewModelWithTemplate(templateModel);

        /*
        Here's a new thing - if R for both Isoniazid  and Rifampicin, then mark/flag as MDR (Multi-Drug Resistant)
        So we'll need to write MDR somewhere
        */

        that.drugsResistanceModel = {
            mdr:false,
            xdr:false
        };

        if ( that.resistant.indexOf('Isoniazid') !== -1 && that.resistant.indexOf('Rifampicin') !== -1 ) {
            that.drugsResistanceModel.mdr = true;   
            /*
            If MDR AND R to both fluoroquinolones and one of the other these 3 (Amikacin, Kanamycin, Capreomycin), then call it XDR (Extensively Drug Resistant) 
            */
            if ( that.resistant.indexOf('Quinolones') ) {
                if ( that.resistant.indexOf('Amikacin') !== -1 || that.resistant.indexOf('Kanamycin') !== -1 || that.resistant.indexOf('Capreomycin') !== -1 ) {
                    that.drugsResistanceModel.xdr = true;
                }
            }
        }

        if ( kTargetSpeciesTB === MykrobeTarget.species ) {
            that.speciesModel = that.species.join(' / ') +' (lineage: '+that.lineage[0]+')';
        }
        else {
            that.speciesModel = that.species.join(' / ');
        }
        return that;
    },

    // replaces each name in the 'contents' with a pair of 'title' and 'class' indicating its state

    _viewModelWithTemplate:function(templateModel_) {
        var that = this,
            model = [],
            contentsModel,
            i,
            o,
            title,
            gene,
            spot,
            mutation,
            mutationDetails,
            found;

        for ( i = 0; i < templateModel_.length; i++ ) {
            o = templateModel_[i];
            contentsModel = [];
            // display in the order provided
            // o['contents'].sort();
            for ( j = 0; j < o['contents'].length; j++ ) {
                title = o['contents'][j];
                found = true;
                c = {
                    "title":title
                }
                if ( that.susceptible.indexOf(title) != -1) {
                    c['class'] = "susceptible";
                    c['tooltip'] = "Susceptible";
                }
                else if ( that.resistant.indexOf(title) != -1) {
                    c['class'] = "resistant";
                    c['tooltip'] = "Resistant";
                }
                else if ( that.inconclusive.indexOf(title) != -1) {
                    c['class'] = "inconclusive";
                    c['tooltip'] = "Inconclusive";
                }
                else if ( that.positive.indexOf(title) != -1) {
                    c['class'] = "positive";
                    c['tooltip'] = "Positive";
                }
                else if ( that.negative.indexOf(title) != -1) {
                    c['class'] = "negative";
                    c['tooltip'] = "Negative";
                }
                else {
                    console.log('not found in results:'+title);
                    found = false;
                }
                if ( found ) {
                    if ( that.inducible.indexOf(title) != -1 ) {
                        c['class'] += " inducible";
                        c['tooltip'] += " (inducible)";
                    }
                    // search through calledVariantsModel
                    mutationDetails = [];
                    for ( gene in that.calledVariants ) {
                        mutation = that.calledVariants[gene];
                        console.log('gene:'+gene);
                        console.log(mutation);
                        if ( mutation['induced_resistance'] == title ) {
                            details = (gene.replace('_',' ') + ' S('+mutation['S_per_cov']+'), R('+mutation['R_per_cov']+')');
                            mutationDetails.push(details);
                        }
                    }
                    if ( mutationDetails.length ) {
                        c['mutation'] = mutationDetails;
                    }

                    geneDetails = [];
                    for ( gene in that.calledGenes ) {
                        spot = that.calledGenes[gene];
                        console.log('gene:'+gene);
                        console.log(spot);
                        if ( spot['induced_resistance'] == title ) {
                            details = (gene.replace('_',' ') + ' '+spot['per_cov']+'%, '+spot['median_cov']);
                            geneDetails.push(details);
                        }
                    }
                    if ( geneDetails.length ) {
                        c['gene'] = geneDetails;
                    }
                    contentsModel.push(c);
                }
            }
            o.contents = contentsModel;
            model.push(o);
        }

        console.log(JSON.stringify(model, null, 4));
        return model;
    }

});