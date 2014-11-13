'use strict';
module.exports = function(grunt) {
    // Load grunt tasks automatically
    require('load-grunt-tasks')(grunt);
    // configurable paths
    var path = require('path'),
    config = {
        app: 'app',
        dist: 'dist',
        tmp: '.tmp',
        // icon: './icon/mac',
        resources: 'resources'
    };
    grunt.initConfig({
        config: config,
        pkg: grunt.file.readJSON(path.join(config.app,'package.json')),
        // Project settings
        clean: {
            dist: {
                files: [{
                    dot: true,
                    src: ['<%= config.dist %>', '<%= config.tmp %>', '<%= nodewebkit.options.build_dir %>', '.grunt', '.sass-cache', '_SpecRunner.html']
                }]
            }
        },
        // Automatically inject Bower components into the HTML file
        bowerInstall: {
            app: {
                src: ['<%= config.app %>/index.html', '<%= config.app %>/about.html'],
                ignorePath: '<%= config.app %>/'
            }
        },
        jshint: {
            options: {
                jshintrc: '.jshintrc'
            },
            files: '<%= config.app %>/js/main.js'
        },
        nodewebkit: {
            options: {
                build_dir: './build', // Where the build version of my node-webkit app is saved
                // credits: './public/credits.html',
                mac_icns: './resources/icon/<%= pkg.targetName %>/icon.icns',
                mac: true,
                win: true,
                linux32: false,
                linux64: false,
            },
            src: '<%= config.dist %>/**/*'
        },
        // Watches files for changes and runs tasks based on the changed files
        watch: {
            js: {
                files: ['<%= config.app %>/js/{,*/}*.js'],
                // tasks: ['jshint'],
                options: {
                    livereload: true
                }
            },
            gruntfile: {
                files: ['Gruntfile.js']
            },
            sass: {
                files: ['<%= config.app %>/css/{,*/}*.{scss,sass}'],
                tasks: ['sass']
            },
            styles: {
                files: ['<%= config.app %>/css/{,*/}*.css'],
                tasks: ['copy:styles']
            },
            livereload: {
                options: {
                    livereload: '<%= connect.options.livereload %>'
                },
                files: ['<%= config.app %>/{,*/}*.html', '.tmp/css/{,*/}*.css', '<%= config.app %>/img/{,*/}*.{gif,jpeg,jpg,png,svg,webp}']
            }
        },
        // The actual grunt server settings
        connect: {
            options: {
                port: 9000,
                livereload: 35729,
                // Change this to '0.0.0.0' to access the server from outside
                hostname: 'localhost'
            },
            livereload: {
                options: {
                    open: true,
                    base: ['.tmp', '<%= config.app %>']
                }
            }
        },
        sass: {
            options: {
                // sassDir: '<%= config.app %>/styles',
                // cssDir: '.tmp/styles',
                // generatedImagesDir: '.tmp/images/generated',
                // imagesDir: '<%= config.app %>/images',
                // javascriptsDir: '<%= config.app %>/scripts',
                // fontsDir: '<%= config.app %>/styles/fonts',
                // importPath: '<%= config.app %>/bower_components',
                // httpImagesPath: '/images',
                // httpGeneratedImagesPath: '/images/generated',
                // httpFontsPath: '/styles/fonts',
                // relativeAssets: false,
                // assetCacheBuster: false
                // expand: true
            },
            dist: {
                files: [{
                    expand: true,
                    cwd: '<%= config.app %>/css/',
                    // src: '{,*/}*.scss',
                    src: 'main.scss',
                    dest: '<%= config.app %>/css/',
                    ext: '.css'
                }]
            }
        },
        exec: {
            plist: {
                cmd: function() {
                    var path = require('path'),
                        config = grunt.config.get(['config']),
                        nodewebkit = grunt.config.get(['nodewebkit']),
                        packageJson = grunt.file.readJSON(config.app + '/package.json'),
                        appPath,
                        command;
                    appPath = path.join(nodewebkit.options.build_dir, packageJson.name, 'osx', packageJson.name + '.app');
                    command = 'resources/mac/plist.sh "' + appPath + '"';
                    grunt.log.writeln(command);
                    return command;
                }
            },

            sign: {
                cmd: function() {
                    var path = require('path'),
                        config = grunt.config.get(['config']),
                        nodewebkit = grunt.config.get(['nodewebkit']),
                        packageJson = grunt.file.readJSON(config.app + '/package.json'),
                        appPath,
                        command;
                    appPath = path.join(nodewebkit.options.build_dir, packageJson.name, 'osx', packageJson.name + '.app');
                    command = 'resources/mac/sign.sh "' + appPath + '"';
                    grunt.log.writeln(command);
                    return command;
                }
            },

            dmg: {
                cmd: function() {
                    var path = require('path'),
                        config = grunt.config.get(['config']),
                        nodewebkit = grunt.config.get(['nodewebkit']),
                        packageJson = grunt.file.readJSON(config.app + '/package.json'),
                        appPath,
                        command;
                    appPath = path.join(nodewebkit.options.build_dir, packageJson.name, 'osx', packageJson.name + '.app');
                    command = 'resources/mac/package.sh "' + packageJson.name + '" "'+appPath+'" "' + packageJson.targetName + '"';
                    grunt.log.writeln(command);
                    return command;
                }
            }
        },
        copy: {
            dist: {
                files: [{
                    expand: true,
                    dot: true,
                    cwd: '<%= config.app %>',
                    dest: '<%= config.dist %>',
                    src: [
                        '*.{ico,png,txt}',
                        'img/**',
                        '{,*/}*.html',
                        'css/*.css',
                        'css/webfonts/{,*/}*.*',
                        'js/vendor/swiffy/runtime.js',
                        'bin/<%= pkg.targetName %>/**'
                    ]
                }]
            },
            styles: {
                expand: true,
                dot: true,
                cwd: '<%= config.app %>/css',
                dest: '<%= config.tmp %>/css/',
                src: '{,*/}*.css'
            }
        },
        // Reads HTML for usemin blocks to enable smart builds that automatically
        // concat, minify and revision files. Creates configurations in memory so
        // additional tasks can operate on them
        useminPrepare: {
            options: {
                dest: '<%= config.dist %>'
            },
            html: ['<%= config.app %>/index.html']
        },

        // Performs rewrites based on rev and the useminPrepare configuration
        usemin: {
            options: {
                assetsDirs: ['<%= config.dist %>', '<%= config.dist %>/img', '<%= config.dist %>/css/webfonts']
            },
            html: ['<%= config.dist %>/{,*/}*.html'],
            css: ['<%= config.dist %>/css/{,*/}*.css']
        },

        uglify: {
            options: {
                compress: {
                    'drop_console':true
                }
            }
        },

        htmlmin: {
            dist: {
                options: {
                    collapseBooleanAttributes: true,
                    collapseWhitespace: false,
                    removeAttributeQuotes: false,
                    removeComments: true,
                    removeCommentsFromCDATA: true,
                    removeEmptyAttributes: true,
                    removeOptionalTags: false,
                    removeRedundantAttributes: false,
                    useShortDoctype: true
                },
                files: [{
                    expand: true,
                    cwd: '<%= config.app %>',
                    src: '{,*/}*.html',
                    dest: '<%= config.dist %>'
                }]
            }
        },
        jasmine: {
            mykrobe: {
                src: [
                    '<%= config.app %>/bower_components/jquery/dist/jquery.js',
                    '<%= config.app %>/js/plugins.js',
                    '<%= config.app %>/js/mykrobe/Model.js'
                ],
                options: {
                    specs: 'spec/*Spec.js',
                    helpers: 'spec/*Helper.js',
                    keepRunner:true
                }
            }
        }
    });

    // generate mac pngs from pdf source
    // then convert into iconset
    /*
    icon_16x16.png
    icon_16x16@2x.png
    icon_32x32.png
    icon_32x32@2x.png
    icon_128x128.png
    icon_128x128@2x.png
    icon_256x256.png
    icon_256x256@2x.png
    icon_512x512.png
    icon_512x512@2x.png
    */
    //http://stackoverflow.com/questions/653380/converting-a-pdf-to-png
    grunt.registerTask('mac-icons', 'Create mac icons', function() {
        var path = require('path'),
            config = grunt.config.get(['config']),
            packageJson = grunt.file.readJSON(path.join(config.app,'package.json'));
        var iconPath = path.join('resources','icon',packageJson['targetName']);
        var sizes = [16, 32, 128, 256, 512];
        var commands = [];
        var counter = 0;
        for (var i = 0; i < sizes.length; i++) {
            var size = sizes[i];
            var sizeRetina = size * 2;
            var filenameBase = 'icon_' + size + "x" + size;
            var filename = filenameBase + '.png';
            var filenameRetina = filenameBase + '@2x.png';  
            // "/Library/Application Support/Adobe/Color/Profiles/Recommended/sRGB Color Space Profile.icm"          
            // var command = 'sips --matchTo "/Library/Application Support/Adobe/Color/Profiles/Recommended/sRGB Color Space Profile.icm" -s format png -z ' + size + ' ' + size + ' --out "'+ path.join(iconPath, 'icon.iconset', filename) + '" "' + path.join(iconPath, 'icon.pdf')+'"';
            // var commandRetina = 'sips --matchTo "/Library/Application Support/Adobe/Color/Profiles/Recommended/sRGB Color Space Profile.icm" -s format png -z ' + sizeRetina + ' ' + sizeRetina + ' --out "'+ path.join(iconPath, 'icon.iconset', filenameRetina) + '" "' + path.join(iconPath, 'icon.pdf')+'"';
            
            var command = 'gs -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r72 -dPDFFitPage=true -g' + size + 'x' + size + ' -sOutputFile="'+ path.join(iconPath, 'icon.iconset', filename) + '" "' + path.join(iconPath, 'icon.pdf')+'"';
            var commandRetina = 'gs -dNOPAUSE -dBATCH -sDEVICE=pngalpha -r72 -dPDFFitPage=true -g' + sizeRetina + 'x' + sizeRetina + ' -sOutputFile="'+ path.join(iconPath, 'icon.iconset', filenameRetina) + '" "' + path.join(iconPath, 'icon.pdf')+'"';           

            grunt.log.writeln(command);
            commands['exec' + counter] = {
                command: command
            };
            counter++;
            commands['exec' + counter] = {
                command: commandRetina
            };
            counter++;
        }
        commands['exec' + counter] = {
            command: 'iconutil -c icns -o "' + path.join(iconPath, 'icon.icns')+'" "' + path.join(iconPath, 'icon.iconset')+'"'
        };
        grunt.config('exec', commands);
        grunt.task.run('exec');
    });
    
    // set options for distribution package.json
    // typically set window.toolbar to false
    grunt.registerTask('packagejsondist', function() {
        var path = require('path'),
            config = grunt.config.get(['config']),
            packageJson = grunt.file.readJSON(config.app + '/package.json');
        packageJson.window.toolbar = false;
        grunt.file.write(path.join(config.dist,'package.json'), JSON.stringify(packageJson,null,2), {
            encoding: 'UTF8'
        });
    });

    grunt.registerTask('set-target', function() {
        var path = require('path'),
            config = grunt.config.get(['config']),
            packageJson = grunt.file.readJSON(path.join(config.app,'package.json')),
            choices = grunt.file.readJSON('targets.json'),
            defaultTargetName = packageJson['targetName'] ? packageJson['targetName'] : choices[0].value;

        grunt.config('prompt', {
            'target': {
                'options': {
                    'questions': [
                        {
                            config: 'targetName',
                            type: 'list',
                            message: 'Select target app',
                            default: defaultTargetName,
                            choices: choices
                        }
                    ],
                    then:function(results){
                        // grunt.log.writeln('hello!');
                        // grunt.log.writeln(JSON.stringify(results));
                        var targetName = results['targetName'];
                        var appName = '';
                        var appVersion = '';
                        var appDisplayVersion = '';
                        for ( var i = 0; i < choices.length; i++ ) {
                            if ( choices[i].value === targetName ) {
                                appName = choices[i].name;
                                appVersion = choices[i].version;
                                appDisplayVersion = choices[i].displayVersion;
                            }
                        }
                        packageJson['targetName'] = targetName;
                        packageJson.name = appName;
                        packageJson.version = appVersion;
                        packageJson.displayVersion = appDisplayVersion;
                        packageJson.window.title = appName;
                        grunt.file.write(path.join(config.app,'package.json'), JSON.stringify(packageJson,null,2), {
                            encoding: 'UTF8'
                        });
                    }
                }
            }
        });
        grunt.task.run('prompt');
        // var validation = grunt.config('validation');
        // grunt.log.writeln('that.targetName:'+that.targetName);
    });

    // add bower packages to html head blocks
    grunt.registerTask('bower-install', ['bowerInstall']);

    // minify, generate app
    grunt.registerTask('test', [
        'jasmine'
    ]);

    // minify assets for final packaging, can be tested using 'nodewebkit dist'
    grunt.registerTask('minify', [
        'sass:dist',
        'copy:dist',
        'useminPrepare',
        'concat',
        'cssmin',    
        'uglify',
        'usemin'
    ]);
    // minify, generate app
    grunt.registerTask('build', [
        'clean', 
        'minify',
        'packagejsondist',
        'nodewebkit'
    ]);
    // sign, create mac dmg
    grunt.registerTask('dist', [
        'build', 
        'exec:plist', 
        'exec:sign', 
        'exec:dmg'
    ]);
};