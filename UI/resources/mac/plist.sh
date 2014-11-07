#!/bin/bash
app="$1"

/usr/libexec/PlistBuddy -c "Delete :CFBundleDocumentTypes" "$app/Contents/Info.plist"

/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes array" "$app/Contents/Info.plist"
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions array" "$app/Contents/Info.plist"
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions:0 string json" "$app/Contents/Info.plist"
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions:1 string bam" "$app/Contents/Info.plist"
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions:2 string gz" "$app/Contents/Info.plist"
/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeExtensions:3 string fastq" "$app/Contents/Info.plist"

/usr/libexec/PlistBuddy -c "Add :CFBundleDocumentTypes:0:CFBundleTypeRole string Viewer" "$app/Contents/Info.plist"
