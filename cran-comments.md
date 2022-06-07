# Re-submission of VIEWpoly (06-07-2022)

## Maintainer comments

Dear CRAN Team,

This is a re-submission of VIEWpoly package. In this version (0.2.0):

- Broken links were fixed
- Download of the images with .tiff format was fixed
- Avoids warnings and error messages displayed in the console during app execution 
- Error messages were improved
- Function documentation were improved
- Aesthetic improvements
- App is also available on shinyapps.io: https://cris-taniguti.shinyapps.io/viewpoly/
- Functional testing added
- Disable the download buttons when image parameters are not reliable
- Title of boxes are now also collapsible

Thank you for reviewing our re-submission!

## Test environments 

The package passed through several checks (with the flag --as-cran) in the following systems:

- Local (Windows 11 v10.0.22000 Build 22000)
- Winbuilder (Windows Server 2008 R2 SP1, R-devel, 32/64 bit)
- GitHub actions (windows-latest release)
- Github actions (macOS-latest release)
- Github actions (ubuntu-20.04 release)
- Github actions (ubuntu-20.04 devel)
- R-hub (Ubuntu Linux 20.04.1 LTS, R-release, GCC)
- R-hub (Fedora Linux, R-devel, clang, gfortran)

## R CMD check results

0 errors | 0 warnings | 1 note

There are one NOTE that appeared when running some checks. One of the directories containing example files exceeds 1Mb. It has a total of 4.5Mb. The examples are important to users explore the app features and file format before input their own data.

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.