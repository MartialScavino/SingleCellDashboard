# Single-cell Dashboard
Shiny dashboard to perform single cell analyses


# Launch the app

Clone this repository, open the `app.R` file and launch every lines of the script.

The `source("R/Dependencies.R")`command will install all the R packages needed.  
If the installation failed for a package, you will be noticed in the R console. Try installing it manually.

For windows users, the installation of the Ucell package is likely to malfunction. 
If this is the case, please try launching this command before : `options(download.file.method = "wininet")`
If it doesn't solve the issue, you can still launch the application with the `runApp('.')` command but the Signature tab will not work properly.

If you have any question or issue, please contact me at : martial.scavino@etu.univ-lyon1.fr
