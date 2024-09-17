REM The following batch file will concatenate all .h and .cpp files

REM To have .h and .cpp pairs next to each other, go through all files and 
REM add an if to check for either .h or .cpp and if true, output file
REM Have not implemented, but should work

echo "Beginning Concatenation of Files"

if exist LISTING.cpp del LISTING.cpp
if exist LISTING.ppx del LISTING.ppx


for %%f in (*.h *.cpp) do (
	echo //############################################################################>> .\LISTING.ppx
	echo //%%f >> .\LISTING.ppx
	echo //############################################################################>> .\LISTING.ppx
	type "%%f"  >> .\LISTING.ppx 
	echo.  >> .\LISTING.ppx 
	echo. >> .\LISTING.ppx 

	echo. >> .\LISTING.ppx
)

ren LISTING.ppx LISTING.cpp

echo "Finishing Concatenation of Files"
