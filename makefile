gcc -o repeatFinderDNA repeatFinderDNA.c -lm
gcc -o repeatFinderProtein repeatFinderProtein.c -lm
chmod 755 repeatFinderDNA && chown www-data:root repeatFinderDNA
chmod 755 repeatFinderProtein && chown www-data:root repeatFinderProtein