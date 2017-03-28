
for /f %%f in ('dir *.exr /b') do exrtotiff.exe %%f %%f.tiff