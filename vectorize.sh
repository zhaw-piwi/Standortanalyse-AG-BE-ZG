

# conda activate geopython

# try it on one raster file

# reclassify into integers:
gdal_calc.py -A data-out/HI-historic/1974-1983.tif --calc="(A<=1400)*1 + (A>1400)*(A<=2100)*2 + (A>2100)*3" --outfile data-temp/1974-1983.tif
# polygonize to geopackage
gdal_polygonize.py data-temp/1974-1983.tif -overwrite data-out/HI_sf.gpkg years1974-1983 type 


rm data-out/HI_sf.gpkg

for raster in find data-out/HI-historic/*.tif ; 
do 
    base_name=$(basename ${raster})
    tempraster=data-temp/$base_name
    year=${base_name%.*}
    echo $base_name
    gdal_calc.py -A $raster --calc="(A<=1400)*1 + (A>1400)*(A<=2100)*2 + (A>2100)*3" --outfile $tempraster
    gdal_polygonize.py $tempraster -overwrite data-out/HI_sf.gpkg years_$year type 

done


rm -r data-temp
