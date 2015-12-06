wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_1997/WL_ALL_1997.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_1998/WL_ALL_1998.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_1999/WL_ALL_1999.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2000/WL_ALL_2000.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2001/WL_ALL_2001.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2002/WL_ALL_2002.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2003/WL_ALL_2003.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2004/WL_ALL_2004.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2005/WL_ALL_2005.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2006/WL_ALL_2006.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2007/WL_ALL_2007.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2008/WL_ALL_2008.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2009/WL_ALL_2009.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2010/WL_ALL_2010.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2011/WL_ALL_2011.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2012/WL_ALL_2012.zip
wget -c http://ne.water.usgs.gov/ogw/hpwlms/wldata/d_2013/WL_ALL_2013.zip

R --no-save --vanilla --slave < calculate_saturated_thickness.R
