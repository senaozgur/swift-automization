from astropy.io import fits
# This is for to check pile up problem. If counts/s> 150 for WT the fit can not be done well.
# Open your source spectral file (.pha)
hdul = fits.open("/cdata1/senaozgur/swift/MAXI_J1348_630/analysis/00011107007-xrt/results/model_soft_35/sw00011107007_grp.pha")
header = hdul[1].header

print("Exposure Time:", header.get("EXPOSURE"))
print("Total Counts:", header.get("TOTCTS"))
print("Count Rate:", header.get("TOTCTS") / header.get("EXPOSURE"))