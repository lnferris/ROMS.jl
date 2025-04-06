using NCDatasets
using Interpolations

"""
    refine_resolution(parent_grid_file::String, child_grid_file::String, dx::Float64, dy::Float64)

Interpolates the parent ROMS grid to a specific uniform spacing (dx, dy) and writes the refined grid to a new NetCDF file.

# Arguments
- `parent_grid_file::String`: Path to the parent grid NetCDF file.
- `child_grid_file::String`: Path to the output child grid NetCDF file.
- `dx::Float64`: Desired uniform grid spacing in the x-direction.
- `dy::Float64`: Desired uniform grid spacing in the y-direction.
"""
function refine_res(parent_grid_file::String, child_grid_file::String, dx::Float64, dy::Float64)
    # Open the parent grid file
    ds_parent = NCDataset(parent_grid_file, "r")

    # Read parent grid variables
    lon_rho = ds_parent["lon_rho"][:,:]
    lat_rho = ds_parent["lat_rho"][:,:]
    h = ds_parent["h"][:,:]

    # Determine the bounds of the parent grid
    lon_min, lon_max = extrema(lon_rho)
    lat_min, lat_max = extrema(lat_rho)

    # Create new uniform grid coordinates
    lon_new = lon_min:dx:lon_max
    lat_new = lat_min:dy:lat_max

    # Create interpolation functions for h, lon_rho, and lat_rho
    itp_h = interpolate((lon_rho[:,1], lat_rho[1,:]), h, Gridded(Linear()))
    itp_lon = interpolate((lon_rho[:,1], lat_rho[1,:]), lon_rho, Gridded(Linear()))
    itp_lat = interpolate((lon_rho[:,1], lat_rho[1,:]), lat_rho, Gridded(Linear()))

    # Initialize arrays for the new grid
    h_new = [itp_h(lon, lat) for lon in lon_new, lat in lat_new]
    lon_rho_new = [itp_lon(lon, lat) for lon in lon_new, lat in lat_new]
    lat_rho_new = [itp_lat(lon, lat) for lon in lon_new, lat in lat_new]

    # Close the parent dataset
    close(ds_parent)

    # Create the child grid NetCDF file
    ds_child = NCDataset(child_grid_file, "c")

    # Define dimensions
    defDim(ds_child, "xi_rho", length(lon_new))
    defDim(ds_child, "eta_rho", length(lat_new))

    # Define variables and write data
    vlon = defVar(ds_child, "lon_rho", Float64, ("xi_rho", "eta_rho"))
    vlat = defVar(ds_child, "lat_rho", Float64, ("xi_rho", "eta_rho"))
    vh = defVar(ds_child, "h", Float64, ("xi_rho", "eta_rho"))

    vlon[:,:] = lon_rho_new
    vlat[:,:] = lat_rho_new
    vh[:,:] = h_new

    # Add global attributes
    attr!(ds_child, "refine_coef", dx)  # Assuming dx = dy for uniform refinement

    # Close the child dataset
    close(ds_child)

    println("Child grid created with uniform spacing dx=$dx, dy=$dy and saved to $child_grid_file")
end
