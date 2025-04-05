Base.getindex(ds::AbstractDataset, name::Symbol) = DatasetVariable(ds, name)

load(dv::DatasetVariable; kwargs...) = load(dv.ds, dv.name; kwargs...)
download(dv::DatasetVariable; kwargs...) = download(dv.ds, dv.name; kwargs...)

function copernicus_marine_resolve(
    product_id, dataset_id;
    asset_name = "timeChunked",
    catalog_url = "https://stac.marine.copernicus.eu/metadata/catalog.stac.json")
    # This function is used to resolve remote URLs.
    cat = STAC.Catalog(catalog_url)
    list_item = collect(keys(cat[product_id].items))
    item_candidates = sort(filter(startswith(dataset_id), list_item))
    # Use the last version by default
    dataset_version_id = item_candidates[end]
    item = cat[product_id].items[dataset_version_id]
    return STAC.href(item.assets[asset_name])
end

"""
    ds = ROMS.CMEMS_zarr(product_id, mapping, cachedir;
                         chunks = 60,
                         time_shift = 0,
                         source = :local,
                         kwargs...)

Returns a structure `ds` that connects to a CMEMS dataset.

`mapping` is a dictionary linking CF standard variable names to the underlying identifiers.
If a mapping entry is, for example,
    :sea_surface_height_above_geoid => "so",
then when `source == :local` the function looks for a file named "so.nc" in `cachedir`.
If `source == :remote`, the function resolves a URL using `copernicus_marine_resolve`.

Arguments:
  - `cachedir`: directory where files are cached locally (or a placeholder when using remote URLs)
  - `source`: either `:local` to use local NetCDF files or `:remote` to use remote URLs.
  - `chunks`: chunking parameter for the underlying dataset.
  - `time_shift`: an integer time shift (in seconds) for the dataset.
"""
function CMEMS_zarr(product_id, mapping, cachedir;
                    chunks = 60,
                    time_shift = 0,
                    source = :local,
                    kwargs...)
    files = DefaultDict{Symbol, String, String}("unknown")
    for (k, v) in mapping
        dataset_id = (length(v) > 1 && v isa Tuple) ? v[end] : v
        if source == :local
            # Build a local file path (e.g. "so.nc") from cachedir
            files[k] = joinpath(cachedir, string(dataset_id, ".nc"))
        elseif source == :remote
            # Resolve the remote URL
            files[k] = copernicus_marine_resolve(product_id, dataset_id; kwargs...)
        else
            error("Unknown source option: $source. Valid options are :local or :remote.")
        end
    end

    # Choose the dataset type based on source
    if source == :local
        dataset_constructor = NCDataset
    elseif source == :remote
        dataset_constructor = ZarrDataset
    else
        error("Unknown source option: $source. Valid options are :local or :remote.")
    end

    dataset = CDMDataset(
        dataset_constructor,
        files,
        cachedir = cachedir,
        options = Dict(:_omitcode => [404, 403]),
        time_shift = time_shift,
        chunks = chunks)

    return dataset
end
