function interp_clim(domain, clim_filename, dataset, timerange;
    padding = 0.5,
    missing_value = -9999.)

    x, y, z_r = domain.lon, domain.lat, domain.z_r
    angle, h, h_u, h_v = domain.angle, domain.h, domain.h_u, domain.h_v

    function wider(x)
        xmin, xmax = extrema(x)
        return xmin - padding, xmax + padding
    end

    query = (
        time = timerange,
        longitude = wider(x),
        latitude = wider(y),
    )

    vars = [
        :sea_surface_height_above_geoid,
        :sea_water_potential_temperature,
        :sea_water_salinity,
        :eastward_sea_water_velocity,
        :northward_sea_water_velocity
    ]

    loaded = Dict{Symbol,Any}()
    coords = Dict{Symbol,Any}()

    for var in vars
        cachefile = joinpath(dataset.cachedir, "cached_$(var)_$(Dates.format(timerange[1], "yyyymmdd"))_$(Dates.format(timerange[end], "yyyymmdd")).nc")

        if isfile(cachefile)
            @info "Loading $var from cache: $cachefile"
            ds = NCDataset(cachefile, "r")
            v = ds[string(var)][:]
            loaded[var] = v
            # Dummy coordinates; replace with real extraction logic if needed
            coords[var] = (x, y, z_r, timerange)
            close(ds)
        else
            @info "Downloading $var"
            v, crd = load(dataset[var]; query...)
            loaded[var] = v
            coords[var] = crd

            @info "Saving $var to cache: $cachefile"
            ds = NCDataset(cachefile, "c")
            ds_def = defVar(ds, string(var), eltype(v), size(v))
            ds_def[:] = v
            close(ds)
        end
    end

    zeta, (zx, zy, zt) = loaded[:sea_surface_height_above_geoid], coords[:sea_surface_height_above_geoid]
    temp, (tx, ty, tz, tt) = loaded[:sea_water_potential_temperature], coords[:sea_water_potential_temperature]
    salt, (sx, sy, sz, st) = loaded[:sea_water_salinity], coords[:sea_water_salinity]
    u, (ux, uy, uz, ut) = loaded[:eastward_sea_water_velocity], coords[:eastward_sea_water_velocity]
    v, (vx, vy, vz, vt) = loaded[:northward_sea_water_velocity], coords[:northward_sea_water_velocity]

    angle = repeat(domain.angle, inner = (1, 1, size(z_r, 3)))

    time = st
    N = length(st)

    ds = ROMS.def_clim(clim_filename, missing_value, size(x, 1), size(x, 2), size(z_r, 3))

    climtime, czeta, cubar, cvbar, cu, cv, ctemp, csalt = ds["time"], ds["zeta"], ds["ubar"], ds["vbar"], ds["u"], ds["v"], ds["temp"], ds["salt"]

    for ni = 1:N
        zz = zeros(size(zeta, 1), size(zeta, 2))
        t = time[ni]
        @info "load $t"
        zv = nomissing(zeta[:, :, ni], NaN)
        sv = nomissing(salt[:, :, :, ni], NaN)
        tv = nomissing(temp[:, :, :, ni], NaN)
        uv = nomissing(u[:, :, :, ni], NaN)
        vv = nomissing(v[:, :, :, ni], NaN)

        @info "interpolate $t"
        zetai = ROMS.model_interp3(zx, zy, zz, zv, x, y, z_r[:, :, end], missing = :ufill)
        salti = ROMS.model_interp3(sx, sy, sz, sv, x, y, z_r, missing = :ufill)
        tempi = ROMS.model_interp3(tx, ty, tz, tv, x, y, z_r, missing = :ufill)

        ui_rc = ROMS.model_interp3(ux, uy, uz, uv, x, y, z_r, missing = :zero)
        vi_rc = ROMS.model_interp3(vx, vy, vz, vv, x, y, z_r, missing = :zero)

        ui_r = cos.(angle) .* ui_rc + sin.(angle) .* vi_rc
        vi_r = -sin.(angle) .* ui_rc + cos.(angle) .* vi_rc

        ui = (ui_r[1:end-1, :, :] + ui_r[2:end, :, :]) / 2
        vi = (vi_r[:, 1:end-1, :] + vi_r[:, 2:end, :]) / 2

        U, = ROMS.vinteg(uv, uz)
        V, = ROMS.vinteg(vv, vz)

        Ui_rc = ROMS.model_interp3(ux, uy, uz, U, x, y, z_r[:, :, end], missing = :zero)
        Vi_rc = ROMS.model_interp3(vx, vy, vz, V, x, y, z_r[:, :, end], missing = :zero)

        Ui_r = cos.(domain.angle) .* Ui_rc + sin.(domain.angle) .* Vi_rc
        Vi_r = -sin.(domain.angle) .* Ui_rc + cos.(domain.angle) .* Vi_rc

        Ui = (Ui_r[1:end-1, :, :] + Ui_r[2:end, :, :]) / 2
        Vi = (Vi_r[:, 1:end-1, :] + Vi_r[:, 2:end, :]) / 2

        ubar = Ui ./ h_u
        vbar = Vi ./ h_v

        ubar2, vbar2 = ROMS.vavg(domain, ui, vi)

        ui .+= (ubar - ubar2)
        vi .+= (vbar - vbar2)

        ubar2c, vbar2c = ROMS.vavg(domain, ui, vi)

        climtime[ni] = t
        czeta[:, :, ni] = zetai
        csalt[:, :, :, ni] = salti
        ctemp[:, :, :, ni] = tempi
        cu[:, :, :, ni] = ui
        cv[:, :, :, ni] = vi
        cubar[:, :, ni] = ubar
        cvbar[:, :, ni] = vbar
    end

    close(ds)
end
