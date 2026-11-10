def find_coord(ds, std_name):
    """This procedure returns the name of the NetCDF coordinate with
    standard name `std_name`, if such a coordinate exists. The
    standard name is only used to know what to search, it is not used
    for the search itself. The search itself is done via a string
    match on the attribute `units`. So the NetCDF variable one looks
    for does not need to have the attribute `std_name`.

    Adapted from NetCDF95.

    """

    if std_name == "longitude":
        required_units = ["degrees_east"]
        exact = True
    elif std_name == "latitude":
        required_units = ["degrees_north"]
        exact = True
    elif std_name == "time":
        required_units = [" since"]
        exact = False
    elif std_name == "plev":
        required_units = ["Pa", "hPa", "millibar", "mbar"]
        exact = True
    else:
        raise ValueError("bad value of std_name")

    found = False

    for my_dim in ds.dimensions:
        if (
            my_dim in ds.variables
            and ds.variables[my_dim].ndim == 1
            and ds.variables[my_dim].dimensions[0] == my_dim
            and "units" in ds.variables[my_dim].ncattrs()
        ):
            # Found a coordinate
            if exact:
                found = ds.variables[my_dim].units in required_units
            else:
                found = required_units[0] in ds.variables[my_dim].units

            if found:
                break
    else:
        my_dim = None

    return my_dim
