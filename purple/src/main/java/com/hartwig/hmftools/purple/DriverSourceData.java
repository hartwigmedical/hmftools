package com.hartwig.hmftools.purple;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;

public class DriverSourceData
{
    public final DriverCatalog DriverData;
    public final Object SourceObject;

    public DriverSourceData(final DriverCatalog driverData, final Object sourceObject)
    {
        DriverData = driverData;
        SourceObject = sourceObject;
    }
}
