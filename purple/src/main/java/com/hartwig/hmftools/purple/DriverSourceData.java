package com.hartwig.hmftools.purple;

import static java.lang.String.format;

import com.hartwig.hmftools.common.driver.DriverCatalog;

public class DriverSourceData
{
    public final DriverCatalog DriverData;
    public final Object SourceObject;

    public DriverSourceData(final DriverCatalog driverData, final Object sourceObject)
    {
        DriverData = driverData;
        SourceObject = sourceObject;
    }

    public String toString()
    {
        return format("gene(%s) type(%s) likelihood(%.2f) location(%s:%s)",
                DriverData.gene(), DriverData.driver(), DriverData.driverLikelihood(), DriverData.chromosome(), DriverData.chromosomeBand());
    }
}
