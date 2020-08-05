package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.drivers.SampleDriverData.DRIVER_TYPE_ALL;

public class DriverPrevData
{
    public final String CancerType;
    public final String Gene;
    public final String DriverType;
    public final double Prevalence;

    public DriverPrevData(final String cancerType, final String gene, final String driverType, final double prevalence)
    {
        CancerType = cancerType;
        Gene = gene;
        DriverType = driverType;
        Prevalence = prevalence;
    }

    public static DriverPrevData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 4)
            return null;

        return new DriverPrevData(items[0], items[1], items[2], Double.parseDouble(items[3]));
    }

    public boolean isTypeAll() { return DriverType.equals(DRIVER_TYPE_ALL); }

    public String toString() { return String.format("ct(%s) gene(%s) type(%s) prev(%.4f)", CancerType, Gene, DriverType, Prevalence); }

}
