package com.hartwig.hmftools.sig_analyser.cup;

import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.DATA_DELIM;
import static com.hartwig.hmftools.sig_analyser.cup.CupSampleDriverData.DRIVER_TYPE_ALL;

public class CupDriverPrevalence
{
    public final String CancerType;
    public final String Gene;
    public final String DriverType;
    public final double Prevalence;

    public CupDriverPrevalence(final String cancerType, final String gene, final String driverType, final double prevalence)
    {
        CancerType = cancerType;
        Gene = gene;
        DriverType = driverType;
        Prevalence = prevalence;
    }

    public static CupDriverPrevalence from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 4)
            return null;

        return new CupDriverPrevalence(items[0], items[1], items[2], Double.parseDouble(items[3]));
    }

    public boolean isTypeAll() { return DriverType.equals(DRIVER_TYPE_ALL); }

    public String toString() { return String.format("ct(%s) gene(%s) type(%s) prev(%.4f)", CancerType, Gene, DriverType, Prevalence); }

}
