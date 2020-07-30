package com.hartwig.hmftools.sig_analyser.cup;

import static com.hartwig.hmftools.sig_analyser.cup.CupConfig.DATA_DELIM;

public class CupSampleDriverData
{
    public final String SampleId;
    public final String Gene;
    public final String DriverType;

    public static final String DRIVER_TYPE_ALL = "ALL";

    public CupSampleDriverData(final String sampleId, final String gene, final String driverType)
    {
        SampleId = sampleId;
        Gene = gene;
        DriverType = driverType;
    }

    public static CupSampleDriverData from(final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        if(items.length != 3)
            return null;

        return new CupSampleDriverData(items[0], items[1], items[2]);
    }

    public boolean isTypeAll() { return DriverType.equals(DRIVER_TYPE_ALL); }
}
