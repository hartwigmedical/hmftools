package com.hartwig.hmftools.redux.jitter;

import java.util.List;

import org.jetbrains.annotations.Nullable;

public class MicrosatelliteSelector
{
    private final List<String> mUnits;
    private final Integer mUnitLengthMin;
    private final Integer mUnitLengthMax; // inclusive

    private MicrosatelliteSelector(@Nullable List<String> units, @Nullable Integer unitLengthMin, @Nullable Integer unitLengthMax)
    {
        mUnits = units;
        mUnitLengthMin = unitLengthMin;
        mUnitLengthMax = unitLengthMax;
    }

    public boolean select(final MicrosatelliteSiteData microsatelliteSiteData)
    {
        // select the type of microsatellite based on the repeat unit
        MicrosatelliteSite msiSite = microsatelliteSiteData.refGenomeMicrosatellite();

        if(mUnits != null)
        {
            return mUnits.contains(msiSite.unitString());
        }
        else if(mUnitLengthMin != null)
        {
            if(mUnitLengthMax == null)
            {
                return msiSite.Unit.length == mUnitLengthMin;
            }
            else
            {
                return msiSite.Unit.length >= mUnitLengthMin && msiSite.Unit.length <= mUnitLengthMax;
            }
        }
        return false;
    }

    public String unitName()
    {
        if(mUnits != null)
        {
            return String.join("/", mUnits);
        }
        else if(mUnitLengthMin != null)
        {
            if(mUnitLengthMax == null)
            {
                return mUnitLengthMin + "bp repeat";
            }
            else
            {
                return mUnitLengthMin + "-" + mUnitLengthMax + "bp repeat";
            }
        }

        return "";
    }

    public static MicrosatelliteSelector fromUnits(List<String> units)
    {
        return new MicrosatelliteSelector(units, null, null);
    }

    public static MicrosatelliteSelector fromUnitLength(int unitLength)
    {
        return new MicrosatelliteSelector(null, unitLength, null);
    }

    public static MicrosatelliteSelector fromUnitLengthRange(int unitLengthMin, int unitLengthMaxInclusive)
    {
        return new MicrosatelliteSelector(null, unitLengthMin, unitLengthMaxInclusive);
    }
}
