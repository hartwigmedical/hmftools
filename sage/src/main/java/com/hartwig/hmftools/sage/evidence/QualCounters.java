package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.round;
import static java.lang.String.format;

public class QualCounters
{
    private int mBaseQualityTotal;
    private int mAltBaseQualityTotal;
    private int mMapQualityTotal;
    private int mAltMapQualityTotal;

    public QualCounters()
    {
        mBaseQualityTotal = 0;
        mAltBaseQualityTotal = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
    }

    public void update(double baseQuality, int mapQuality, boolean supportsAlt)
    {
        mBaseQualityTotal += (int)round(baseQuality);
        mMapQualityTotal += mapQuality;

        if(supportsAlt)
        {
            mAltBaseQualityTotal += (int)round(baseQuality);
            mAltMapQualityTotal += mapQuality;
        }
    }

    public int baseQualityTotal() { return mBaseQualityTotal; }
    public int altBaseQualityTotal() { return mAltBaseQualityTotal; }
    public int mapQualityTotal() { return mMapQualityTotal; }
    public int altMapQualityTotal() { return mAltMapQualityTotal; }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                mBaseQualityTotal, mAltBaseQualityTotal, mMapQualityTotal, mAltMapQualityTotal);
    }
}
