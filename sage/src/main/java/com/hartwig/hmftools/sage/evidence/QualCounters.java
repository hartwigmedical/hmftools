package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.round;
import static java.lang.String.format;

import com.hartwig.hmftools.sage.quality.QualityScores;

public class QualCounters
{
    private int mBaseQualityTotal;
    private int mAltBaseQualityTotal;
    private double mModifiedBaseQualityTotal;

    private int mMapQualityTotal;
    private int mAltMapQualityTotal;
    private double mModifiedMapQualityTotal;

    public QualCounters()
    {
        mBaseQualityTotal = 0;
        mAltBaseQualityTotal = 0;
        mModifiedBaseQualityTotal = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
        mModifiedMapQualityTotal = 0;
    }

    public void update(final QualityScores qualityScores, int mapQuality, boolean supportsAlt)
    {
        mBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        mMapQualityTotal += mapQuality;
        mModifiedBaseQualityTotal += qualityScores.ModifiedBaseQuality;
        mModifiedMapQualityTotal += qualityScores.ModifiedMapQuality;

        if(supportsAlt)
        {
            mAltBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
            mAltMapQualityTotal += mapQuality;
        }
    }

    public int baseQualityTotal() { return mBaseQualityTotal; }
    public int altBaseQualityTotal() { return mAltBaseQualityTotal; }
    public double modifiedBaseQualityTotal() { return mModifiedBaseQualityTotal; }
    public int mapQualityTotal() { return mMapQualityTotal; }
    public int altMapQualityTotal() { return mAltMapQualityTotal; }
    public double modifiedMapQualityTotal() { return mModifiedMapQualityTotal; }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                mBaseQualityTotal, mAltBaseQualityTotal, mMapQualityTotal, mAltMapQualityTotal);
    }
}
