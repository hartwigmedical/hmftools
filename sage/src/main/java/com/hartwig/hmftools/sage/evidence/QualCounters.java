package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.round;
import static java.lang.String.format;

import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.quality.QualityScores;

public class QualCounters
{
    private int mBaseQualityTotal;
    private int mAltBaseQualityTotal;
    private double mModifiedAltBaseQualityTotal;

    private int mMapQualityTotal;
    private int mAltMapQualityTotal;
    private double mModifiedAltMapQualityTotal;
    private int mLowQualAltSupportCount;

    public QualCounters()
    {
        mBaseQualityTotal = 0;
        mAltBaseQualityTotal = 0;
        mModifiedAltBaseQualityTotal = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
        mModifiedAltMapQualityTotal = 0;
        mLowQualAltSupportCount = 0;
    }

    public void update(final QualityScores qualityScores, int mapQuality, final ReadContextMatch matchType)
    {
        mBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        mMapQualityTotal += mapQuality;

        if(matchType.FullAltSupport)
        {
            mModifiedAltBaseQualityTotal += qualityScores.ModifiedBaseQuality;
            mModifiedAltMapQualityTotal += qualityScores.ModifiedMapQuality;
        }

        if(matchType.SupportsAlt)
        {
            mAltBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
            mAltMapQualityTotal += mapQuality;
        }
    }

    public void update(final QualityScores qualityScores)
    {
        mAltBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        mLowQualAltSupportCount += 1;
    }

    public int baseQualityTotal() { return mBaseQualityTotal; }
    public int altBaseQualityTotal() { return mAltBaseQualityTotal; }
    public double modifiedAltBaseQualityTotal() { return mModifiedAltBaseQualityTotal; }
    public int mapQualityTotal() { return mMapQualityTotal; }
    public int altMapQualityTotal() { return mAltMapQualityTotal; }
    public double altModifiedMapQualityTotal() { return mModifiedAltMapQualityTotal; }
    public int lowQualAltSupportCount() { return mLowQualAltSupportCount; }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                mBaseQualityTotal, mAltBaseQualityTotal, mMapQualityTotal, mAltMapQualityTotal);
    }
}
