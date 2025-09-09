package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.isMediumBaseQual;

import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.quality.QualityScores;

public class QualCounters
{
    private int mAltBaseQualityTotal;

    private int mRecalibratedBaseQualityTotal;
    private int mRecalibratedAltBaseQualityTotal;

    private double mRecalibratedStrongAltBaseQualityTotal;
    private double mRecalibratedStrongAltMediumBaseQualityTotal;

    private double mModifiedAltBaseQualityTotal;

    private int mMapQualityTotal;
    private int mAltMapQualityTotal;
    private double mModifiedAltMapQualityTotal;

    private int mLowQualAltSupportCount;

    public QualCounters()
    {
        mRecalibratedBaseQualityTotal = 0;
        mRecalibratedAltBaseQualityTotal = 0;
        mAltBaseQualityTotal = 0;
        mModifiedAltBaseQualityTotal = 0;
        mRecalibratedStrongAltBaseQualityTotal = 0;
        mRecalibratedStrongAltMediumBaseQualityTotal = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
        mModifiedAltMapQualityTotal = 0;
        mLowQualAltSupportCount = 0;
    }

    public void update(final QualityScores qualityScores, int mapQuality, final ReadContextMatch matchType)
    {
        mRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        mMapQualityTotal += mapQuality;

        if(matchType.FullAltSupport)
        {
            mModifiedAltBaseQualityTotal += qualityScores.ModifiedBaseQuality;

            mRecalibratedStrongAltBaseQualityTotal += qualityScores.RecalibratedBaseQuality;

            if(isMediumBaseQual(qualityScores.CalcBaseQuality))
                mRecalibratedStrongAltMediumBaseQualityTotal += qualityScores.RecalibratedBaseQuality;

            mModifiedAltMapQualityTotal += qualityScores.ModifiedMapQuality;
        }

        if(matchType.SupportsAlt)
        {
            mRecalibratedAltBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
            mAltBaseQualityTotal += (int)round(qualityScores.CalcBaseQuality);
            mAltMapQualityTotal += mapQuality;
        }
    }

    public void update(final QualityScores qualityScores)
    {
        mAltBaseQualityTotal += (int)round(qualityScores.CalcBaseQuality);
        mLowQualAltSupportCount += 1;
    }

    public int baseQualityTotal() { return mRecalibratedBaseQualityTotal; }
    public int altBaseQualityTotal() { return mAltBaseQualityTotal; }

    public int recalibratedAltBaseQualityTotal() { return mRecalibratedAltBaseQualityTotal; }
    public double recalibratedStrongAltMediumBaseQualityTotal() { return mRecalibratedStrongAltMediumBaseQualityTotal; }
    public double recalibratedStrongAltBaseQualityTotal() { return mRecalibratedStrongAltBaseQualityTotal; }

    public double modifiedAltBaseQualityTotal() { return mModifiedAltBaseQualityTotal; }

    public int mapQualityTotal() { return mMapQualityTotal; }
    public int altMapQualityTotal() { return mAltMapQualityTotal; }
    public double altModifiedMapQualityTotal() { return mModifiedAltMapQualityTotal; }
    public int lowQualAltSupportCount() { return mLowQualAltSupportCount; }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                mRecalibratedBaseQualityTotal, mRecalibratedAltBaseQualityTotal, mMapQualityTotal, mAltMapQualityTotal);
    }
}
