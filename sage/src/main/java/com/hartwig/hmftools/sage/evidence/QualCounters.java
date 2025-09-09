package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.quality.QualityCalculator.isMediumBaseQual;

import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.quality.QualityScores;

public class QualCounters
{
    private int mAltSeqTechBaseQualityTotal;

    private int mRecalibratedBaseQualityTotal;
    private int mAltRecalibratedBaseQualityTotal;
    private double mStrongAltRecalibratedBaseQualityTotal;
    private double mStrongAltMediumRecalibratedBaseQualityTotal;

    private double mAltFinalBaseQualityTotal;

    private int mMapQualityTotal;
    private int mAltMapQualityTotal;
    private double mAltFinalMapQualityTotal;

    private int mLowQualAltSupportCount;

    public QualCounters()
    {
        mRecalibratedBaseQualityTotal = 0;
        mAltRecalibratedBaseQualityTotal = 0;
        mAltSeqTechBaseQualityTotal = 0;
        mAltFinalBaseQualityTotal = 0;
        mStrongAltRecalibratedBaseQualityTotal = 0;
        mStrongAltMediumRecalibratedBaseQualityTotal = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
        mAltFinalMapQualityTotal = 0;
        mLowQualAltSupportCount = 0;
    }

    public void update(final QualityScores qualityScores, int mapQuality, final ReadContextMatch matchType)
    {
        mRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        mMapQualityTotal += mapQuality;

        if(matchType.FullAltSupport)
        {
            mAltFinalBaseQualityTotal += qualityScores.FinalBaseQuality;

            mStrongAltRecalibratedBaseQualityTotal += qualityScores.RecalibratedBaseQuality;

            if(isMediumBaseQual(qualityScores.SeqTechBaseQuality))
                mStrongAltMediumRecalibratedBaseQualityTotal += qualityScores.RecalibratedBaseQuality;

            mAltFinalMapQualityTotal += qualityScores.FinalMapQuality;
        }

        if(matchType.SupportsAlt)
        {
            mAltRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
            mAltSeqTechBaseQualityTotal += (int)round(qualityScores.SeqTechBaseQuality);
            mAltMapQualityTotal += mapQuality;
        }
    }

    public void update(final QualityScores qualityScores)
    {
        mAltSeqTechBaseQualityTotal += (int)round(qualityScores.SeqTechBaseQuality);
        mLowQualAltSupportCount += 1;
    }

    public int recalibratedBaseQualityTotal() { return mRecalibratedBaseQualityTotal; }
    public int altSeqTechBaseQualityTotal() { return mAltSeqTechBaseQualityTotal; }

    public int altRecalibratedBaseQualityTotal() { return mAltRecalibratedBaseQualityTotal; }
    public double strongAltMediumRecalibratedBaseQualityTotal() { return mStrongAltMediumRecalibratedBaseQualityTotal; }
    public double strongAltRecalibratedBaseQualityTotal() { return mStrongAltRecalibratedBaseQualityTotal; }

    public double altFinalBaseQualityTotal() { return mAltFinalBaseQualityTotal; }

    public int mapQualityTotal() { return mMapQualityTotal; }
    public int altMapQualityTotal() { return mAltMapQualityTotal; }
    public double altFinalMapQualityTotal() { return mAltFinalMapQualityTotal; }
    public int lowQualAltSupportCount() { return mLowQualAltSupportCount; }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                mRecalibratedBaseQualityTotal, mAltRecalibratedBaseQualityTotal, mMapQualityTotal, mAltMapQualityTotal);
    }
}
