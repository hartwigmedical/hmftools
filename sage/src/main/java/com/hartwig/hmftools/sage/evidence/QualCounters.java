package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageConstants.MAX_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.bqr.BqrRegionReader.extractReadType;
import static java.lang.Math.round;
import static java.lang.String.format;

import com.hartwig.hmftools.common.qual.BqrReadType;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.quality.QualityScores;
import htsjdk.samtools.SAMRecord;

public class QualCounters
{
    private int mRecalibratedBaseQualityTotal;
    private int mAltRecalibratedBaseQualityTotal;
    private int mStrongSimplexAltRecalibratedBaseQualityTotal;
    private int mStrongDuplexAltRecalibratedBaseQualityTotal;
    private int mAltBaseQualityTotal;
    private double mModifiedAltBaseQualityTotal;
    private int mMapQualityTotal;
    private int mAltMapQualityTotal;
    private double mModifiedAltMapQualityTotal;
    private int mLowQualAltSupportCount;

    public QualCounters()
    {
        mRecalibratedBaseQualityTotal = 0;
        mAltRecalibratedBaseQualityTotal = 0;
        mStrongSimplexAltRecalibratedBaseQualityTotal = 0;
        mStrongDuplexAltRecalibratedBaseQualityTotal = 0;
        mAltBaseQualityTotal = 0;
        mModifiedAltBaseQualityTotal = 0;
        mMapQualityTotal = 0;
        mAltMapQualityTotal = 0;
        mModifiedAltMapQualityTotal = 0;
        mLowQualAltSupportCount = 0;
    }

    public void update(final SAMRecord record, final QualityScores qualityScores, final ReadContextMatch matchType, final int readIndex)
    {
        int mapQuality = record.getMappingQuality();
        mRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        mMapQualityTotal += mapQuality;

        if(matchType.FullAltSupport)
        {
            mModifiedAltBaseQualityTotal += qualityScores.ModifiedBaseQuality;
            mModifiedAltMapQualityTotal += qualityScores.ModifiedMapQuality;
            BqrReadType readType = extractReadType(record, SequencingType.SBX, record.getBaseQualities()[readIndex]);
            if(readType.equals(BqrReadType.DUAL))
                mStrongDuplexAltRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
            else
                mStrongSimplexAltRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
        }

        if(matchType.SupportsAlt)
        {
            mAltRecalibratedBaseQualityTotal += (int)round(qualityScores.RecalibratedBaseQuality);
            mAltBaseQualityTotal += Math.min((int)round(qualityScores.CalcBaseQuality), MAX_RAW_BASE_QUAL);
            mAltMapQualityTotal += mapQuality;
        }
    }

    public void update(final QualityScores qualityScores)
    {
        mAltBaseQualityTotal += Math.min((int)round(qualityScores.CalcBaseQuality), MAX_RAW_BASE_QUAL);
        mLowQualAltSupportCount += 1;
    }

    public int baseQualityTotal() { return mRecalibratedBaseQualityTotal; }
    public int altRecalibratedBaseQualityTotal() { return mAltRecalibratedBaseQualityTotal; }
    public int strongSimplexAltRecalibratedBaseQualityTotal() { return mStrongSimplexAltRecalibratedBaseQualityTotal; }
    public int strongDuplexAltRecalibratedBaseQualityTotal() { return mStrongDuplexAltRecalibratedBaseQualityTotal; }

    public int altBaseQualityTotal() { return mAltBaseQualityTotal; }
    public double modifiedAltBaseQualityTotal() { return mModifiedAltBaseQualityTotal; }
    public int mapQualityTotal() { return mMapQualityTotal; }
    public int altMapQualityTotal() { return mAltMapQualityTotal; }
    public double altModifiedMapQualityTotal() { return mModifiedAltMapQualityTotal; }
    public int lowQualAltSupportCount() { return mLowQualAltSupportCount; }

    public String toString()
    {
        return format("baseQualTotal(%d alt=%d) mapQualTotal(%d alt=%d)",
                mRecalibratedBaseQualityTotal, mAltRecalibratedBaseQualityTotal, mMapQualityTotal, mAltMapQualityTotal);
    }
}
