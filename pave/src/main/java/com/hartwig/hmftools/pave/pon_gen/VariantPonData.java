package com.hartwig.hmftools.pave.pon_gen;

import static java.lang.Math.max;
import static java.lang.String.format;

import java.util.Comparator;

import com.hartwig.hmftools.common.pathogenic.Pathogenicity;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.pon.MultiPonStatus;

public class VariantPonData extends SimpleVariant
{
    // counts collected across all samples
    private int mSampleCount;
    private int mTotalReadCount;
    private int mMaxSampleReadCount;

    private int mRepeatCount;
    private Pathogenicity mClinvarPathogenicity;
    private boolean mSomaticHotspot;
    private boolean mGermlineHotspot;
    private boolean mInCodingRegion;
    private boolean mInBasePonCache;
    private boolean mInMultipleCaches;

    public VariantPonData(final String chromosome, final int position, final String ref, final String alt)
    {
        super(chromosome, position, ref, alt);

        mRepeatCount = 0;
        mSampleCount = 0;
        mTotalReadCount = 0;
        mMaxSampleReadCount = 0;
        mClinvarPathogenicity = null;
        mInCodingRegion = false;
        mSomaticHotspot = false;
        mGermlineHotspot = false;
        mInBasePonCache = false;
        mInMultipleCaches = false;
    }

    public int sampleCount() { return mSampleCount; }
    public int totalReadCount() { return mTotalReadCount; }
    public int maxSampleReadCount() { return mMaxSampleReadCount; }

    public void addSampleData(int readSupport)
    {
        ++mSampleCount;

        mMaxSampleReadCount = max(mMaxSampleReadCount, readSupport);
        mTotalReadCount += readSupport;
    }

    public void setSampleCount(int sampleCount) { mSampleCount = sampleCount; }
    public void setMaxSampleReadCount(int readCount) { mMaxSampleReadCount = readCount; }
    public void setTotalReadCount(int readCount) { mTotalReadCount = readCount; }

    public void setClinvarPathogenicity(final Pathogenicity pathogenicity) { mClinvarPathogenicity = pathogenicity; }
    public Pathogenicity clinvarPathogenicity() { return mClinvarPathogenicity; }
    public boolean isClinvarPathogenic() { return mClinvarPathogenicity != null && mClinvarPathogenicity.isPathogenic(); }

    public void setRepeatCount(int repeatCount) { mRepeatCount = repeatCount; }
    public int repeatCount() { return mRepeatCount; }

    public void markInBasePonCache() { mInBasePonCache = true; }
    public boolean inBasePonCache() { return mInBasePonCache; }

    public void markInMultiplePonCaches() { mInMultipleCaches = true; }
    public boolean inMultipleCaches() { return mInMultipleCaches; }

    public MultiPonStatus multiPonStatus()
    {
        if(mInBasePonCache && mInMultipleCaches)
            return MultiPonStatus.MULTI;

        if(mInBasePonCache)
            return MultiPonStatus.BASE;

        return MultiPonStatus.ARTEFACT;
    }

    public void markSomaticHotspot() { mSomaticHotspot = true; }
    public boolean isSomaticHotspot() { return mSomaticHotspot; }

    public void markGermlineHotspot() { mGermlineHotspot = true; }
    public boolean isGermlineHotspot() { return mGermlineHotspot; }

    public void markInCodingRegion() { mInCodingRegion = true; }
    public boolean inCodingRegion() { return mInCodingRegion; }

    public String toString()
    {
        return format("var(%s) samples(%d)", super.toString(), mSampleCount);
    }

    public static class VariantSorter implements Comparator<VariantPonData>
    {
        // note does not include chromosome since assumes just sorting within the same region
        public int compare(final VariantPonData first, final VariantPonData second)
        {
            // sort on position and then ref and alt
            if(first.Position != second.Position)
                return first.Position < second.Position ? -1 : 1;

            if(!first.Ref.equals(second.Ref))
                return first.Ref.compareTo(second.Ref);

            return first.Alt.compareTo(second.Alt);
        }
    }
}
