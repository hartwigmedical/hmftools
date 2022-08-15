package com.hartwig.hmftools.purple.gene;

import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class GeneCopyNumberBuilder implements RegionZipperHandler<PurpleCopyNumber,HmfExonRegion>
{
    private final ImmutableGeneCopyNumber.Builder mBuilder;

    private double mMinCopyNumber;
    private double mMinMinorAllelePloidy;
    private double mMaxCopyNumber;
    private int mSomaticCount;

    private PurpleCopyNumber mPrevious;

    private double mPreviousCopyNumber;
    private HmfExonRegion mExon;
    private PurpleCopyNumber mCopyNumber;

    private int mMinRegions;
    private long mMinRegionStart;
    private long mMinRegionEnd;
    private SegmentSupport mMinRegionStartSupport;
    private SegmentSupport mMinRegionEndSupport;
    private CopyNumberMethod mMinRegionMethod;

    public GeneCopyNumberBuilder(final HmfTranscriptRegion gene)
    {
        mMinCopyNumber = Double.MAX_VALUE;
        mMinMinorAllelePloidy = Double.MAX_VALUE;
        mMaxCopyNumber = -Double.MAX_VALUE;
        mSomaticCount = 0;

        mPrevious = null;

        mPreviousCopyNumber = -Double.MAX_VALUE;
        mExon = null;
        mCopyNumber = null;

        mMinRegions = 0;
        mMinRegionStart = 0;
        mMinRegionEnd = 0;
        mMinRegionStartSupport = SegmentSupport.NONE;
        mMinRegionEndSupport = SegmentSupport.NONE;
        mMinRegionMethod = CopyNumberMethod.UNKNOWN;

        mBuilder = ImmutableGeneCopyNumber.builder()
                .from(gene)
                .minRegionStart(gene.start())
                .minRegionEnd(gene.end())
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE);
    }

    @Override
    public void enterChromosome(final String chromosome)
    {
        // IGNORE
    }

    @Override
    public void primary(final PurpleCopyNumber copyNumber)
    {
        mCopyNumber = copyNumber;
        if(mExon != null)
        {
            addOverlap(mExon, this.mCopyNumber);
        }
    }

    @Override
    public void secondary(final HmfExonRegion exon)
    {
        mExon = exon;
        if(mCopyNumber != null)
        {
            addOverlap(this.mExon, mCopyNumber);
        }
    }

    private void addOverlap(final HmfExonRegion exon, final PurpleCopyNumber copyNumber)
    {
        long overlap = exon.overlappingBases(copyNumber);
        if(overlap > 0)
        {
            double currentCopyNumber = copyNumber.averageTumorCopyNumber();

            mMaxCopyNumber = Math.max(mMaxCopyNumber, currentCopyNumber);
            mMinMinorAllelePloidy = Math.min(mMinMinorAllelePloidy, copyNumber.minorAlleleCopyNumber());

            if(!Doubles.equal(currentCopyNumber, mPreviousCopyNumber))
            {
                switch(copyNumber.method())
                {
                    case GERMLINE_HOM_DELETION:
                    case GERMLINE_HET2HOM_DELETION:
                        break;
                    default:
                        mSomaticCount++;
                }
            }

            if(isUnprocessedCopyNumberRegion(copyNumber))
            {
                if(Doubles.lessThan(currentCopyNumber, mMinCopyNumber))
                {
                    mMinRegions = 1;
                    mMinCopyNumber = currentCopyNumber;
                    mMinRegionStart = copyNumber.start();
                    mMinRegionStartSupport = copyNumber.segmentStartSupport();
                    mMinRegionEnd = copyNumber.end();
                    mMinRegionEndSupport = copyNumber.segmentEndSupport();
                    mMinRegionMethod = copyNumber.method();

                }
                else if(Doubles.equal(currentCopyNumber, mMinCopyNumber))
                {
                    mMinRegionEnd = copyNumber.end();
                    mMinRegionEndSupport = copyNumber.segmentEndSupport();
                    mMinRegionMethod = copyNumber.method();

                    if(!Doubles.equal(currentCopyNumber, mPreviousCopyNumber))
                    {
                        mMinRegions++;
                    }
                }
            }

            mPreviousCopyNumber = currentCopyNumber;
            mPrevious = copyNumber;
        }
    }

    private boolean isUnprocessedCopyNumberRegion(final PurpleCopyNumber copyNumber)
    {
        return mPrevious == null || !mPrevious.equals(copyNumber);
    }

    @NotNull
    public GeneCopyNumber build()
    {
        return mBuilder.maxCopyNumber(mMaxCopyNumber)
                .minRegionStartSupport(mMinRegionStartSupport)
                .minRegionEndSupport(mMinRegionEndSupport)
                .minRegionMethod(mMinRegionMethod)
                .minRegionStart(mMinRegionStart)
                .minRegionEnd(mMinRegionEnd)
                .minCopyNumber(mMinCopyNumber)
                .somaticRegions(mSomaticCount)
                .minRegions(mMinRegions)
                .minMinorAlleleCopyNumber(mMinMinorAllelePloidy)
                .build();
    }
}
