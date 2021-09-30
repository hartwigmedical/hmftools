package com.hartwig.hmftools.purple.gene;

import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class GeneCopyNumberBuilder implements RegionZipperHandler<PurpleCopyNumber,HmfExonRegion>
{
    private final ImmutableGeneCopyNumber.Builder mBuilder;

    private double mMinCopyNumber = Double.MAX_VALUE;
    private double mMinMinorAllelePloidy = Double.MAX_VALUE;
    private double mMaxCopyNumber = -Double.MAX_VALUE;
    private int mSomaticCount;
    private int mHomCount;
    private int mHet2HomCount;

    private PurpleCopyNumber mPrevious;

    private double mPreviousCopyNumber = -Double.MAX_VALUE;
    private HmfExonRegion mExon;
    private PurpleCopyNumber mCopyNumber;

    private int mMinRegions = 0;
    private long mMinRegionStart = 0;
    private long mMinRegionEnd = 0;
    private SegmentSupport mMinRegionStartSupport = SegmentSupport.NONE;
    private SegmentSupport mMinRegionEndSupport = SegmentSupport.NONE;
    private CopyNumberMethod mMinRegionMethod = CopyNumberMethod.UNKNOWN;

    public GeneCopyNumberBuilder(final HmfTranscriptRegion gene)
    {
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
        this.mCopyNumber = copyNumber;
        if(mExon != null)
        {
            addOverlap(mExon, this.mCopyNumber);
        }
    }

    @Override
    public void secondary(final HmfExonRegion exon)
    {
        this.mExon = exon;
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
                        mHomCount++;
                        break;
                    case GERMLINE_HET2HOM_DELETION:
                        mHet2HomCount++;
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
                .germlineHomRegions(mHomCount)
                .germlineHet2HomRegions(mHet2HomCount)
                .minRegions(mMinRegions)
                .minMinorAlleleCopyNumber(mMinMinorAllelePloidy)
                .build();
    }
}
