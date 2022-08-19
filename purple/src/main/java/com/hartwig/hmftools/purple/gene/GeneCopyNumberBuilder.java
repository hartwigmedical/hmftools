package com.hartwig.hmftools.purple.gene;

import static java.lang.Math.min;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

public class GeneCopyNumberBuilder
{
    private final TranscriptData mTransData;
    private final List<PurpleCopyNumber> mCopyNumbers;
    private final ImmutableGeneCopyNumber.Builder mBuilder;

    private double mMinCopyNumber;
    private double mMinMinorAllelePloidy;
    private double mMaxCopyNumber;
    private int mSomaticCount;

    private PurpleCopyNumber mPrevious;

    private double mPreviousCopyNumber;
    private ExonData mExon;
    private PurpleCopyNumber mCopyNumber;

    private int mMinRegions;
    private long mMinRegionStart;
    private long mMinRegionEnd;
    private SegmentSupport mMinRegionStartSupport;
    private SegmentSupport mMinRegionEndSupport;
    private CopyNumberMethod mMinRegionMethod;

    public static List<GeneCopyNumber> createGeneCopyNumbers(
            final RefGenomeVersion refGenomeVersion, final EnsemblDataCache geneTransCache, final List<PurpleCopyNumber> copyNumbers)
    {
        final List<GeneCopyNumber> result = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrString = refGenomeVersion.versionedChromosome(chromosome.toString());

            List<GeneData> geneDataList = geneTransCache.getChrGeneDataMap().get(chrString);

            List<PurpleCopyNumber> chromosomeCopyNumbers = copyNumbers.stream()
                    .filter(x -> x.chromosome().equals(chrString)).collect(Collectors.toList());

            for(GeneData geneData : geneDataList)
            {
                List<TranscriptData> transDataList = geneTransCache.getTranscripts(geneData.GeneId);

                for(TranscriptData tranData : transDataList)
                {
                    final GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(geneData, tranData, chromosomeCopyNumbers);
                    GeneCopyNumber geneCopyNumber2 = builder.create();

                    if(geneCopyNumber2.totalRegions() > 0)
                        result.add(geneCopyNumber2);
                }
            }
        }

        return result;
    }

    public GeneCopyNumberBuilder(final GeneData geneData, final TranscriptData transData, final List<PurpleCopyNumber> copyNumbers)
    {
        mCopyNumbers = copyNumbers;
        mTransData = transData;

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
                .geneName(geneData.GeneName)
                .transName(transData.TransName)
                .isCanonical(transData.IsCanonical)
                .chromosomeBand(geneData.KaryotypeBand)
                .chromosome(geneData.Chromosome)
                .start(geneData.GeneStart)
                .end(geneData.GeneEnd)
                .minRegionStart(geneData.GeneStart)
                .minRegionEnd(geneData.GeneEnd)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEndSupport(SegmentSupport.NONE);
    }

    public GeneCopyNumber create()
    {
        int cnIndex = 0;
        int exonIndex = 0;

        while(cnIndex < mCopyNumbers.size() || exonIndex < mTransData.exons().size())
        {
            PurpleCopyNumber copyNumber = cnIndex < mCopyNumbers.size() ? mCopyNumbers.get(cnIndex) : null;
            ExonData exonData = exonIndex < mTransData.exons().size() ? mTransData.exons().get(exonIndex) : null;

            if(copyNumber == null || (exonData != null && exonData.Start < copyNumber.start()))
            {
                handleExon(exonData);
                ++exonIndex;
            }
            else
            {
                handleCopyNumber(copyNumber);
                ++cnIndex;
            }
        }

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

    private void handleCopyNumber(final PurpleCopyNumber copyNumber)
    {
        mCopyNumber = copyNumber;
        if(mExon != null)
        {
            addOverlap(mExon, this.mCopyNumber);
        }
    }

    private void handleExon(final ExonData exon)
    {
        mExon = exon;
        if(mCopyNumber != null)
        {
            addOverlap(mExon, mCopyNumber);
        }
    }

    private void addOverlap(final ExonData exon, final PurpleCopyNumber copyNumber)
    {
        int minEnd = min(exon.End, copyNumber.end());
        int maxStart = Math.max(exon.Start, copyNumber.start());
        int overlap = Math.max(0, 1 + minEnd - maxStart);

        if(overlap > 0)
        {
            double currentCopyNumber = copyNumber.averageTumorCopyNumber();

            mMaxCopyNumber = Math.max(mMaxCopyNumber, currentCopyNumber);
            mMinMinorAllelePloidy = min(mMinMinorAllelePloidy, copyNumber.minorAlleleCopyNumber());

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
}
