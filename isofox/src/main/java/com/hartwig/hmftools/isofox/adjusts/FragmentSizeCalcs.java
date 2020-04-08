package com.hartwig.hmftools.isofox.adjusts;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_FREQUENCY;
import static com.hartwig.hmftools.isofox.exp_rates.ExpectedRatesGenerator.FL_LENGTH;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.BamFragmentAllocator;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// calculate fragment length distribution for a sample
// this can be done either independently from fragment counting or lengths can be registering during that process

public class FragmentSizeCalcs
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final List<int[]> mFragmentLengths;
    private final List<int[]> mFragmentLengthsByGene;
    private int mMaxReadLength;

    private BufferedWriter mFragLengthWriter;

    private EnsemblGeneData mCurrentGeneData;
    private List<TranscriptData> mCurrentTransDataList;
    private int mCurrentReadCount;
    private int mTotalReadCount;
    private int mProcessedFragments;

    private PerformanceCounter mPerfCounter;

    public FragmentSizeCalcs(final IsofoxConfig config, final EnsemblDataCache geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        mBamSlicer = new BamSlicer(BamFragmentAllocator.DEFAULT_MIN_MAPPING_QUALITY, true);

        mCurrentGeneData = null;
        mCurrentTransDataList = null;
        mCurrentReadCount = 0;
        mTotalReadCount = 0;
        mProcessedFragments = 0;
        mMaxReadLength = 0;

        mFragLengthWriter = writer;

        mFragmentLengths = Lists.newArrayList();
        mFragmentLengthsByGene = Lists.newArrayList();
        mPerfCounter = new PerformanceCounter("FragLengthDist");
    }

    public final List<int[]> getFragmentLengths() { return mFragmentLengths; }
    public final int getMaxReadLength() { return mMaxReadLength; }

    private static final int MIN_GENE_LENGTH = 1000;
    private static final int MAX_GENE_LENGTH = 1000000;
    private static final int MAX_GENE_TRANS = 20;
    private static final int MAX_TRAN_EXONS = 20;
    private static final int MAX_GENE_READ_COUNT = 1000; // to avoid impact of highly enriched genes

    private static int FRAG_LENGTH_CAP = 3000; // to prevent map blowing out in size

    public void calcSampleFragmentSize(final String chromosome, final List<EnsemblGeneData> geneDataList, int requiredFragCount)
    {
        if (geneDataList == null || geneDataList.isEmpty())
            return;

        if (requiredFragCount == 0)
            return;

        // walk through each chromosome, ignoring any gene which overlaps the previous gene
        ISF_LOGGER.info("calculating fragment size for chromosome({}) geneCount({})", chromosome, geneDataList.size());

        List<long[]> excludedRegions = generateExcludedRegions(chromosome);
        long lastGeneEnd = 0;

        for (int i = 0; i < geneDataList.size(); ++i)
        {
            EnsemblGeneData geneData = geneDataList.get(i);

            if(mConfig.ExcludedGeneIds.contains(geneData.GeneId) || mConfig.EnrichedGeneIds.contains(geneData.GeneId))
                continue;

            if (geneData.GeneStart < lastGeneEnd)
                continue;

            long geneLength = geneData.length();

            if (geneLength < MIN_GENE_LENGTH || geneLength > MAX_GENE_LENGTH)
                continue;

            if(excludedRegions.stream().anyMatch(x -> positionsOverlap(x[SE_START], x[SE_END], geneData.GeneStart, geneData.GeneEnd)))
                continue;

            mCurrentTransDataList = mGeneTransCache.getTranscripts(geneData.GeneId).stream()
                    .filter(x -> x.exons().size() <= MAX_TRAN_EXONS)
                    .collect(Collectors.toList());

            if (mCurrentTransDataList.isEmpty() || mCurrentTransDataList.size() > MAX_GENE_TRANS)
                continue;

            if(i > 0 && (i % 100) == 0)
            {
                ISF_LOGGER.debug("chromosome({}) processed {} genes, lastGenePos({}) fragCount({}) totalReads({})",
                        chromosome, i, lastGeneEnd, mProcessedFragments, mTotalReadCount);
            }

            lastGeneEnd = geneData.GeneEnd;

            mPerfCounter.start();

            mCurrentReadCount = 0;
            mCurrentGeneData = geneData;
            List<GenomeRegion> regions = Lists.newArrayList(GenomeRegions.create(chromosome, geneData.GeneStart, geneData.GeneEnd));

            mBamSlicer.slice(mSamReader, regions, this::processBamRead);

            mPerfCounter.stop();

            if(mConfig.FragmentLengthsByGene)
            {
                writeFragmentLengths(mFragLengthWriter, mFragmentLengthsByGene, geneData);
                mFragmentLengthsByGene.clear();
            }

            if (mProcessedFragments >= requiredFragCount)
            {
                ISF_LOGGER.debug("max fragment length samples reached: {}", mProcessedFragments);
                break;
            }
        }

        ISF_LOGGER.debug("chromosome({}) processing complete", chromosome);
    }

    private List<long[]> generateExcludedRegions(final String chromosome)
    {
        // create a buffer around the enriched gene to avoid excessive reads in this vacinity
        int buffer = 100000;

        final List<long[]> excludedRegions = Lists.newArrayList();
        for(final String geneId : mConfig.EnrichedGeneIds)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataById(geneId);
            if(geneData == null)
                continue;

            if(geneData.Chromosome.equals(chromosome))
            {
                excludedRegions.add(new long[] { geneData.GeneStart - buffer, geneData.GeneEnd + buffer});
            }
        }

        return excludedRegions;
    }

    private void processBamRead(@NotNull final SAMRecord record)
    {
        ++mCurrentReadCount;
        ++mTotalReadCount;

        // ISF_LOGGER.trace("currentGene({}:{}) read count({})", mCurrentGeneData.GeneId, mCurrentGeneData.GeneName, mCurrentReadCount);

        if(mTotalReadCount > 0 && (mTotalReadCount % 10000) == 0)
        {
            ISF_LOGGER.trace("currentGene({}:{}) totalReads({})", mCurrentGeneData.GeneId, mCurrentGeneData.GeneName, mTotalReadCount);
        }

        if(mCurrentReadCount >= MAX_GENE_READ_COUNT)
        {
            ISF_LOGGER.trace("currentGene({}:{}) reached max read count", mCurrentGeneData.GeneId, mCurrentGeneData.GeneName);
            mBamSlicer.haltProcessing();
            return;
        }

        if(!isCandidateRecord(record))
            return;

        mMaxReadLength = max(mMaxReadLength, record.getReadLength());

        long posStart = record.getStart();
        long posEnd = record.getEnd();

        for(final TranscriptData transData : mCurrentTransDataList)
        {
            if(transData.exons().stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.ExonStart, x.ExonEnd)))
                return;
        }

        addFragmentLength(record, mFragmentLengths);

        if(mConfig.FragmentLengthsByGene)
            addFragmentLength(record, mFragmentLengthsByGene);
    }

    private boolean isCandidateRecord(final SAMRecord record)
    {
        if(!record.getFirstOfPairFlag())
            return false;

        // ignore translocations and inversions
        if(!record.getMateReferenceName().equals(record.getReferenceName()) || record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag())
            return false;

        // ignore split or unmapped reads
        if(record.getCigar() == null || record.getCigar().containsOperator(CigarOperator.N) || !record.getCigar().containsOperator(CigarOperator.M))
            return false;

        return true;
    }

    private synchronized void addFragmentLength(final SAMRecord record, final List<int[]> fragmentLengths)
    {
        int fragmentLength = min(abs(record.getInferredInsertSize()), FRAG_LENGTH_CAP);

        if(fragmentLength == 0)
            return;

        int index = 0;
        boolean exists = false;
        while(index < fragmentLengths.size())
        {
            final int[] fragLengthCount = fragmentLengths.get(index);

            if(fragLengthCount[FL_LENGTH] < fragmentLength)
            {
                ++index;
                continue;
            }

            if(fragLengthCount[FL_LENGTH] == fragmentLength)
            {
                ++fragLengthCount[FL_FREQUENCY];
                exists = true;
            }

            break;
        }

        if(!exists)
        {
            int[] newFragLengthCount = { fragmentLength, 1 };
            fragmentLengths.add(index, newFragLengthCount);
        }

        ++mProcessedFragments;
    }

    public static void setConfigFragmentLengthData(final IsofoxConfig config, int maxReadLength, final List<int[]> fragmentLengths)
    {
        if(maxReadLength > 0)
        {
            config.ReadLength = maxReadLength;
            ISF_LOGGER.info("max read length({}) set", maxReadLength);
        }
        else
        {
            ISF_LOGGER.warn("max read length not determined from fragment length calcs");
        }

        final List<int[]> lengthFrequencies = config.FragmentLengthData;

        int currentRangeMin = 0;
        int currentRangeMax = 0;

        for(int i = 0; i < lengthFrequencies.size(); ++i)
        {
            int[] lengthFrequency = lengthFrequencies.get(i);

            currentRangeMin = (i == 0) ? 0 : currentRangeMax + 1;

            if(i == lengthFrequencies.size() - 1)
            {
                currentRangeMax = config.MaxFragmentLength - 1;
            }
            else
            {
                int[] nextLengthFrequency = lengthFrequencies.get(i + 1);
                currentRangeMax = (lengthFrequency[FL_LENGTH] + nextLengthFrequency[FL_LENGTH]) / 2;
            }

            int lengthCount = 0;

            for (final int[] fragLengthCount : fragmentLengths)
            {
                if(fragLengthCount[FL_LENGTH] >= currentRangeMin && fragLengthCount[FL_LENGTH] <= currentRangeMax)
                {
                    lengthCount += fragLengthCount[FL_FREQUENCY];
                }
            }

            lengthFrequency[FL_FREQUENCY] = lengthCount;
            ISF_LOGGER.info("fragmentLength({}) frequency({})", lengthFrequency[FL_LENGTH], lengthCount);
        }
    }

    public static List<Double> calcPercentileData(final List<int[]> fragmentLengths, final List<Double> percentiles)
    {
        final List<Double> percentileLengths = Lists.newArrayList();

        double totalFragments = fragmentLengths.stream().mapToLong(x -> x[FL_FREQUENCY]).sum();

        for(Double percentile : percentiles)
        {
            long currentTotal = 0;
            int prevLength = 0;

            for (final int[] fragLengthCount : fragmentLengths)
            {
                double nextPercTotal = (currentTotal + fragLengthCount[FL_FREQUENCY]) / totalFragments;

                if(nextPercTotal >= percentile)
                {
                    double percLength = prevLength > 0 ? (prevLength + fragLengthCount[FL_LENGTH]) * 0.5 : fragLengthCount[FL_LENGTH];
                    percentileLengths.add(percLength);
                    break;
                }

                currentTotal += fragLengthCount[FL_FREQUENCY];
                prevLength = fragLengthCount[FL_LENGTH];
            }
        }

        return percentileLengths;
    }

    public static void mergeData(final List<int[]> fragmentLengths, final FragmentSizeCalcs other)
    {
        for(final int[] lengthData : other.getFragmentLengths())
        {
            int fragmentLength = lengthData[FL_LENGTH];
            int frequency = lengthData[FL_FREQUENCY];

            int index = 0;
            boolean exists = false;
            while (index < fragmentLengths.size())
            {
                final int[] fragLengthCount = fragmentLengths.get(index);

                if (fragLengthCount[FL_LENGTH] < fragmentLength)
                {
                    ++index;
                    continue;
                }

                if (fragLengthCount[FL_LENGTH] == fragmentLength)
                {
                    fragLengthCount[FL_FREQUENCY] += frequency;
                    exists = true;
                }

                break;
            }

            if (!exists)
            {
                int[] newFragLengthCount = { fragmentLength, frequency };
                fragmentLengths.add(index, newFragLengthCount);
            }
        }
    }

    public static BufferedWriter createFragmentLengthWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.FragmentLengthsByGene ?
                    config.formOutputFile("frag_length_by_gene.csv") : config.formOutputFile("frag_length.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            if (config.FragmentLengthsByGene)
            {
                writer.write("GeneId,GeneName,Chromosome,");
            }

            writer.write("FragmentLength,Count");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fragment length file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeFragmentLengths(
            final BufferedWriter writer, final List<int[]> fragmentLengths, final EnsemblGeneData geneData)
    {
        if(writer == null)
            return;

        if(fragmentLengths.isEmpty())
            return;

        try
        {
            for (final int[] fragLengthCount : fragmentLengths)
            {
                if(geneData != null)
                {
                    writer.write(String.format("%s,%s,%s,",
                            geneData.GeneId, geneData.GeneName, geneData.Chromosome));
                }

                writer.write(String.format("%d,%d", fragLengthCount[FL_LENGTH], fragLengthCount[FL_FREQUENCY]));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fragment length file: {}", e.toString());
        }

    }

}
