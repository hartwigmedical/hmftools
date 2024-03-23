package com.hartwig.hmftools.isofox.adjusts;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.ChromosomeTaskExecutor.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.SINGLE_MAP_QUALITY;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.FragmentTracker;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// calculate fragment length distribution for a sample
// this can be done either independently from fragment counting or lengths can be registering during that process

public class FragmentSizeCalcs implements Callable
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private String mChromosome;
    private final List<GeneData> mGeneDataList;
    private int mRequiredFragCount;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final List<FragmentSize> mFragmentLengths;
    private final List<FragmentSize> mFragmentLengthsByGene;
    private int mMaxReadLength;

    private BufferedWriter mGeneWriter;

    private String mCurrentGenes;
    private final int[] mCurrentGenesRange;
    private List<TranscriptData> mCurrentTransDataList;
    private int mCurrentFragmentCount;
    private int mTotalFragmentCount;
    private int mProcessedFragments;
    private final FragmentTracker mFragmentTracker;

    private static final int MIN_GENE_LENGTH = 1000;
    private static final int MAX_GENE_LENGTH = 1000000;
    private static final int MAX_GENE_TRANS = 50;
    private static final int MAX_GENE_FRAGMENT_COUNT = 5000; // to avoid impact of highly enriched genes
    private static final int MAX_FRAGMENT_LENGTH = 5000; // ignored beyond this

    private PerformanceCounter mPerfCounter;

    public FragmentSizeCalcs(
            final IsofoxConfig config, final EnsemblDataCache geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(SINGLE_MAP_QUALITY, false, false, false);

        mCurrentGenes = "";
        mCurrentGenesRange = new int[SE_PAIR];
        mCurrentTransDataList = Lists.newArrayList();
        mCurrentFragmentCount = 0;
        mTotalFragmentCount = 0;
        mProcessedFragments = 0;
        mMaxReadLength = 0;
        mFragmentTracker = new FragmentTracker();

        mGeneWriter = writer;

        mGeneDataList = Lists.newArrayList();
        mFragmentLengths = Lists.newArrayList();
        mFragmentLengthsByGene = Lists.newArrayList();

        mChromosome = "";
        mRequiredFragCount = 0;

        mPerfCounter = new PerformanceCounter("FragLengthDist");
    }

    public final List<FragmentSize> getFragmentLengths() { return mFragmentLengths; }
    public final int getMaxReadLength() { return mMaxReadLength; }
    public final PerformanceCounter getPerformanceCounter() { return mPerfCounter; }

    public void initialise(final String chromosome, final List<GeneData> geneDataList, int requiredFragCount)
    {
        mChromosome = chromosome;
        mGeneDataList.clear();
        mGeneDataList.addAll(geneDataList);
        mRequiredFragCount = requiredFragCount;
    }

    @Override
    public Long call()
    {
        calcSampleFragmentSize();
        return (long)0;
    }

    private void calcSampleFragmentSize()
    {
        if(mGeneDataList.isEmpty() || mRequiredFragCount == 0)
        {
            ISF_LOGGER.error("chromosome({}) fragment size uninitialised", mChromosome);
            return;
        }

        // walk through each chromosome, taking groups of overlapping genes together
        ISF_LOGGER.debug("calculating fragment size for chromosome({}) geneCount({})", mChromosome, mGeneDataList.size());

        final List<GeneData> overlappingGenes = Lists.newArrayList();
        int currentGeneIndex = 0;
        int nextLogCount = 100;

        while(currentGeneIndex < mGeneDataList.size())
        {
            currentGeneIndex = findNextOverlappingGenes(mGeneDataList, currentGeneIndex, overlappingGenes);

            if(overlappingGenes.stream().anyMatch(x -> mConfig.Filters.EnrichedGeneIds.contains(x.GeneId)))
                continue;

            mCurrentTransDataList.clear();

            mFragmentTracker.clear();
            mCurrentGenesRange[SE_START] = 0;
            mCurrentGenesRange[SE_END] = 0;

            for(int i = 0; i < overlappingGenes.size(); ++i)
            {
                GeneData geneData = overlappingGenes.get(i);

                mCurrentGenesRange[SE_START] = i == 0 ? geneData.GeneStart : min(geneData.GeneStart, mCurrentGenesRange[SE_START]);
                mCurrentGenesRange[SE_END] = i == 0 ? geneData.GeneEnd : max(geneData.GeneEnd, mCurrentGenesRange[SE_END]);

                mCurrentTransDataList.addAll(mGeneTransCache.getTranscripts(geneData.GeneId));
            }

            if(mCurrentTransDataList.isEmpty() || mCurrentTransDataList.size() > MAX_GENE_TRANS)
                continue;

            int geneLength = mCurrentGenesRange[SE_END] - mCurrentGenesRange[SE_START];

            if(geneLength < MIN_GENE_LENGTH || geneLength > MAX_GENE_LENGTH)
                continue;

            if(mChromosome.equals(mConfig.Filters.ExcludedRegion.Chromosome)
            && positionsOverlap(
                    mConfig.Filters.ExcludedRegion.start(), mConfig.Filters.ExcludedRegion.end(),
                    mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END]))
            {
                continue;
            }

            if(currentGeneIndex >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.debug("chromosome({}) processed {} genes, fragCount({}) totalReads({})",
                        mChromosome, currentGeneIndex, mProcessedFragments, mTotalFragmentCount);
            }

            mPerfCounter.start();

            mCurrentFragmentCount = 0;
            mCurrentGenes = overlappingGenes.get(0).GeneName;

            ChrBaseRegion sliceRegion = new ChrBaseRegion(mChromosome, mCurrentGenesRange);

            ISF_LOGGER.trace("chromosome({}) gene({} index={}) fragCount({}) nextRegion({})",
                    mChromosome, mCurrentGenes, currentGeneIndex, mProcessedFragments, sliceRegion);

            try
            {
                mBamSlicer.slice(mSamReader, sliceRegion, this::processBamRead);
            }
            catch(Exception e)
            {
                ISF_LOGGER.error("chromosome({}) geneIndex({}) fragCount({}) currentRegion({}) error slicing bam: {}",
                        mChromosome, currentGeneIndex, mProcessedFragments, sliceRegion, e.toString());
                return;
            }

            mPerfCounter.stop();

            if(mConfig.WriteFragmentLengthsByGene)
            {
                String genesName = overlappingGenes.get(0).GeneName;
                for(int i = 1; i < min(overlappingGenes.size(), 10); ++i)
                {
                    genesName += ITEM_DELIM + overlappingGenes.get(i).GeneName;
                }

                writeGeneFragmentLengths(mGeneWriter, mFragmentLengthsByGene, genesName, overlappingGenes.size(), mChromosome, mCurrentGenesRange);
                mFragmentLengthsByGene.clear();
            }

            if(mProcessedFragments >= mRequiredFragCount)
            {
                ISF_LOGGER.debug("chromosome({}) max fragment length samples reached: {}", mChromosome, mProcessedFragments);
                break;
            }
        }

        ISF_LOGGER.debug("chromosome({}) processing complete", mChromosome);

        mFragmentTracker.clear();
        mCurrentTransDataList.clear();

        if(mConfig.WriteFragmentLengthsByGene)
        {
            final List<Double> fragLengths = FragmentSizeCalcs.calcPercentileData(mFragmentLengths, Lists.newArrayList(0.05, 0.5, 0.95));

            if(!fragLengths.isEmpty())
            {
                ISF_LOGGER.info("chromosome({}) frag lengths 5th({}) 50th({}) 95th({})",
                        mChromosome, fragLengths.get(0), fragLengths.get(1), fragLengths.get(2));
            }
        }
    }

    private void processBamRead(@NotNull final SAMRecord read)
    {
        // cull invalid reads without waiting for the paired read
        if(!isCandidateRecord(read))
            return;

        ++mCurrentFragmentCount;
        ++mTotalFragmentCount;

        if(mCurrentFragmentCount >= MAX_GENE_FRAGMENT_COUNT)
        {
            ISF_LOGGER.debug("currentGenes({}) reached max fragment count", mCurrentGenes);
            mBamSlicer.haltProcessing();
            return;
        }

        mMaxReadLength = max(mMaxReadLength, read.getReadLength());

        final SAMRecord otherRead = (SAMRecord)mFragmentTracker.checkRead(read.getReadName(), read);

        if(otherRead == null)
            return;

        addFragmentLength(read, mFragmentLengths);

        if(mConfig.WriteFragmentLengthsByGene)
        {
            addFragmentLength(read, mFragmentLengthsByGene);
        }
    }

    private boolean isCandidateRecord(final SAMRecord record)
    {
        int fragmentLength = abs(record.getInferredInsertSize());
        if(fragmentLength > MAX_FRAGMENT_LENGTH)
            return false;

        // ignore translocations and inversions
        if(!record.getMateReferenceName().equals(record.getReferenceName()) || mateNegativeStrand(record) == record.getReadNegativeStrandFlag())
            return false;

        // ignore split and soft-clipped reads above the read length
        if(record.getCigar().containsOperator(CigarOperator.N) || !record.getCigar().containsOperator(CigarOperator.M))
            return false;

        int readLength = max(mMaxReadLength, mConfig.ReadLength);

        if(readLength > 0 && fragmentLength > readLength && record.getCigar().containsOperator(CigarOperator.S))
        {
            // permit small soft-clips up to a point and then none for longer fragments
            if(record.getCigar().getCigarElements().stream().anyMatch(x -> x.getOperator() == CigarOperator.S && x.getLength() > 2))
                return false;
        }

        // both reads must fall in the current gene
        int otherStartPos = record.getMateAlignmentStart();
        if(!positionWithin(otherStartPos, mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END]))
            return false;

        // reads cannot cover any part of an exon
        int posStart = record.getStart();
        int posEnd = record.getEnd();

        for(final TranscriptData transData : mCurrentTransDataList)
        {
            if(transData.exons().stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.Start, x.End)))
                return false;
        }

        return true;
    }

    private void addFragmentLength(final SAMRecord record, final List<FragmentSize> fragmentLengths)
    {
        int fragmentLength = getLengthBucket(abs(record.getInferredInsertSize()));

        if(fragmentLength == 0)
            return;

        int index = 0;
        boolean exists = false;
        while(index < fragmentLengths.size())
        {
            final FragmentSize fragLengthCount = fragmentLengths.get(index);

            if(fragLengthCount.Length < fragmentLength)
            {
                ++index;
                continue;
            }

            if(fragLengthCount.Length == fragmentLength)
            {
                ++fragLengthCount.Frequency;
                exists = true;
            }

            break;
        }

        if(!exists)
        {
            FragmentSize newFragLengthCount = new FragmentSize(fragmentLength, 1);
            fragmentLengths.add(index, newFragLengthCount);
        }

        ++mProcessedFragments;
    }

    private int getLengthBucket(int fragmentLength)
    {
        // round to nearest unit up to 1000, then 10s up to 3000 then 100s
        if(fragmentLength < 1000)
            return fragmentLength;

        if(fragmentLength < 3000)
            return 10 * (int)round(fragmentLength/10.0);

        return 100 * (int)round(fragmentLength/100.0);
    }

    public static void setConfigFragmentLengthData(final IsofoxConfig config, final List<FragmentSize> fragmentLengths)
    {
        final List<FragmentSize> lengthFrequencies = config.FragmentSizeData;

        int currentRangeMin = 0;
        int currentRangeMax = 0;

        for(int i = 0; i < lengthFrequencies.size(); ++i)
        {
            FragmentSize lengthFrequency = lengthFrequencies.get(i);

            currentRangeMin = (i == 0) ? 0 : currentRangeMax + 1;

            if(i == lengthFrequencies.size() - 1)
            {
                currentRangeMax = config.MaxFragmentLength - 1;
            }
            else
            {
                FragmentSize nextLengthFrequency = lengthFrequencies.get(i + 1);
                currentRangeMax = (lengthFrequency.Length + nextLengthFrequency.Length) / 2;
            }

            int lengthCount = 0;

            for (final FragmentSize fragLengthCount : fragmentLengths)
            {
                if(fragLengthCount.Length >= currentRangeMin && fragLengthCount.Length <= currentRangeMax)
                {
                    lengthCount += fragLengthCount.Frequency;
                }
            }

            lengthFrequency.Frequency = lengthCount;
            ISF_LOGGER.info("fragmentLength({}) frequency({})", lengthFrequency.Length, lengthCount);
        }
    }

    public static List<Double> calcPercentileData(final List<FragmentSize> fragmentLengths, final List<Double> percentiles)
    {
        final List<Double> percentileLengths = Lists.newArrayList();

        double totalFragments = fragmentLengths.stream().mapToLong(x -> x.Frequency).sum();

        for(Double percentile : percentiles)
        {
            int currentTotal = 0;
            int prevLength = 0;

            for (final FragmentSize fragLengthCount : fragmentLengths)
            {
                double nextPercTotal = (currentTotal + fragLengthCount.Frequency) / totalFragments;

                if(nextPercTotal >= percentile)
                {
                    double percLength = prevLength > 0 ? (prevLength + fragLengthCount.Length) * 0.5 : fragLengthCount.Length;
                    percentileLengths.add(percLength);
                    break;
                }

                currentTotal += fragLengthCount.Frequency;
                prevLength = fragLengthCount.Length;
            }
        }

        return percentileLengths;
    }

    public static void mergeData(final List<FragmentSize> fragmentLengths, final FragmentSizeCalcs other)
    {
        for(FragmentSize otherLengthData : other.getFragmentLengths())
        {
            int fragmentLength = otherLengthData.Length;

            int index = 0;
            boolean exists = false;
            while (index < fragmentLengths.size())
            {
                final FragmentSize fragLengthCount = fragmentLengths.get(index);

                if(fragLengthCount.Length < fragmentLength)
                {
                    ++index;
                    continue;
                }

                if(fragLengthCount.Length == fragmentLength)
                {
                    fragLengthCount.Frequency += otherLengthData.Frequency;
                    exists = true;
                }

                break;
            }

            if(!exists)
            {
                FragmentSize newLengthData = new FragmentSize(fragmentLength, otherLengthData.Frequency);
                fragmentLengths.add(index, newLengthData);
            }
        }
    }

    public static BufferedWriter createGeneFragmentLengthWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("frag_length_by_gene.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneNames,GeneCount,Chromosome,GenesStart,GenesEnd,");
            writer.write("FragmentLength,Count");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene fragment length file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeGeneFragmentLengths(
            final BufferedWriter writer, final List<FragmentSize> fragmentLengths, final String geneNames, int geneCount,
            final String chromosome, final int[] genesRegion)
    {
        if(writer == null)
            return;

        if(fragmentLengths.isEmpty())
            return;

        try
        {
            for (FragmentSize fragLengthCount : fragmentLengths)
            {
                writer.write(String.format("%s,%d,%s,%d,%d",
                        geneNames, geneCount, chromosome, genesRegion[SE_START], genesRegion[SE_END]));

                writer.write(String.format(",%d,%d", fragLengthCount.Length, fragLengthCount.Frequency));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene fragment length file: {}", e.toString());
        }
    }

    public static void writeFragmentLengths(final IsofoxConfig config, final List<FragmentSize> fragmentLengths)
    {
        if(fragmentLengths.isEmpty())
            return;

        try
        {
            final String outputFileName = config.formOutputFile("frag_length.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("FragmentLength,Count");
            writer.newLine();

            for (final FragmentSize fragLengthData : fragmentLengths)
            {
                writer.write(String.format("%d,%d", fragLengthData.Length, fragLengthData.Frequency));
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write fragment length file: {}", e.toString());
        }
    }

}
