package com.hartwig.hmftools.isofox.adjusts;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.BamFragmentReader.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.IsofoxConstants.ENRICHED_GENE_BUFFER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.utils.sv.SvRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
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
    private final List<EnsemblGeneData> mGeneDataList;
    private int mRequiredFragCount;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final List<int[]> mFragmentLengths;
    private final List<int[]> mFragmentLengthsByGene;
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

    private static final int FD_LENGTH = 0;
    private static final int FD_FREQUENCY = 1;

    private final List<Integer> mSoftClipLengthBuckets;

    private PerformanceCounter mPerfCounter;

    public FragmentSizeCalcs(
            final IsofoxConfig config, final EnsemblDataCache geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true, true);

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

        mSoftClipLengthBuckets = Lists.newArrayList(1,2,3,5,10,25,50);
        
        mChromosome = "";
        mRequiredFragCount = 0;

        mPerfCounter = new PerformanceCounter("FragLengthDist");
    }

    public final List<int[]> getFragmentLengths() { return mFragmentLengths; }
    public final int getMaxReadLength() { return mMaxReadLength; }
    public final PerformanceCounter getPerformanceCounter() { return mPerfCounter; }
    public final List<Integer> getSoftClipLengthBuckets() { return mSoftClipLengthBuckets; }
    
    public void initialise(final String chromosome, final List<EnsemblGeneData> geneDataList, int requiredFragCount)
    {
        mChromosome = chromosome;
        mGeneDataList.clear();
        mGeneDataList.addAll(geneDataList);
        mRequiredFragCount = requiredFragCount;
    }

    @Override
    public Long call()
    {
        mPerfCounter.start();
        calcSampleFragmentSize();
        mPerfCounter.stop();
        return (long)0;
    }

    private void calcSampleFragmentSize()
    {
        if (mGeneDataList.isEmpty() || mRequiredFragCount == 0)
        {
            ISF_LOGGER.error("chromosome({}) fragment size uninitialised", mChromosome);
            return;
        }

        // walk through each chromosome, taking groups of overlapping genes together
        ISF_LOGGER.info("calculating fragment size for chromosome({}) geneCount({})", mChromosome, mGeneDataList.size());

        List<int[]> excludedRegions = generateExcludedRegions();

        final List<EnsemblGeneData> overlappingGenes = Lists.newArrayList();
        int currentGeneIndex = 0;
        int nextLogCount = 100;

        while(currentGeneIndex < mGeneDataList.size())
        {
            currentGeneIndex = findNextOverlappingGenes(mGeneDataList, currentGeneIndex, overlappingGenes);

            if(overlappingGenes.stream().anyMatch(x -> mConfig.EnrichedGeneIds.contains(x.GeneId)))
                continue;

            mCurrentTransDataList.clear();

            mFragmentTracker.clear();
            mCurrentGenesRange[SE_START] = 0;
            mCurrentGenesRange[SE_END] = 0;

            for (int i = 0; i < overlappingGenes.size(); ++i)
            {
                EnsemblGeneData geneData = overlappingGenes.get(i);

                mCurrentGenesRange[SE_START] = i == 0 ? geneData.GeneStart : min(geneData.GeneStart, mCurrentGenesRange[SE_START]);
                mCurrentGenesRange[SE_END] = i == 0 ? geneData.GeneEnd : max(geneData.GeneEnd, mCurrentGenesRange[SE_END]);

                mCurrentTransDataList.addAll(mGeneTransCache.getTranscripts(geneData.GeneId));
            }

            if (mCurrentTransDataList.isEmpty() || mCurrentTransDataList.size() > MAX_GENE_TRANS)
                continue;

            int geneLength = mCurrentGenesRange[SE_END] - mCurrentGenesRange[SE_START];

            if (geneLength < MIN_GENE_LENGTH || geneLength > MAX_GENE_LENGTH)
                continue;

            if(excludedRegions.stream().anyMatch(x -> positionsOverlap(x[SE_START], x[SE_END], mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END])))
                continue;

            if(currentGeneIndex >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.debug("chromosome({}) processed {} genes, fragCount({}) totalReads({})",
                        mChromosome, currentGeneIndex, mProcessedFragments, mTotalFragmentCount);
            }

            mPerfCounter.start();

            mCurrentFragmentCount = 0;
            mCurrentGenes = overlappingGenes.get(0).GeneName;

            final List<SvRegion> regions = Lists.newArrayList(new SvRegion(mChromosome, mCurrentGenesRange));

            ISF_LOGGER.trace("chromosome({}) gene({} index={}) fragCount({}) nextRegion({})",
                    mChromosome, mCurrentGenes, currentGeneIndex, mProcessedFragments, regions.get(0).toString());

            try
            {
                mBamSlicer.slice(mSamReader, regions, this::processBamRead);
            }
            catch(Exception e)
            {
                ISF_LOGGER.error("chromosome({}) geneIndex({}) fragCount({}) currentRegion({}) error slicing bam: {}",
                        mChromosome, currentGeneIndex, mProcessedFragments, regions.get(0).toString(), e.toString());
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

            if (mProcessedFragments >= mRequiredFragCount)
            {
                ISF_LOGGER.debug("chromosome({}) max fragment length samples reached: {}", mChromosome, mProcessedFragments);
                break;
            }
        }

        ISF_LOGGER.debug("chromosome({}) processing complete", mChromosome);
    }

    private List<int[]> generateExcludedRegions()
    {
        // create a buffer around the enriched gene to avoid excessive reads in this vicinity
        final List<int[]> excludedRegions = Lists.newArrayList();
        for(final String geneId : mConfig.EnrichedGeneIds)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataById(geneId);
            if(geneData == null)
                continue;

            if(geneData.Chromosome.equals(mChromosome))
            {
                excludedRegions.add(new int[] { geneData.GeneStart - ENRICHED_GENE_BUFFER, geneData.GeneEnd + ENRICHED_GENE_BUFFER});
            }
        }

        return excludedRegions;
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
        if(!record.getMateReferenceName().equals(record.getReferenceName()) || record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag())
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
            if(transData.exons().stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.ExonStart, x.ExonEnd)))
                return false;
        }

        return true;
    }

    private int getLengthBucket(int fragmentLength)
    {
        if(fragmentLength < 1000)
            return fragmentLength;

        if(fragmentLength < 3000)
            return 10 * (int)round(fragmentLength/10.0);

        return 100 * (int)round(fragmentLength/100.0);
    }

    private int getSoftClipBucketIndex(int scLength)
    {
        if(scLength == 0)
            return -1;

        for(int i = 0; i < mSoftClipLengthBuckets.size(); ++i)
        {
            if(scLength <= mSoftClipLengthBuckets.get(i))
                return i;
        }

        return mSoftClipLengthBuckets.size() - 1;
    }

    private void addFragmentLength(final SAMRecord record, final List<int[]> fragmentLengths)
    {
        int fragmentLength = getLengthBucket(abs(record.getInferredInsertSize()));

        if(fragmentLength == 0)
            return;

        int scLength = record.getCigar().getCigarElements().stream()
                .filter(x -> x.getOperator() == CigarOperator.S).mapToInt(x -> x.getLength()).sum();

        int scIndex = getSoftClipBucketIndex(scLength);

        int index = 0;
        boolean exists = false;
        while(index < fragmentLengths.size())
        {
            final int[] fragLengthCount = fragmentLengths.get(index);

            if(fragLengthCount[FD_LENGTH] < fragmentLength)
            {
                ++index;
                continue;
            }

            if(fragLengthCount[FD_LENGTH] == fragmentLength)
            {
                ++fragLengthCount[FD_FREQUENCY];

                if(scIndex >= 0)
                    ++fragLengthCount[scIndex];

                exists = true;
            }

            break;
        }

        if(!exists)
        {
            int[] newFragLengthCount = new int[2 + mSoftClipLengthBuckets.size()];
            newFragLengthCount[FD_LENGTH] = fragmentLength;
            newFragLengthCount[FD_FREQUENCY] = 1;

            if(scIndex >= 0)
                ++newFragLengthCount[scIndex];

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
                currentRangeMax = (lengthFrequency[FD_LENGTH] + nextLengthFrequency[FD_LENGTH]) / 2;
            }

            int lengthCount = 0;

            for (final int[] fragLengthCount : fragmentLengths)
            {
                if(fragLengthCount[FD_LENGTH] >= currentRangeMin && fragLengthCount[FD_LENGTH] <= currentRangeMax)
                {
                    lengthCount += fragLengthCount[FD_FREQUENCY];
                }
            }

            lengthFrequency[FD_FREQUENCY] = lengthCount;
            ISF_LOGGER.info("fragmentLength({}) frequency({})", lengthFrequency[FD_LENGTH], lengthCount);
        }
    }

    public static List<Double> calcPercentileData(final List<int[]> fragmentLengths, final List<Double> percentiles)
    {
        final List<Double> percentileLengths = Lists.newArrayList();

        double totalFragments = fragmentLengths.stream().mapToLong(x -> x[FD_FREQUENCY]).sum();

        for(Double percentile : percentiles)
        {
            int currentTotal = 0;
            int prevLength = 0;

            for (final int[] fragLengthCount : fragmentLengths)
            {
                double nextPercTotal = (currentTotal + fragLengthCount[FD_FREQUENCY]) / totalFragments;

                if(nextPercTotal >= percentile)
                {
                    double percLength = prevLength > 0 ? (prevLength + fragLengthCount[FD_LENGTH]) * 0.5 : fragLengthCount[FD_LENGTH];
                    percentileLengths.add(percLength);
                    break;
                }

                currentTotal += fragLengthCount[FD_FREQUENCY];
                prevLength = fragLengthCount[FD_LENGTH];
            }
        }

        return percentileLengths;
    }

    public static void mergeData(final List<int[]> fragmentLengths, final FragmentSizeCalcs other)
    {
        for(final int[] otherLengthData : other.getFragmentLengths())
        {
            int fragmentLength = otherLengthData[FD_LENGTH];

            int index = 0;
            boolean exists = false;
            while (index < fragmentLengths.size())
            {
                final int[] fragLengthCount = fragmentLengths.get(index);

                if (fragLengthCount[FD_LENGTH] < fragmentLength)
                {
                    ++index;
                    continue;
                }

                if (fragLengthCount[FD_LENGTH] == fragmentLength)
                {
                    for(int i = FD_FREQUENCY; i < fragLengthCount.length; ++i)
                    {
                        fragLengthCount[i] += otherLengthData[i];
                    }
                    exists = true;
                }

                break;
            }

            if (!exists)
            {
                int[] newLengthData = new int[otherLengthData.length];

                for(int i = 0; i < otherLengthData.length; ++i)
                {
                    newLengthData[i] = otherLengthData[i];
                }
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
            final BufferedWriter writer, final List<int[]> fragmentLengths, final String geneNames, int geneCount,
            final String chromosome, final int[] genesRegion)
    {
        if(writer == null)
            return;

        if(fragmentLengths.isEmpty())
            return;

        try
        {
            for (final int[] fragLengthCount : fragmentLengths)
            {
                writer.write(String.format("%s,%d,%s,%d,%d",
                        geneNames, geneCount, chromosome, genesRegion[SE_START], genesRegion[SE_END]));

                writer.write(String.format(",%d,%d", fragLengthCount[FD_LENGTH], fragLengthCount[FD_FREQUENCY]));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene fragment length file: {}", e.toString());
        }
    }

    public static void writeFragmentLengths(final IsofoxConfig config, final List<int[]> fragmentLengths, final List<Integer> scLengths)
    {
        if(fragmentLengths.isEmpty())
            return;

        try
        {
            final String outputFileName = config.formOutputFile("frag_length.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("FragmentLength,Count");

            for(Integer scLength : scLengths)
            {
                writer.write(String.format(",Sc%d", scLength));
            }

            writer.newLine();

            for (final int[] fragLengthData : fragmentLengths)
            {
                writer.write(String.format("%d,%d", fragLengthData[FD_LENGTH], fragLengthData[FD_FREQUENCY]));

                for(int i = 2; i < fragLengthData.length; ++i)
                {
                    writer.write(String.format(",%d", fragLengthData[i]));
                }

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
