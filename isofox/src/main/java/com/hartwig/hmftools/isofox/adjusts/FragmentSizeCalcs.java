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
import static com.hartwig.hmftools.isofox.ChromosomeGeneTask.findNextOverlappingGenes;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
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

public class FragmentSizeCalcs
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final List<int[]> mFragmentLengths;
    private final List<int[]> mFragmentLengthsByGene;
    private int mMaxReadLength;

    private BufferedWriter mGeneWriter;

    private String mCurrentGenes;
    private final long[] mCurrentGenesRange;
    private List<TranscriptData> mCurrentTransDataList;
    private int mCurrentFragmentCount;
    private int mTotalFragmentCount;
    private int mProcessedFragments;
    private final FragmentTracker mFragmentTracker;

    private static final int MIN_GENE_LENGTH = 1000;
    private static final int MAX_GENE_LENGTH = 1000000;
    private static final int MAX_GENE_TRANS = 50;
    private static final int MAX_TRAN_EXONS = 20;
    private static final int MAX_GENE_FRAGMENT_COUNT = 5000; // to avoid impact of highly enriched genes

    private PerformanceCounter mPerfCounter;

    public FragmentSizeCalcs(final IsofoxConfig config, final EnsemblDataCache geneTransCache, final BufferedWriter writer)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, true);

        mCurrentGenes = "";
        mCurrentGenesRange = new long[SE_PAIR];
        mCurrentTransDataList = Lists.newArrayList();
        mCurrentFragmentCount = 0;
        mTotalFragmentCount = 0;
        mProcessedFragments = 0;
        mMaxReadLength = 0;
        mFragmentTracker = new FragmentTracker();

        mGeneWriter = writer;

        mFragmentLengths = Lists.newArrayList();
        mFragmentLengthsByGene = Lists.newArrayList();
        mPerfCounter = new PerformanceCounter("FragLengthDist");
    }

    public final List<int[]> getFragmentLengths() { return mFragmentLengths; }
    public final int getMaxReadLength() { return mMaxReadLength; }

    public void calcSampleFragmentSize(final String chromosome, final List<EnsemblGeneData> geneDataList, int requiredFragCount)
    {
        if (geneDataList == null || geneDataList.isEmpty())
            return;

        if (requiredFragCount == 0)
            return;

        // walk through each chromosome, taking groups of overlapping genes together
        ISF_LOGGER.info("calculating fragment size for chromosome({}) geneCount({})", chromosome, geneDataList.size());

        List<long[]> excludedRegions = generateExcludedRegions(chromosome);

        final List<EnsemblGeneData> overlappingGenes = Lists.newArrayList();
        int currentGeneIndex = 0;
        int nextLogCount = 100;

        while(currentGeneIndex < geneDataList.size())
        {
            currentGeneIndex = findNextOverlappingGenes(geneDataList, currentGeneIndex, overlappingGenes);

            mCurrentTransDataList.clear();
            mFragmentTracker.clear();
            mCurrentGenesRange[SE_START] = 0;
            mCurrentGenesRange[SE_END] = 0;

            for (int i = 0; i < overlappingGenes.size(); ++i)
            {
                EnsemblGeneData geneData = overlappingGenes.get(i);

                if (mConfig.EnrichedGeneIds.contains(geneData.GeneId))
                {
                    mCurrentTransDataList.clear();
                    break;
                }

                mCurrentGenesRange[SE_START] = i == 0 ? geneData.GeneStart : min(geneData.GeneStart, mCurrentGenesRange[SE_START]);
                mCurrentGenesRange[SE_END] = i == 0 ? geneData.GeneEnd : max(geneData.GeneEnd, mCurrentGenesRange[SE_END]);

                mCurrentTransDataList.addAll(mGeneTransCache.getTranscripts(geneData.GeneId).stream()
                        .filter(x -> x.exons().size() <= MAX_TRAN_EXONS).collect(Collectors.toList()));
            }

            if (mCurrentTransDataList.isEmpty() || mCurrentTransDataList.size() > MAX_GENE_TRANS)
                continue;

            long geneLength = mCurrentGenesRange[SE_END] - mCurrentGenesRange[SE_START];

            if (geneLength < MIN_GENE_LENGTH || geneLength > MAX_GENE_LENGTH)
                continue;

            if(excludedRegions.stream().anyMatch(x -> positionsOverlap(x[SE_START], x[SE_END], mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END])))
                continue;

            if(currentGeneIndex >= nextLogCount)
            {
                nextLogCount += 100;
                ISF_LOGGER.debug("chromosome({}) processed {} genes, fragCount({}) totalReads({})",
                        chromosome, currentGeneIndex, mProcessedFragments, mTotalFragmentCount);
            }

            mPerfCounter.start();

            mCurrentFragmentCount = 0;
            mCurrentGenes = overlappingGenes.get(0).GeneName;
            List<GenomeRegion> regions = Lists.newArrayList(GenomeRegions.create(chromosome, mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END]));

            mBamSlicer.slice(mSamReader, regions, this::processBamRead);

            mPerfCounter.stop();

            if(mConfig.WriteFragmentLengthsByGene)
            {
                String genesName = overlappingGenes.get(0).GeneName;
                for(int i = 1; i < min(overlappingGenes.size(), 10); ++i)
                {
                    genesName += ";" + overlappingGenes.get(i).GeneName;
                }

                writeGeneFragmentLengths(mGeneWriter, mFragmentLengthsByGene, genesName, overlappingGenes.size(), chromosome, mCurrentGenesRange);
                mFragmentLengthsByGene.clear();
            }

            if (mProcessedFragments >= requiredFragCount)
            {
                ISF_LOGGER.debug("chromosome({}) max fragment length samples reached: {}", chromosome, mProcessedFragments);
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

    private void processBamRead(@NotNull final SAMRecord read)
    {
        // cull invalid reads without waiting for the paired read
        if(!isCandidateRecord(read))
            return;

        mMaxReadLength = max(mMaxReadLength, read.getReadLength());

        // reads cannot cover any part of an exon

        long posStart = read.getStart();
        long posEnd = read.getEnd();

        for(final TranscriptData transData : mCurrentTransDataList)
        {
            if(transData.exons().stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x.ExonStart, x.ExonEnd)))
                return;
        }

        final SAMRecord otherRead = (SAMRecord)mFragmentTracker.checkRead(read.getReadName(), read);

        if(otherRead == null)
            return;

        ++mCurrentFragmentCount;
        ++mTotalFragmentCount;

        if(mCurrentFragmentCount >= MAX_GENE_FRAGMENT_COUNT)
        {
            ISF_LOGGER.trace("currentGenes({}) reached max fragment count", mCurrentGenes);
            mBamSlicer.haltProcessing();
            return;
        }

        addFragmentLength(read, mFragmentLengths);

        if(mConfig.WriteFragmentLengthsByGene)
        {
            addFragmentLength(read, mFragmentLengthsByGene);

            if(abs(read.getInferredInsertSize()) > FRAGMENT_LENGTH_CAP)
            {
                // Time,Genes,Chromosome,Read1Start,Read1End,Read2Start,Read2End,InsertSize,Cigar1,Cigar2,ReadId
                ISF_LOGGER.info(String.format("LONG_FRAG:%s,%s,%d,%d,%d,%d,%d,%s,%s,%s",
                        mCurrentGenes, read.getReferenceName(), read.getStart(), read.getEnd(), otherRead.getStart(), otherRead.getEnd(),
                        abs(read.getInferredInsertSize()), read.getCigar().toString(), otherRead.getCigar().toString(), read.getReadName()));
            }
        }
    }

    private boolean isCandidateRecord(final SAMRecord record)
    {
        // ignore translocations and inversions
        if(!record.getMateReferenceName().equals(record.getReferenceName()) || record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag())
            return false;

        // ignore split or unmapped reads
        if(record.getCigar() == null || record.getCigar().containsOperator(CigarOperator.N) || !record.getCigar().containsOperator(CigarOperator.M))
            return false;

        // both reads must fall in the current gene
        long otherStartPos = record.getMateAlignmentStart();
        if(!positionWithin(otherStartPos, mCurrentGenesRange[SE_START], mCurrentGenesRange[SE_END]))
            return false;

        return true;
    }

    private static final int FRAGMENT_LENGTH_CAP = 10000;

    private int getLengthBucket(int fragmentLength)
    {
        if(fragmentLength < 1000)
            return fragmentLength;

        if(fragmentLength < 3000)
            return 10 * (int)round(fragmentLength/10.0);

        if(fragmentLength < FRAGMENT_LENGTH_CAP)
            return 100 * (int)round(fragmentLength/100.0);

        return FRAGMENT_LENGTH_CAP;
    }

    private void addFragmentLength(final SAMRecord record, final List<int[]> fragmentLengths)
    {
        int fragmentLength = getLengthBucket(abs(record.getInferredInsertSize()));

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

    public static BufferedWriter createGeneFragmentLengthWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.WriteFragmentLengthsByGene ?
                    config.formOutputFile("frag_length_by_gene.csv") : config.formOutputFile("frag_length.csv");

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
            final String chromosome, final long[] genesRegion)
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

                writer.write(String.format(",%d,%d", fragLengthCount[FL_LENGTH], fragLengthCount[FL_FREQUENCY]));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene fragment length file: {}", e.toString());
        }
    }

    public static void writeFragmentLengths(final IsofoxConfig config, final List<int[]> fragmentLengths)
    {
        if(fragmentLengths.isEmpty())
            return;

        try
        {
            final String outputFileName = config.formOutputFile("frag_length.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("FragmentLength,Count");
            writer.newLine();

            for (final int[] fragLengthCount : fragmentLengths)
            {
                writer.write(String.format("%d,%d", fragLengthCount[FL_LENGTH], fragLengthCount[FL_FREQUENCY]));
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
