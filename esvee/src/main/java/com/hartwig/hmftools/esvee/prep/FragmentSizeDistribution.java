package com.hartwig.hmftools.esvee.prep;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.inferredInsertSize;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FragmentLengthBounds.INVALID;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FRAG_LENGTH_1_STD_DEV_PERCENTILE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FRAG_LENGTH_DIST_MAX_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FRAG_LENGTH_DIST_MIN_QUAL;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FRAG_LENGTH_DIST_PERCENTILE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.FRAG_LENGTH_DIST_SAMPLE_SIZE;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.FRAGMENT_LENGTH_DIST;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.types.LengthFrequency;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FragmentSizeDistribution
{
    private final PrepConfig mConfig;
    private final List<LengthFrequency> mLengthFrequencies;
    private int mMaxReadLength;
    private boolean mHasPairedReads;

    public FragmentSizeDistribution(final PrepConfig config)
    {
        mConfig = config;
        mLengthFrequencies = Lists.newArrayList();
        mMaxReadLength = 0;
        mHasPairedReads = false;
    }
    
    public void run()
    {
        SV_LOGGER.info("calculating fragment size distribution");

        // for(HumanChromosome chromosome : HumanChromosome.values())
        List<String> sampledChromosomes = Lists.newArrayList();

        if(mConfig.SpecificChrRegions.hasFilters())
        {
            sampledChromosomes.addAll(mConfig.SpecificChrRegions.Chromosomes);
        }
        else
        {
            for(int i = 1; i <= 10; ++i)
            {
                sampledChromosomes.add(mConfig.RefGenVersion.versionedChromosome(String.valueOf(i)));
            }
        }

        List<ChromosomeTask> chrTasks = sampledChromosomes.stream().map(x -> new ChromosomeTask(x)).collect(Collectors.toList());

        final List<Callable<Long>> callableList = chrTasks.stream().collect(Collectors.toList());
        boolean validExecution = TaskExecutor.executeTasks(callableList, mConfig.Threads);

        if(!validExecution)
            return;

        // merge results from all chromosomes
        for(final ChromosomeTask chrTask : chrTasks)
        {
            mergeDistributions(mLengthFrequencies, chrTask.lengthFrequencies());
            mMaxReadLength = max(mMaxReadLength, chrTask.maxReadLength());
            mHasPairedReads |= chrTask.hasPairedRead();
        }

        SV_LOGGER.info("maxReadLength({}) {}", mMaxReadLength, !mHasPairedReads ? "unpaired reads" : "");

        if(!mLengthFrequencies.isEmpty())
        {
            int minLength = mLengthFrequencies.get(0).Length;
            int maxLength = mLengthFrequencies.get(mLengthFrequencies.size() - 1).Length;

            SV_LOGGER.debug("fragment size distribution complete: min({}) max({})", minLength, maxLength);
        }
        else
        {
            SV_LOGGER.debug("no fragment lengths recorded");

        }

        if(mConfig.WriteTypes.contains(FRAGMENT_LENGTH_DIST))
            writeDistribution();
    }

    public int maxReadLength() { return mMaxReadLength; }
    public boolean hasPairedReads() { return mHasPairedReads; }

    public FragmentLengthBounds calculateFragmentLengthBounds() { return calculateFragmentLengthBounds(mLengthFrequencies); }

    public static FragmentLengthBounds calculateFragmentLengthBounds(final List<LengthFrequency> lengthFrequencies)
    {
        if(lengthFrequencies.isEmpty())
            return INVALID;

        int totalFragments = lengthFrequencies.stream().mapToInt(x -> x.Frequency).sum();
        int cumulativeTotal = 0;
        int requiredMinTotal = (int)floor(totalFragments * (1 - FRAG_LENGTH_DIST_PERCENTILE));
        int requiredStdDevTotal = (int)floor(totalFragments * FRAG_LENGTH_1_STD_DEV_PERCENTILE);
        int requiredMaxTotal = (int)floor(totalFragments * FRAG_LENGTH_DIST_PERCENTILE);
        int medianFragment = totalFragments / 2;

        int lowerBound = 0;
        int upperBound = 0;
        int stdDevLength = 0;
        int median = 0;

        for(LengthFrequency lengthData : lengthFrequencies)
        {
            if(lowerBound == 0 && lengthData.Frequency + cumulativeTotal >= requiredMinTotal)
            {
                lowerBound = lengthData.Length;
            }

            if(stdDevLength == 0 && lengthData.Frequency + cumulativeTotal >= requiredStdDevTotal)
            {
                stdDevLength = lengthData.Length;
            }

            if(median == 0 && lengthData.Frequency + cumulativeTotal >= medianFragment)
            {
                median = lengthData.Length;
            }

            if(upperBound == 0 && lengthData.Frequency + cumulativeTotal >= requiredMaxTotal)
            {
                upperBound = lengthData.Length;
            }
            else
            {
                cumulativeTotal += lengthData.Frequency;
            }
        }

        double stdDeviation = median - stdDevLength;

        return new FragmentLengthBounds(lowerBound, upperBound, median, stdDeviation);
    }

    private void mergeDistributions(final List<LengthFrequency> lengthFrequencies, final List<LengthFrequency> otherFrequencies)
    {
        for(LengthFrequency otherLengthData : otherFrequencies)
        {
            int fragmentLength = otherLengthData.Length;

            int index = 0;
            boolean exists = false;
            while(index < lengthFrequencies.size())
            {
                final LengthFrequency fragLengthCount = lengthFrequencies.get(index);

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
                LengthFrequency newLengthData = new LengthFrequency(fragmentLength, otherLengthData.Frequency);
                lengthFrequencies.add(index, newLengthData);
            }
        }
    }

    private class ChromosomeTask implements Callable<Long>
    {
        private final String mChromosome;
        private int mProcessedReads;

        private final BamSlicer mBamSlicer;
        private final SamReader mSamReader;
        private final List<LengthFrequency> mLengthFrequencies;
        private int mMaxReadLength;
        private boolean mHasPairedReads;

        public ChromosomeTask(final String chromosome)
        {
            mChromosome = chromosome;
            mProcessedReads = 0;
            mMaxReadLength = 0;
            mHasPairedReads = false;

            // only run on the first sample if more than 1 are loaded
            mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.bamFile()));
            mBamSlicer = new BamSlicer(FRAG_LENGTH_DIST_MIN_QUAL, false, false, false);

            mLengthFrequencies = Lists.newArrayList();
        }

        public List<LengthFrequency> lengthFrequencies() { return mLengthFrequencies; }
        public int maxReadLength() { return mMaxReadLength; }
        public boolean hasPairedRead() { return mHasPairedReads; }

        @Override
        public Long call()
        {
            // slice a fixed region from each chromosome

            ChrBaseRegion region;

            if(!mConfig.SpecificChrRegions.Regions.isEmpty())
            {
                region = mConfig.SpecificChrRegions.Regions.stream().filter(x -> x.Chromosome.equals(mChromosome)).findFirst().orElse(null);

                if(region == null)
                    return (long)1;
            }
            else
            {
                region = new ChrBaseRegion(mChromosome, 1_000_000, 10_000_000);
            }

            mBamSlicer.slice(mSamReader, region, this::processBamRead);

            return (long)0;
        }

        private void processBamRead(final SAMRecord record)
        {
            // cull invalid reads without waiting for the paired read
            if(!isCandidateRecord(record))
                return;

            mMaxReadLength = max(mMaxReadLength, record.getReadBases().length);
            mHasPairedReads |= record.getReadPairedFlag();

            ++mProcessedReads;

            if(mProcessedReads >= FRAG_LENGTH_DIST_SAMPLE_SIZE)
            {
                SV_LOGGER.debug("chr({}) reached max fragment count", mChromosome);
                mBamSlicer.haltProcessing();
                return;
            }

            addFragmentLength(record);
        }

        private boolean isCandidateRecord(final SAMRecord record)
        {
            if(record.getDuplicateReadFlag())
                return false;

            int fragmentLength = abs(inferredInsertSize(record));
            if(fragmentLength > FRAG_LENGTH_DIST_MAX_LENGTH)
                return false;

            int readLength = record.getReadBases().length;

            if(!record.getReadPairedFlag() && fragmentLength == 0)
                fragmentLength = readLength;

            if(fragmentLength < readLength)
                return false;

            String alignedCigar = format("%dM", readLength);

            // only fully aligned reads
            if(!record.getCigarString().equals(alignedCigar))
                return false;

            // without too many mismatches
            Integer nmCount = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
            if(nmCount != null && nmCount.intValue() > 2)
                return false;

            if(record.getReadPairedFlag())
            {
                if(record.getSecondOfPairFlag())
                    return false;

                // ignore translocations and inversions
                if(!record.getMateReferenceName().equals(record.getReferenceName()))
                    return false;

                if(mateNegativeStrand(record) == record.getReadNegativeStrandFlag())
                    return false;

                if(record.isSecondaryOrSupplementary())
                    return false;

                String mateCigar = record.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

                if(mateCigar != null && !mateCigar.equals(alignedCigar))
                    return false;
            }

            return true;
        }

        private void addFragmentLength(final SAMRecord record)
        {
            int fragmentLength = getLengthBucket(abs(inferredInsertSize(record)));

            if(fragmentLength <= 0)
                return;

            int index = 0;
            boolean exists = false;
            while(index < mLengthFrequencies.size())
            {
                final LengthFrequency fragLengthCount = mLengthFrequencies.get(index);

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
                LengthFrequency newFragLengthCount = new LengthFrequency(fragmentLength, 1);
                mLengthFrequencies.add(index, newFragLengthCount);
            }
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
    }

    private void writeDistribution()
    {
        try
        {
            final String outputFileName = mConfig.formFilename(FRAGMENT_LENGTH_DIST);

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("FragmentLength\tCount");
            writer.newLine();

            for(LengthFrequency lengthFrequency : mLengthFrequencies)
            {
                // cap any fragmemt distribution entry at the observed read length to avoid the use of trimmed fragments impacting it
                if(mHasPairedReads &&  mMaxReadLength > 0 && lengthFrequency.Length < mMaxReadLength)
                    continue;

                writer.write(format("%d\t%d", lengthFrequency.Length, lengthFrequency.Frequency));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write fragment length file: {}", e.toString());
        }
    }

    public static FragmentLengthBounds loadFragmentLengthBounds(final String filename)
    {
        try
        {
            List<LengthFrequency> lengthFrequencies = Lists.newArrayList();

            List<String> lines = Files.readAllLines(Paths.get(filename));
            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM);

                lengthFrequencies.add(new LengthFrequency(Integer.parseInt(values[0]), Integer.parseInt(values[1])));
            }

            return calculateFragmentLengthBounds(lengthFrequencies);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to read fragment length file: {}", e.toString());
            return INVALID;
        }
    }
}
