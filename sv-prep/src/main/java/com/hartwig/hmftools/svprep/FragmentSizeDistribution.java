package com.hartwig.hmftools.svprep;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConstants.FRAG_LENGTH_DIST_MAX_LENGTH;
import static com.hartwig.hmftools.svprep.SvConstants.FRAG_LENGTH_DIST_MIN_QUAL;
import static com.hartwig.hmftools.svprep.SvConstants.FRAG_LENGTH_DIST_PERCENTILE;
import static com.hartwig.hmftools.svprep.SvConstants.FRAG_LENGTH_DIST_SAMPLE_SIZE;
import static com.hartwig.hmftools.svprep.WriteType.FRAGMENT_LENGTH_DIST;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FragmentSizeDistribution
{
    private final SvConfig mConfig;
    private final List<LengthFrequency> mLengthFrequencies;
    private int mMaxReadLength;
    
    public FragmentSizeDistribution(final SvConfig config)
    {
        mConfig = config;
        mLengthFrequencies = Lists.newArrayList();
        mMaxReadLength = 0;
    }
    
    public void run()
    {
        SV_LOGGER.info("calculating fragment size distribution");

        final List<ChromosomeTask> chrTasks = Lists.newArrayList();

        // for(HumanChromosome chromosome : HumanChromosome.values())
        for(int i = 1; i <= 10; ++i)
        {
            String chromosome = String.valueOf(i);
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome);

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            ChromosomeTask chrTask = new ChromosomeTask(chromosomeStr);
            chrTasks.add(chrTask);
        }

        boolean validExecution = TaskExecutor.executeTasks(chrTasks, mConfig.Threads);

        if(!validExecution)
            return;

        // merge results from all chromosomes
        for(final ChromosomeTask chrTask : chrTasks)
        {
            mergeDistributions(mLengthFrequencies, chrTask.lengthFrequencies());
            mMaxReadLength = max(mMaxReadLength, chrTask.maxReadLength());
        }

        SV_LOGGER.info("maxReadLength({})", mMaxReadLength);

        if(mLengthFrequencies.isEmpty())
        {
            SV_LOGGER.debug("no fragment lengths recorded");
            return;
        }

        int minLength = mLengthFrequencies.get(0).Length;
        int maxLength = mLengthFrequencies.get(mLengthFrequencies.size() - 1).Length;

        SV_LOGGER.info("fragment size distribution complete: min({}) max({})", minLength, maxLength);

        if(mConfig.WriteTypes.contains(FRAGMENT_LENGTH_DIST))
            writeDistribution();
    }

    public int maxReadLength() { return mMaxReadLength; }

    public int[] calculatePercentileLengths()
    {
        int[] percentileLengths = { 0, 0 };

        if(mLengthFrequencies.isEmpty())
            return percentileLengths;

        int totalFragments = mLengthFrequencies.stream().mapToInt(x -> x.Frequency).sum();
        int cumulativeTotal = 0;
        int requiredMinTotal = (int)floor(totalFragments * (1 - FRAG_LENGTH_DIST_PERCENTILE));
        int requiredMaxTotal = (int)floor(totalFragments * FRAG_LENGTH_DIST_PERCENTILE);

        for(LengthFrequency lengthData : mLengthFrequencies)
        {
            if(percentileLengths[0] == 0 && lengthData.Frequency + cumulativeTotal >= requiredMinTotal)
            {
                percentileLengths[0] = lengthData.Length;
            }

            if(lengthData.Frequency + cumulativeTotal >= requiredMaxTotal)
            {
                percentileLengths[1] = lengthData.Length;
                break;
            }
            else
            {
                cumulativeTotal += lengthData.Frequency;
            }
        }

        return percentileLengths;
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
                writer.write(String.format("%d\t%d", lengthFrequency.Length, lengthFrequency.Frequency));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write fragment length file: {}", e.toString());
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

        public ChromosomeTask(final String chromosome)
        {
            mChromosome = chromosome;
            mProcessedReads = 0;
            mMaxReadLength = 0;

            mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
            mBamSlicer = new BamSlicer(FRAG_LENGTH_DIST_MIN_QUAL, false, false, false);

            mLengthFrequencies = Lists.newArrayList();
        }

        public List<LengthFrequency> lengthFrequencies() { return mLengthFrequencies; }
        public int maxReadLength() { return mMaxReadLength; }

        @Override
        public Long call()
        {
            // slice a fixed region from each chromosome
            ChrBaseRegion region = !mConfig.SpecificChrRegions.Regions.isEmpty() ?
                mConfig.SpecificChrRegions.Regions.get(0) : new ChrBaseRegion(mChromosome, 1_000_000, 10_000_000);

            mBamSlicer.slice(mSamReader, region, this::processBamRead);

            return (long)0;
        }

        private void processBamRead(final SAMRecord record)
        {
            // cull invalid reads without waiting for the paired read
            if(!isCandidateRecord(record))
                return;

            mMaxReadLength = max(mMaxReadLength, record.getReadBases().length);

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
            boolean isPaired = record.getReadPairedFlag();

            if(isPaired)
            {
                if(record.getSecondOfPairFlag())
                    return false;

                int fragmentLength = abs(record.getInferredInsertSize());
                if(fragmentLength > FRAG_LENGTH_DIST_MAX_LENGTH)
                    return false;

                // ignore translocations and inversions
                if(!record.getMateReferenceName().equals(record.getReferenceName()))
                    return false;

                if(mateNegativeStrand(record) == record.getReadNegativeStrandFlag())
                    return false;

                if(record.isSecondaryOrSupplementary())
                    return false;
            }

            if(record.getCigar().containsOperator(CigarOperator.N) || record.getCigar().containsOperator(CigarOperator.S))
                return false;

            return true;
        }

        private void addFragmentLength(final SAMRecord record)
        {
            int fragmentLength = getLengthBucket(abs(record.getInferredInsertSize()));

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

    private static class LengthFrequency
    {
        public final int Length;
        public int Frequency;

        public LengthFrequency(final int length, final int frequency)
        {
            Length = length;
            Frequency = frequency;
        }

        public String toString() { return String.format("length(%d) freq(%d)", Length, Frequency); }
    }

}
