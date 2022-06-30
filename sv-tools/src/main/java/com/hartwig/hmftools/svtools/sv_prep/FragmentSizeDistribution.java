package com.hartwig.hmftools.svtools.sv_prep;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.svtools.sv_prep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.FRAG_LENGTH_DIST_SAMPLE_SIZE;
import static com.hartwig.hmftools.svtools.sv_prep.SvConstants.MAX_FRAGMENT_LENGTH;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FragmentSizeDistribution
{
    private final SvConfig mConfig;
    private final List<LengthFrequency> mLengthFrequencies;
    
    public FragmentSizeDistribution(final SvConfig config)
    {
        mConfig = config;
        mLengthFrequencies = Lists.newArrayList();
    }
    
    public void run()
    {
        SV_LOGGER.info("calculating fragment size distribution");

        final List<ChromosomeTask> chrTasks = Lists.newArrayList();

        // for(HumanChromosome chromosome : HumanChromosome.values())
        for(int i = 1; i <= 10; ++i)
        {
            String chromosome = String.valueOf(i);
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosomeStr))
                continue;

            ChromosomeTask chrTask = new ChromosomeTask(chromosomeStr);
            chrTasks.add(chrTask);
        }

        final List<Callable> callableList = chrTasks.stream().collect(Collectors.toList());
        boolean validExecution = TaskExecutor.executeTasks(callableList, mConfig.Threads);

        if(!validExecution)
            return;

        // merge results from all chromosomes
        for(final ChromosomeTask chrTask : chrTasks)
        {
            mergeDistributions(mLengthFrequencies, chrTask.lengthFrequencies());
        }

        final int[] lengthRange = calcFragmentLengthRange();

        SV_LOGGER.info("fragment size distribution complete: min({}) max({})", lengthRange[0], lengthRange[1]);
    }

    public final int[] calcFragmentLengthRange()
    {
        int minLength = mLengthFrequencies.get(0).Length;
        int maxLength = mLengthFrequencies.get(mLengthFrequencies.size() - 1).Length;
        return new int[] { minLength, maxLength };
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

    private class ChromosomeTask implements Callable
    {
        private final String mChromosome;
        private int mProcessedReads;

        private final BamSlicer mBamSlicer;
        private final SamReader mSamReader;
        private final List<LengthFrequency> mLengthFrequencies;

        public ChromosomeTask(final String chromosome)
        {
            mChromosome = chromosome;
            mProcessedReads = 0;

            mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
            mBamSlicer = new BamSlicer(60, false, false, false);

            mLengthFrequencies = Lists.newArrayList();
        }

        public List<LengthFrequency> lengthFrequencies() { return mLengthFrequencies; }

        @Override
        public Long call()
        {
            ChrBaseRegion region = new ChrBaseRegion(mChromosome, 1000_000, 20_000_000);
            mBamSlicer.slice(mSamReader, Lists.newArrayList(region), this::processBamRead);

            return (long)0;
        }

        private void processBamRead(@NotNull final SAMRecord read)
        {
            // cull invalid reads without waiting for the paired read
            if(!isCandidateRecord(read))
                return;

            ++mProcessedReads;

            if(mProcessedReads >= FRAG_LENGTH_DIST_SAMPLE_SIZE)
            {
                SV_LOGGER.debug("chr({}) reached max fragment count", mChromosome);
                mBamSlicer.haltProcessing();
                return;
            }

            addFragmentLength(read);
        }

        private boolean isCandidateRecord(final SAMRecord record)
        {
            if(!record.getFirstOfPairFlag())
                return false;

            int fragmentLength = abs(record.getInferredInsertSize());
            if(fragmentLength > MAX_FRAGMENT_LENGTH)
                return false;

            // ignore translocations and inversions
            if(!record.getMateReferenceName().equals(record.getReferenceName()) || record.getMateNegativeStrandFlag() == record.getReadNegativeStrandFlag())
                return false;

            if(record.isSecondaryOrSupplementary())
                return false;

            // ignore split and soft-clipped reads above the read length
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

    private class LengthFrequency
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
