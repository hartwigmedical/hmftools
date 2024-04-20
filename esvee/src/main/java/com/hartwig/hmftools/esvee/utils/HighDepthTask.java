package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.utils.HighDepthConfig.HIGH_DEPTH_REGION_MAX_GAP;
import static com.hartwig.hmftools.esvee.utils.HighDepthFinder.writeHighDepthRegions;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class HighDepthTask implements Callable
{
    private final HighDepthConfig mConfig;
    private final String mChromosome;
    private final BufferedWriter mWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final int[] mBaseDepth;

    private ChrBaseRegion mCurrentPartition;
    private final PerformanceCounter mPerfCounter;
    private int mRecordCounter;
    private int mHighDepthRegionCounter;

    public HighDepthTask(final String chromosome, final HighDepthConfig config, final BufferedWriter writer)
    {
        mConfig = config;
        mChromosome = chromosome;
        mWriter = writer;

        mBamSlicer = new BamSlicer(0, false, true, false);

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));
        mBaseDepth = new int[mConfig.PartitionSize];
        mCurrentPartition = null;

        mPerfCounter = new PerformanceCounter("Slice");
        mRecordCounter = 0;
        mHighDepthRegionCounter = 0;
    }

    @Override
    public Long call()
    {
        RefGenomeCoordinates refGenomeCoords = mConfig.RefGenVersion == V37 ?
                RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        int chromosomeLength = refGenomeCoords.length(stripChrPrefix(mChromosome));

        List<ChrBaseRegion> partitions = Lists.newArrayList();

        if(!mConfig.SpecificRegions.isEmpty())
        {
            mConfig.SpecificRegions.stream().filter(x -> x.Chromosome.equals(mChromosome)).forEach(x -> partitions.add(x));
        }
        else
        {
            for(int i = 0; ; i++)
            {
                int start = 1 + i * mConfig.PartitionSize;
                int end = min(start + mConfig.PartitionSize - 1, chromosomeLength);
                partitions.add(new ChrBaseRegion(mChromosome, start, end));

                if(end >= chromosomeLength)
                    break;
            }
        }

        SV_LOGGER.info("chr({}) processing {} partitions", mChromosome, partitions.size());

        int processed = 0;
        for(ChrBaseRegion partition : partitions)
        {
            processPartition(partition);

            ++processed;

            if((processed % 100) == 0)
            {
                SV_LOGGER.info("chr({}) processed {} partitions", mChromosome, processed);
            }
        }

        SV_LOGGER.info("chr({}) processing complete, totalReads({}) highDepthRegions({})",
                mChromosome, mRecordCounter, mHighDepthRegionCounter);
        mPerfCounter.logStats();

        return (long)0;
    }

    private void processPartition(final ChrBaseRegion partition)
    {
        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            mBaseDepth[i] = 0;
        }

        mCurrentPartition = partition;

        mPerfCounter.start();

        mBamSlicer.slice(mSamReader, mCurrentPartition, this::processSamRecord);

        findHighDepthRegions();

        mPerfCounter.stop();
    }

    private void processSamRecord(final SAMRecord record)
    {
        ++mRecordCounter;

        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();
        int baseStart = max(readStart - mCurrentPartition.start(), 0);
        int baseEnd = min(readEnd - mCurrentPartition.start(), mBaseDepth.length - 1);

        for(int i = baseStart; i <= baseEnd; ++i)
        {
            ++mBaseDepth[i];
        }
    }

    private void findHighDepthRegions()
    {
        List<HighDepthRegion> highDepthRegions = Lists.newArrayList();

        HighDepthRegion currentRegion = null;

        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            int position = mCurrentPartition.start() + i;
            int baseDepth = mBaseDepth[i];

            if(baseDepth >= mConfig.HighDepthThreshold)
            {
                if(currentRegion == null)
                {
                    currentRegion = new HighDepthRegion(new ChrBaseRegion(mChromosome, position, position));
                    currentRegion.DepthMin = baseDepth;
                    currentRegion.DepthMax = baseDepth;
                    highDepthRegions.add(currentRegion);
                }
                else
                {
                    // extend the region
                    currentRegion.setEnd(position);
                    currentRegion.DepthMax = max(currentRegion.DepthMax, baseDepth);
                }
            }
            else
            {
                if(currentRegion == null)
                    continue;

                if(position - currentRegion.end() < HIGH_DEPTH_REGION_MAX_GAP) // continue checking but don't extend the region
                    continue;

                // end this region
                currentRegion = null;
            }
        }

        if(!highDepthRegions.isEmpty())
        {
            writeHighDepthRegions(mWriter, highDepthRegions);
            mHighDepthRegionCounter += highDepthRegions.size();
        }
    }
}
