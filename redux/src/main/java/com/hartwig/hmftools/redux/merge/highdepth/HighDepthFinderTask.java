package com.hartwig.hmftools.redux.merge.highdepth;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.redux.merge.highdepth.HighDepthFinderConfig.MD_LOGGER;
import static com.hartwig.hmftools.redux.merge.highdepth.HighDepthFinder.writeHighDepthRegions;

import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class HighDepthFinderTask implements Callable
{
    private final HighDepthFinderConfig mConfig;
    private final String mChromosome;
    private final BufferedWriter mWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final int[] mBaseDepth;

    private ChrBaseRegion mCurrentPartition;
    private HighDepthRegion mLastRegion;
    private final PerformanceCounter mPerfCounter;
    private int mRecordCounter;
    private int mHighDepthRegionCounter;

    public HighDepthFinderTask(final String chromosome, final HighDepthFinderConfig config, final BufferedWriter writer)
    {
        mConfig = config;
        mChromosome = chromosome;
        mWriter = writer;

        mBamSlicer = new BamSlicer(0, !config.DropDuplicates, true, false);

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));
        mBaseDepth = new int[mConfig.PartitionSize];
        mCurrentPartition = null;
        mLastRegion = null;

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

        List<ChrBaseRegion> chrRegions;
        if(mConfig.SpecificRegions.isEmpty())
        {
            chrRegions = Lists.newArrayList(new ChrBaseRegion(mChromosome, 1, chromosomeLength));
        }
        else
        {
            chrRegions = mConfig.SpecificRegions.stream().filter(x -> x.Chromosome.equals(mChromosome)).collect(Collectors.toList());
            ChrBaseRegion.checkMergeOverlaps(chrRegions, true);
        }

        List<ChrBaseRegion> partitions = Lists.newArrayList();
        for(ChrBaseRegion chrRegion : chrRegions)
        {
            for(int i = 0; ; i++)
            {
                int start = chrRegion.start() + i * mConfig.PartitionSize;
                int end = min(start + mConfig.PartitionSize - 1, chrRegion.end());
                partitions.add(new ChrBaseRegion(mChromosome, start, end));

                if(end >= chrRegion.end())
                    break;
            }
        }

        MD_LOGGER.info("chr({}) processing {} partitions", mChromosome, partitions.size());

        int processed = 0;
        for(ChrBaseRegion partition : partitions)
        {
            processPartition(partition);

            ++processed;

            if((processed % 100) == 0)
            {
                MD_LOGGER.info("chr({}) processed {} partitions", mChromosome, processed);
            }
        }

        if(mLastRegion != null)
        {
            writeHighDepthRegions(mWriter, List.of(mLastRegion));
            mHighDepthRegionCounter++;
        }

        MD_LOGGER.info("chr({}) processing complete, totalReads({}) highDepthRegions({})",
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
        ArrayList<HighDepthRegion> highDepthRegions = Lists.newArrayList();

        // Create regions left to right.
        HighDepthRegion currentRegion = null;
        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            int position = mCurrentPartition.start() + i;
            int baseDepth = mBaseDepth[i];

            if(baseDepth >= mConfig.InitHighDepthThreshold)
            {
                if(currentRegion == null)
                {
                    currentRegion = new HighDepthRegion(new ChrBaseRegion(mChromosome, position, position));
                    currentRegion.DepthMin = baseDepth;
                    currentRegion.DepthMax = baseDepth;
                    currentRegion.BaseVolume = baseDepth;
                    highDepthRegions.add(currentRegion);
                }
                else
                {
                    // extend the region
                    currentRegion.DepthMin = Math.min(currentRegion.DepthMin, baseDepth);
                    currentRegion.DepthMax = Math.max(currentRegion.DepthMax, baseDepth);
                    currentRegion.BaseVolume += baseDepth;
                    currentRegion.setEnd(position);
                }
            }
            else
            {
                if(currentRegion == null)
                    continue;

                if(baseDepth >= mConfig.FinalHighDepthThreshold)
                {
                    // extend the region
                    currentRegion.DepthMin = Math.min(currentRegion.DepthMin, baseDepth);
                    currentRegion.DepthMax = Math.max(currentRegion.DepthMax, baseDepth);
                    currentRegion.BaseVolume += baseDepth;
                    currentRegion.setEnd(position);
                    continue;
                }

                // end this region
                currentRegion = null;
            }
        }

        if(highDepthRegions.isEmpty())
        {
            if(mLastRegion != null)
            {
                writeHighDepthRegions(mWriter, List.of(mLastRegion));
                mHighDepthRegionCounter++;
                mLastRegion = null;
            }
            return;
        }

        // backfill regions
        for(HighDepthRegion region : highDepthRegions)
        {
            for(int position = region.start() - 1; position >= mCurrentPartition.start(); --position)
            {
                final int baseDepth = mBaseDepth[position - mCurrentPartition.start()];
                if(baseDepth < mConfig.FinalHighDepthThreshold)
                {
                    break;
                }

                // extend the region
                region.DepthMin = Math.min(region.DepthMin, baseDepth);
                region.DepthMax = Math.max(region.DepthMax, baseDepth);
                region.BaseVolume += baseDepth;
                region.setStart(position);
            }
        }

        // is the last region from the last partition and the first region from this partition contiguous?
        if(mLastRegion != null)
        {
            if(mLastRegion.end() == highDepthRegions.get(0).start() - 1)
            {
                highDepthRegions.get(0).merge(mLastRegion);
                mLastRegion = null;
            }
        }

        final HighDepthRegion nextLastRegion = highDepthRegions.get(highDepthRegions.size() - 1);
        highDepthRegions.remove(highDepthRegions.size() - 1);

        if(mLastRegion != null)
        {
            writeHighDepthRegions(mWriter, List.of(mLastRegion));
            mHighDepthRegionCounter++;
        }

        mLastRegion = nextLastRegion;

        if(!highDepthRegions.isEmpty())
        {
            writeHighDepthRegions(mWriter, highDepthRegions);
            mHighDepthRegionCounter += highDepthRegions.size();
        }
    }
}
