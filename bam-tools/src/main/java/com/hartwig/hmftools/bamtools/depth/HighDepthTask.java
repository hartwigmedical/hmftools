package com.hartwig.hmftools.bamtools.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.depth.HighDepthFinder.writeHighDepthRegions;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;

import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class HighDepthTask implements Callable<Void>
{
    private final FinderConfig mConfig;
    private final String mChromosome;
    private final BufferedWriter mWriter;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final int[] mBaseDepth;

    private ChrBaseRegion mCurrentPartition;
    private int mRecordCounter;
    private int mHighDepthRegionCounter;

    public HighDepthTask(final String chromosome, final FinderConfig config, final BufferedWriter writer)
    {
        mConfig = config;
        mChromosome = chromosome;
        mWriter = writer;

        mBamSlicer = new BamSlicer(0, false, true, false);

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));
        mBaseDepth = new int[mConfig.PartitionSize];
        mCurrentPartition = null;

        mRecordCounter = 0;
        mHighDepthRegionCounter = 0;
    }

    @Override
    public Void call()
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

        BT_LOGGER.info("chr({}) processing {} partitions", mChromosome, partitions.size());

        int processed = 0;
        for(ChrBaseRegion partition : partitions)
        {
            processPartition(partition);

            ++processed;

            if((processed % 100) == 0)
            {
                BT_LOGGER.info("chr({}) processed {} partitions", mChromosome, processed);
            }
        }

        BT_LOGGER.info("chr({}) processing complete, totalReads({}) highDepthRegions({})",
                mChromosome, mRecordCounter, mHighDepthRegionCounter);

        return null;
    }

    private void processPartition(final ChrBaseRegion partition)
    {
        for(int i = 0; i < mBaseDepth.length; ++i)
        {
            mBaseDepth[i] = 0;
        }

        mCurrentPartition = partition;

        mBamSlicer.slice(mSamReader, mCurrentPartition, this::processSamRecord);

        findHighDepthRegions();
    }

    private void processSamRecord(final SAMRecord record)
    {
        if(record.getMappingQuality() < mConfig.MinMapQual)
            return;

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

        long regionDepthTotal = 0;
        long depthBelowThreshold = 0;

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
                    regionDepthTotal = baseDepth;

                    highDepthRegions.add(currentRegion);
                }
                else
                {
                    // extend the region
                    currentRegion.setEnd(position);
                    currentRegion.DepthMax = max(currentRegion.DepthMax, baseDepth);
                    regionDepthTotal += baseDepth;

                    regionDepthTotal += depthBelowThreshold;
                    depthBelowThreshold = 0;
                }
            }
            else
            {
                if(currentRegion == null)
                    continue;

                depthBelowThreshold += baseDepth;

                if(position - currentRegion.end() < mConfig.MaxRegionGap) // continue checking but don't extend the region
                    continue;

                // end this region
                currentRegion.DepthAvg = (int)Math.round(regionDepthTotal / (double)currentRegion.baseLength());
                currentRegion = null;
                regionDepthTotal = 0;
                depthBelowThreshold = 0;
            }
        }

        if(currentRegion != null)
            currentRegion.DepthAvg = (int)Math.round(regionDepthTotal / (double)currentRegion.baseLength());

        if(!highDepthRegions.isEmpty())
        {
            writeHighDepthRegions(mWriter, highDepthRegions);
            mHighDepthRegionCounter += highDepthRegions.size();
        }
    }
}
