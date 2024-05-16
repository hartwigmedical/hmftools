package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.createBamSlicer;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.CombinedStats;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final String mChromosome;
    private final PrepConfig mConfig;
    private final SpanningReadCache mSpanningReadCache;
    private final ResultsWriter mWriter;
    private final CombinedStats mCombinedStats;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;
    private final Queue<ChrBaseRegion> mPartitions;

    private final int mPartitionCount;

    public PartitionThread(
            final String chromosome, final PrepConfig config, final Queue<ChrBaseRegion> partitions,
            final SpanningReadCache spanningReadCache, final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mChromosome = chromosome;
        mConfig = config;
        mSpanningReadCache = spanningReadCache;
        mWriter = writer;
        mCombinedStats = combinedStats;
        mPartitions = partitions;

        mPartitionCount = partitions.size();

        mSamReaders = Lists.newArrayList();

        for(String bamFile : mConfig.BamFiles)
        {
            mSamReaders.add(
                    SamReaderFactory.makeDefault().validationStringency(mConfig.BamStringency)
                            .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFile)));
        }

        mBamSlicer = createBamSlicer();

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion partition = mPartitions.remove();

                int processedCount = mPartitionCount - mPartitions.size();

                PartitionSlicer slicer = new PartitionSlicer(
                        0, partition, mConfig, mSamReaders, mBamSlicer, mSpanningReadCache, mWriter, mCombinedStats);

                slicer.run();

                if((processedCount % 10) == 0)
                {
                    SV_LOGGER.info("chromosome({}) processed {} partitions", mChromosome, mPartitions.size());
                }
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Throwable e)
            {
                SV_LOGGER.error("thread execution error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        try
        {
            for(SamReader reader : mSamReaders)
            {
                reader.close();
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
