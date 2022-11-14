package com.hartwig.hmftools.svprep.reads;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.svprep.CombinedStats;
import com.hartwig.hmftools.svprep.ExistingJunctionCache;
import com.hartwig.hmftools.svprep.ResultsWriter;
import com.hartwig.hmftools.svprep.SpanningReadCache;
import com.hartwig.hmftools.svprep.SvConfig;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final String mChromosome;
    private final SvConfig mConfig;
    private final SpanningReadCache mSpanningReadCache;
    private final ResultsWriter mWriter;
    private final CombinedStats mCombinedStats;
    private final ExistingJunctionCache mExistingJunctionCache;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final Queue<PartitionTask> mPartitions;

    public PartitionThread(
            final String chromosome, final SvConfig config, final Queue<PartitionTask> partitions,
            final SpanningReadCache spanningReadCache, final ExistingJunctionCache existingJunctionCache,
            final ResultsWriter writer, final CombinedStats combinedStats)
    {
        mChromosome = chromosome;
        mConfig = config;
        mSpanningReadCache = spanningReadCache;
        mExistingJunctionCache = existingJunctionCache;
        mWriter = writer;
        mCombinedStats = combinedStats;
        mPartitions = partitions;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mBamSlicer.setKeepUnmapped();
        mBamSlicer.setKeepHardClippedSecondaries();

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mPartitions.remove();

                PartitionSlicer slicer = new PartitionSlicer(
                        partition.TaskId, partition.Region, mConfig, mSamReader, mBamSlicer,
                        mSpanningReadCache, mExistingJunctionCache, mWriter, mCombinedStats);

                boolean logAndGc = partition.TaskId > 0 && (partition.TaskId % 10) == 0;

                if(logAndGc)
                {
                    SV_LOGGER.debug("chromosome({}) processing partition({}), remaining({})",
                            mChromosome, partition.TaskId, mPartitions.size());
                }

                slicer.run();

                if(logAndGc)
                    System.gc();
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                SV_LOGGER.error("thread execution error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
