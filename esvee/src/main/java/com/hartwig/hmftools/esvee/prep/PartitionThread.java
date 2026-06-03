package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.createBamSlicer;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.common.saga.SagaMatcherFactory;
import com.hartwig.hmftools.esvee.prep.types.CombinedStats;

import org.jetbrains.annotations.Nullable;

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

    @Nullable
    private final SagaMatcherFactory mSagaMatcherFactory;

    private final int mPartitionCount;

    public PartitionThread(
            final String chromosome, final PrepConfig config, final Queue<ChrBaseRegion> partitions,
            final SpanningReadCache spanningReadCache, @Nullable final SagaMatcherFactory sagaMatcherFactory, final ResultsWriter writer,
            final CombinedStats combinedStats)
    {
        mChromosome = chromosome;
        mConfig = config;
        mSpanningReadCache = spanningReadCache;
        mWriter = writer;
        mCombinedStats = combinedStats;
        mPartitions = partitions;
        mSagaMatcherFactory = sagaMatcherFactory;

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
            ChrBaseRegion partition;
            try
            {
                partition = mPartitions.remove();
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }

            try
            {
                int processedCount = mPartitionCount - mPartitions.size();

                PartitionSlicer slicer = new PartitionSlicer(
                        partition, mConfig, mSamReaders, mBamSlicer, mSpanningReadCache, mSagaMatcherFactory, mWriter, mCombinedStats);

                slicer.run();

                if((processedCount % 10) == 0)
                {
                    SV_LOGGER.debug("chromosome({}) processed {} partitions", mChromosome, mPartitions.size());
                }
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
