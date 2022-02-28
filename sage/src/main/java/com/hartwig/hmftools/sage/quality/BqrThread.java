package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.PartitionTask;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BqrThread extends Thread
{
    private final IndexedFastaSequenceFile mRefGenome;
    private final SageConfig mConfig;
    private final SamReader mBamReader;

    private final Queue<PartitionTask> mRegions;
    private final BaseQualityResults mResults;

    public BqrThread(
            final SageConfig config, final IndexedFastaSequenceFile refGenome, final String bamFile,
            final Queue<PartitionTask> regions, final BaseQualityResults results)
    {
        mRefGenome = refGenome;
        mConfig = config;
        mRegions = regions;
        mResults = results;

        mBamReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefGenome))
                .open(new File(bamFile));

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                PartitionTask partition = mRegions.remove();

                BaseQualityRegionCounter regionCounter = new BaseQualityRegionCounter(
                        mConfig, mBamReader, mRefGenome, partition.Partition, mResults);

                if(partition.TaskId > 0 && (partition.TaskId % 100) == 0)
                {
                    SG_LOGGER.debug("base-qual regions assigned({}) remaining({})",
                            partition.TaskId, mRegions.size());
                }

                regionCounter.run();
            }
            catch(NoSuchElementException e)
            {
                SG_LOGGER.trace("all tasks complete");
                break;
            }
        }
    }

}
