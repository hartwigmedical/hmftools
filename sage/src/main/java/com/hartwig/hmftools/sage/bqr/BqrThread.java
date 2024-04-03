package com.hartwig.hmftools.sage.bqr;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
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

    private final Queue<ChrBaseRegion> mRegions;
    private final BaseQualityResults mResults;
    private final int mRegionCount;

    private final Map<String, List<Integer>> mKnownVariantMap;

    private final BqrRegionReader mRegionCounter; // will be reused for each region

    public BqrThread(
            final SageConfig config, final IndexedFastaSequenceFile refGenome, final String bamFile,
            final Queue<ChrBaseRegion> regions, final BaseQualityResults results, final BqrRecordWriter recordWriter,
            final Map<String, List<Integer>> knownVariantMap)
    {
        mRefGenome = refGenome;
        mConfig = config;
        mRegions = regions;
        mRegionCount = regions.size();
        mResults = results;
        mKnownVariantMap = knownVariantMap;

        mBamReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.BamStringency)
                .referenceSource(new ReferenceSource(mRefGenome))
                .open(new File(bamFile));

        mRegionCounter = new BqrRegionReader(mConfig, mBamReader, mRefGenome, mResults, recordWriter);

        start();
    }

    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion partition = mRegions.remove();

                Set<Integer> knownPositions;
                if(!mKnownVariantMap.isEmpty())
                {
                    knownPositions = Sets.newHashSet();
                    List<Integer> snpPositions = mKnownVariantMap.get(partition.Chromosome);

                    if(snpPositions != null)
                        snpPositions.stream().filter(x -> partition.containsPosition(x)).forEach(x -> knownPositions.add(x));
                }
                else
                {
                    knownPositions = Collections.emptySet();
                }

                mRegionCounter.initialise(partition, knownPositions);

                mRegionCounter.run();

                int processed = mRegionCount - mRegions.size();

                if((processed % 100) == 0)
                {
                    SG_LOGGER.debug("base-qual regions processed({}) remaining({})", processed, mRegions.size());
                }

            }
            catch(NoSuchElementException e)
            {
                SG_LOGGER.trace("all tasks complete");
                break;
            }
        }
    }
}
