package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.Callable;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;

public class SampleGermlineGeneTask implements Callable<Void>
{
    private final int mTaskId;
    private final String mPurpleDir;
    private final List<String> mSampleIds;
    private final EventCounts mEventCounts = new EventCounts(new TreeMap<>());

    public SampleGermlineGeneTask(int taskId, final String purpleDir)
    {
        mTaskId = taskId;
        mPurpleDir = purpleDir;
        mSampleIds = Lists.newArrayList();
    }

    public List<String> getSampleIds()
    {
        return mSampleIds;
    }

    public EventCounts getEventCounts()
    {
        return mEventCounts;
    }

    @Override
    public Void call()
    {
        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);
            try
            {
                mEventCounts.mergeAdd(new EventCounts(processSampleUsingAmpDelData(sampleId)));
            }
            catch(Exception e)
            {
                PPL_LOGGER.error("could not process sample({}): {}", sampleId, e.toString());
            }

            if(i > 0 && (i % 100) == 0)
            {
                PPL_LOGGER.debug("{}: processed {} samples", mTaskId, i);
            }
        }

        PPL_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());

        return null;
    }

    private SortedMap<HumanChromosome, ChromosomeRegionCounts> processSampleUsingAmpDelData(String sampleId) throws IOException
    {
        RegionEventExtractor extractor = new RegionEventExtractor(GermlineAmpDel.generateFilename(mPurpleDir, sampleId));
        final ListMultimap<HumanChromosome, RegionGeneEvents> events = extractor.events();
        SortedMap<HumanChromosome, ChromosomeRegionCounts> result = new TreeMap<>();
        for(HumanChromosome chromosome : events.keySet())
        {
            List<RegionGeneEvents> regionGeneEvents = events.get(chromosome);
            result.put(chromosome, new ChromosomeRegionCounts(chromosome, regionGeneEvents));
        }
        return result;
    }
}

