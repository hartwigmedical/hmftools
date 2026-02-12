package com.hartwig.hmftools.redux.ms_sites;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ms_sites.MsFinderConfig.MIN_MAPPABILITY;
import static com.hartwig.hmftools.redux.ms_sites.MsFinderConfig.MIN_TARGET_SITE_COUNT;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.redux.jitter.MicrosatelliteSite;

public class SiteDownsampler
{
    private final MsFinderConfig mConfig;
    private final Multimap<UnitRepeatKey,MicrosatelliteSite> mAllMsSites;
    private final Map<String,List<BaseRegion>> mExomeRegions;

    // this needs to be thread safe
    private final Multimap<UnitRepeatKey, MicrosatelliteSite> mDownSampledMsSites;

    public SiteDownsampler(
            final MsFinderConfig config, final Multimap<UnitRepeatKey,MicrosatelliteSite> allMicrosatelliteSites,
            final Map<String,List<BaseRegion>> exomeRegions)
    {
        mConfig = config;
        mAllMsSites = allMicrosatelliteSites;
        mExomeRegions = exomeRegions;

        mDownSampledMsSites = Multimaps.synchronizedListMultimap(ArrayListMultimap.create());
    }

    public Multimap<UnitRepeatKey,MicrosatelliteSite> downSampledMsSites() { return mDownSampledMsSites; }

    private static int downsampleCount(UnitRepeatKey unitRepeatKey)
    {
        int rawValue = MsFinderConfig.DOWNSAMPLE_FACTOR / (int) Math.pow(2, (unitRepeatKey.NumRepeats + unitRepeatKey.Key.length * 2));
        return Math.max(MIN_TARGET_SITE_COUNT, rawValue);
    }

    public void downsampleSites() throws ExecutionException, InterruptedException
    {
        RD_LOGGER.info("down-sampling sites");

        // filter the microsatellites such that each type of (unit, length) is approximately the target count
        int numDigits = Integer.toString(mConfig.Threads - 1).length();

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("thread-%0" + numDigits + "d").build();
        ExecutorService executorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        // first we need to decide how much to downsample by. We do so by working out how many sites need to be kept
        for(UnitRepeatKey unitRepeatKey : mAllMsSites.keySet())
        {
            int downsampleCountPerType = downsampleCount(unitRepeatKey);

            RD_LOGGER.trace("[{}] filtering microsatellite sites, target count per type = {}",
                    unitRepeatKey, downsampleCountPerType);

            Runnable task = () -> downSampleMsType(downsampleCountPerType, unitRepeatKey);
            futures.add(CompletableFuture.runAsync(task, executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        RD_LOGGER.debug("filtered {} microsatellite sites down to {}", mAllMsSites.size(), mDownSampledMsSites.size());
    }

    private static List<MicrosatelliteSite> downsampleList(List<MicrosatelliteSite> originalList, int numItemsToKeep, Random seed)
    {
        if(numItemsToKeep >= originalList.size())
            return originalList;
        else if(numItemsToKeep <= 0)
            return Collections.emptyList();

        Collections.shuffle(originalList, seed);
        return originalList.subList(0, numItemsToKeep);
    }

    private void downSampleMsType(final int targetCountPerType, final UnitRepeatKey unitRepeatKey)
    {
        Collection<MicrosatelliteSite> allList = mAllMsSites.get(unitRepeatKey);

        RD_LOGGER.trace("[{}] filtering {} microsatellite sites", unitRepeatKey, allList.size());

        List<MicrosatelliteSite> sitesOutsideBedRegions = new ArrayList<>();
        List<MicrosatelliteSite> sitesInsideBedRegions = new ArrayList<>();

        // we filter to sites satisfying min mappability constraints
        // we then downsample to our target downsample count, prioritising sites inside the bed regions + 950 base buffer
        // if sites inside our bed regions already exceed target downsample count, we downsample these and take nothing from off-target

        int bedRegionExpansion = MsFinderConfig.BED_REGION_EXPANSION;
        double mappabilityCutoff = MIN_MAPPABILITY;

        for(MicrosatelliteSite msSite : allList)
        {
            if(msSite == null)
                continue;

            if(msSite.mappability() < mappabilityCutoff)
                continue;

            boolean isInBed = false;

            for(BaseRegion region : mExomeRegions.get(msSite.chromosome()))
            {
                if(positionsOverlap(
                        msSite.Region.start(), msSite.Region.end(),
                        region.start() - bedRegionExpansion, region.end() + bedRegionExpansion))
                {
                    isInBed = true;
                    break;
                }
            }
            if(isInBed)
            {
                sitesInsideBedRegions.add(msSite);
            }
            else
            {
                sitesOutsideBedRegions.add(msSite);
            }
        }

        // now decide how many to filter out from the rest
        int numSitesInBed = sitesInsideBedRegions.size();
        Random randomSeed = new Random(0);
        Collection<MicrosatelliteSite> filteredList = mDownSampledMsSites.get(unitRepeatKey);

        int sitesOutsideBedToKeep = targetCountPerType - numSitesInBed;
        List<MicrosatelliteSite> downsampledSitesInsideBedRegions = downsampleList(sitesInsideBedRegions, targetCountPerType, randomSeed);
        List<MicrosatelliteSite> downsampledSitesOutsideBedRegions =
                downsampleList(sitesOutsideBedRegions, sitesOutsideBedToKeep, randomSeed);
        filteredList.addAll(downsampledSitesInsideBedRegions);
        filteredList.addAll(downsampledSitesOutsideBedRegions);

        if(filteredList.size() < allList.size())
        {
            RD_LOGGER.debug("[{}] sites({} filtered={}) exome(orig={} filtered={}) outside-exome({} filtered={})",
                    unitRepeatKey, allList.size(), filteredList.size(), numSitesInBed, downsampledSitesInsideBedRegions.size(),
                    sitesOutsideBedRegions.size(), downsampledSitesOutsideBedRegions.size());
        }
    }
}
