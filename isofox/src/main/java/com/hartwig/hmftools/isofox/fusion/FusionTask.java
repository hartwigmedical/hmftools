package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FS_DOWNSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FS_UPSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionTask implements Callable
{
    private final int mTaskId;
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final Map<String,Map<Integer,BaseDepth>> mChrGeneDepthMap;
    private final List<FusionFragment> mAllFragments;
    private final Set<String> mReadIds;

    private int mNextFusionId;
    private final Map<String,List<FusionReadData>> mFusionCandidates; // keyed by the chromosome pair
    private final Map<String,List<FusionReadData>> mFusionsByGene; // keyed by the locationId
    private final Map<String,List<FusionFragment>> mDiscordantFragments; // keyed by the chromosome pair
    private final Map<String,List<FusionFragment>> mSingleGeneFragments; // keyed by chromosome, since single-sided

    private final FusionWriter mFusionWriter;

    private final PerformanceCounter[] mPerfCounters;

    private static final int PC_DISCOVERY = 0;
    private static final int PC_ANNOTATE = 1;
    private static final int PC_DEPTH = 2;

    public FusionTask(
            int taskId, final IsofoxConfig config, final EnsemblDataCache geneTransCache,
            final Map<String,Map<Integer,BaseDepth>> geneDepthMap, final List<FusionFragment> fragments, final FusionWriter fusionWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mGeneTransCache = geneTransCache;
        mChrGeneDepthMap = geneDepthMap;

        mAllFragments = Lists.newArrayList(fragments);

        mNextFusionId = 0;
        mFusionCandidates = Maps.newHashMap();
        mFusionsByGene = Maps.newHashMap();
        mDiscordantFragments = Maps.newHashMap();
        mSingleGeneFragments = Maps.newHashMap();
        mReadIds = Sets.newHashSet();

        mFusionWriter = fusionWriter;

        mPerfCounters = new PerformanceCounter[PC_DEPTH + 1];
        mPerfCounters[PC_DISCOVERY] = new PerformanceCounter("Discovery");
        mPerfCounters[PC_ANNOTATE] = new PerformanceCounter("Annotate");
        mPerfCounters[PC_DEPTH] = new PerformanceCounter("ReadDepth");
   }

    public final Map<String,List<FusionReadData>> getFusionCandidates() { return mFusionCandidates; }
    public final Map<String,List<FusionFragment>> getUnfusedFragments() { return mDiscordantFragments; }

    @Override
    public Long call()
    {
        formInitialFusions();

        calculateReadDepth();

        annotateFusions();

        writeData();

        if(mAllFragments.size() > 1000)
            ISF_LOGGER.info("{}: fusion task complete", mTaskId);

        return (long)1;
    }

    public final PerformanceCounter[] getPerfCounters() { return mPerfCounters; }

    private void formInitialFusions()
    {
        mPerfCounters[PC_DISCOVERY].start();;

        int junctioned = 0;
        for (FusionFragment fragment : mAllFragments)
        {
            mReadIds.add(fragment.readId());

            if (fragment.type() == MATCHED_JUNCTION)
            {
                ++junctioned;
                createOrUpdateFusion(fragment);
            }
            else if(fragment.isSingleGene())
            {
                cacheSingleGeneFragment(fragment);
            }
            else
            {
                cacheDiscordantFragment(fragment);
            }
        }

        ISF_LOGGER.info("{}: chimeric fragments({} unfused={} junc={}) fusions(loc={} gene={} total={})",
                mTaskId, mAllFragments.size(), mDiscordantFragments.values().stream().mapToInt(x -> x.size()).sum(), junctioned,
                mFusionCandidates.size(), mFusionsByGene.size(), mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum());

        mPerfCounters[PC_DISCOVERY].stop();
    }

    private void calculateReadDepth()
    {
        mPerfCounters[PC_DEPTH].start();

        ISF_LOGGER.info("{}: calculating junction depth", mTaskId);

        for (List<FusionReadData> fusions : mFusionCandidates.values())
        {
            for(FusionReadData fusionData : fusions)
            {
                for (int se = SE_START; se <= SE_END; ++se)
                {
                    final Map<Integer,BaseDepth> baseDepthMap = mChrGeneDepthMap.get(fusionData.chromosomes()[se]);

                    if(baseDepthMap == null)
                        continue;

                    final BaseDepth baseDepth = baseDepthMap != null ? baseDepthMap.get(fusionData.geneCollections()[se]) : null;

                    if(baseDepth == null)
                    {
                        ISF_LOGGER.error("fusion({}) gc({}) missing base depth", fusionData, fusionData.geneCollections()[se]);
                        continue;
                    }

                    int baseDepthCount = baseDepth.depthAtBase(fusionData.junctionPositions()[se]);
                    fusionData.getReadDepth()[se] += baseDepthCount;
                }
            }
        }

        mPerfCounters[PC_DEPTH].stop();
    }

    private void annotateFusions()
    {
        mPerfCounters[PC_ANNOTATE].start();

        ISF_LOGGER.debug("{}: marking related fusions", mTaskId);
        markRelatedFusions();

        // assign any discordant reads
        ISF_LOGGER.debug("{}: assigning unfused fragments", mTaskId);
        assignDiscordantFragments();
        assignRealignedFragments();

        mPerfCounters[PC_ANNOTATE].stop();
    }

    private void writeData()
    {
        // write results
        // ISF_LOGGER.debug("{}: writing results", mTaskId);
        mFusionWriter.writeFusionData(mFusionCandidates);
        mFusionWriter.writeUnfusedFragments(mDiscordantFragments);
    }

    private void createOrUpdateFusion(final FusionFragment fragment)
    {
        // scenarios:
        // 1. New fusion with correct splice-junction support - may or may not match a known transcript and exon
        // 2. Potential discordant or realigned fragment

        // fusions will be stored in a map keyed by their location pair (chromosome + geneCollectionId)
        List<FusionReadData> fusions = mFusionCandidates.get(fragment.locationPair());

        if(fusions == null)
        {
            fusions = Lists.newArrayList();
            mFusionCandidates.put(fragment.locationPair(), fusions);
        }

        for(final FusionReadData fusionData : fusions)
        {
            if(fusionData.junctionMatch(fragment))
            {
                fusionData.addFusionFragment(fragment);
                return;
            }
        }

        int fusionId = mTaskId * 10000 + mNextFusionId++;
        final FusionReadData fusionData = new FusionReadData(fusionId, fragment);
        setGeneData(fusionData);
        fusions.add(fusionData);

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final String chrGene = formLocation(fragment.chromosomes()[se], fragment.geneCollections()[se]);
            List<FusionReadData> fusionsByGene = mFusionsByGene.get(chrGene);

            if(fusionsByGene == null)
                mFusionsByGene.put(chrGene, Lists.newArrayList(fusionData));
            else
                fusionsByGene.add(fusionData);
        }
    }

    private void cacheDiscordantFragment(final FusionFragment fragment)
    {
        List<FusionFragment> fragments = mDiscordantFragments.get(fragment.locationPair());

        if(fragments == null)
            mDiscordantFragments.put(fragment.locationPair(), Lists.newArrayList(fragment));
        else
            fragments.add(fragment);
    }

    private void cacheSingleGeneFragment(final FusionFragment fragment)
    {
        final String locationId = formLocation(fragment.chromosomes()[SE_START], fragment.geneCollections()[SE_START]);
        List<FusionFragment> fragments = mSingleGeneFragments.get(locationId);

        if(fragments == null)
            mSingleGeneFragments.put(locationId, Lists.newArrayList(fragment));
        else
            fragments.add(fragment);
    }

    private void setGeneData(final FusionReadData fusionData)
    {
        // get the genes supporting the splice junction in the terms of an SV (ie lower chromosome and lower position first)
        final List<List<EnsemblGeneData>> genesByPosition = Lists.newArrayList(Lists.newArrayList(), Lists.newArrayList());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<String> spliceGeneIds = fusionData.getSampleFragment().getGeneIds(se);

            if(!spliceGeneIds.isEmpty())
            {
                genesByPosition.set(se, spliceGeneIds.stream()
                        .map(x -> mGeneTransCache.getGeneDataById(x)).collect(Collectors.toList()));
            }

            // purge any invalid transcript exons when the splice junction is known
            final List<TranscriptData> transDataList = Lists.newArrayList();
            spliceGeneIds.forEach(x -> transDataList.addAll(mGeneTransCache.getTranscripts(x)));

            for(FusionFragment fragment : fusionData.getAllFragments())
            {
                fragment.validateTranscriptExons(transDataList, se);
            }
        }

        // organise genes by strand based on the orientations around the splice junction
        // a positive orientation implies either an upstream +ve strand gene or a downstream -ve strand gene
        final byte[] sjOrientations = fusionData.junctionOrientations();

        boolean foundCandidates = false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int upstreamIndex = se;
            int downstreamIndex = switchIndex(se);

            // for the start being the upstream gene, if the orientation is +1 then require a +ve strand gene,
            // and so the end is the downstream gene and will set the orientation as required

            final List<EnsemblGeneData> upstreamGenes = genesByPosition.get(upstreamIndex).stream()
                    .filter(x -> x.Strand == sjOrientations[upstreamIndex]).collect(Collectors.toList());

            final List<EnsemblGeneData> downstreamGenes = genesByPosition.get(downstreamIndex).stream()
                    .filter(x -> x.Strand == -sjOrientations[downstreamIndex]).collect(Collectors.toList());

            if(!upstreamGenes.isEmpty() && !downstreamGenes.isEmpty())
            {
                if(foundCandidates)
                {
                    // both combinations have possible gene-pairings
                    ISF_LOGGER.debug("fusion({}) has multiple gene pairings by strand and orientation", fusionData.toString());
                    fusionData.setIncompleteData();
                    break;
                }

                foundCandidates = true;
                fusionData.setStreamData(upstreamGenes, downstreamGenes, se == SE_START);
            }
        }

        // mark donor-acceptor types whether strands are known or not
        final byte[] geneStrands = fusionData.getGeneStrands();

        for(FusionFragment fragment : fusionData.getAllFragments())
        {
            fragment.setJunctionTypes(mConfig.RefFastaSeqFile, geneStrands) ;
        }

        fusionData.cacheTranscriptData();
    }

    private void markRelatedFusions()
    {
        for( Map.Entry<String,List<FusionReadData>> entry : mFusionCandidates.entrySet())
        {
            // final String locationId = entry.getKey();
            final List<FusionReadData> fusions = entry.getValue();

            if(fusions.size() == 1)
                continue;

            // annotate similar fusions for post-run analysis
            for(int i = 0; i < fusions.size() - 1; ++i)
            {
                FusionReadData fusion1 = fusions.get(i);

                // firstly the special case of a spliced fusion matching its unspliced fusion
                boolean isSpliced = fusion1.isKnownSpliced();
                boolean isUnspliced = fusion1.isUnspliced();

                final List<TransExonRef> upRefs1 = fusion1.getTransExonRefsByStream(FS_UPSTREAM);
                final List<TransExonRef> downRefs1 = fusion1.getTransExonRefsByStream(FS_DOWNSTREAM);

                for(int j = i + 1; j < fusions.size() - 1; ++j)
                {
                    FusionReadData fusion2 = fusions.get(j);

                    if(isSpliced && fusion2.isUnspliced())
                    {
                        if(hasTranscriptExonMatch(upRefs1, fusion2.getTransExonRefsByStream(FS_UPSTREAM))
                                && hasTranscriptExonMatch(downRefs1, fusion2.getTransExonRefsByStream(FS_DOWNSTREAM), -1))
                        {
                            fusion1.addRelatedFusion(fusion2.id());
                            fusion2.addRelatedFusion(fusion1.id());
                            continue;
                        }
                    }
                    else if(isUnspliced && fusion2.isKnownSpliced())
                    {
                        if(hasTranscriptExonMatch(upRefs1, fusion2.getTransExonRefsByStream(FS_UPSTREAM))
                                && hasTranscriptExonMatch(fusion2.getTransExonRefsByStream(FS_DOWNSTREAM), downRefs1, -1))
                        {
                            fusion1.addRelatedFusion(fusion2.id());
                            fusion2.addRelatedFusion(fusion1.id());
                            continue;
                        }
                    }

                    boolean isSimilar = false;

                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        if (abs(fusion1.junctionPositions()[se] - fusion2.junctionPositions()[se]) <= POSITION_REALIGN_DISTANCE)
                        {
                            isSimilar = true;
                            break;
                        }
                    }

                    if(isSimilar)
                    {
                        fusion1.addRelatedFusion(fusion2.id());
                        fusion2.addRelatedFusion(fusion1.id());
                    }
                }
            }
        }
    }

    private static final int POSITION_REALIGN_DISTANCE = 20;

    private void assignDiscordantFragments()
    {
        for(Map.Entry<String,List<FusionFragment>> entry : mDiscordantFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionCandidates.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final Set<FusionFragment> allocatedFragments = Sets.newHashSet();

            for (FusionFragment fragment : fragments)
            {
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                for(FusionReadData fusionData : fusions)
                {
                    if(!fusionData.isValid())
                        continue;

                    if(fusionData.canAddUnfusedFragment(fragment, mConfig.MaxFragmentLength))
                    {
                        if(fusionData.isRelignedFragment(fragment))
                        {
                            fragment.setType(REALIGNED);
                        }
                        else
                        {
                            fragment.setType(DISCORDANT);
                        }

                        fusionData.addFusionFragment(fragment);
                        allocatedFragments.add(fragment);
                    }
                }
            }

            allocatedFragments.forEach(x -> fragments.remove(x));
        }
    }

    private void assignRealignedFragments()
    {
        for(Map.Entry<String,List<FusionFragment>> entry : mSingleGeneFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionsByGene.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final Set<FusionFragment> allocatedFragments = Sets.newHashSet();

            for (FusionFragment fragment : fragments)
            {
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                for(FusionReadData fusionData : fusions)
                {
                    if(!fusionData.isValid())
                        continue;

                    if(fusionData.isRelignedFragment(fragment))
                    {
                        fragment.setType(REALIGNED);
                        fusionData.addFusionFragment(fragment);
                        allocatedFragments.add(fragment);
                    }
                }
            }

            allocatedFragments.forEach(x -> fragments.remove(x));
        }
    }

    public void clearState()
    {
        mNextFusionId = 0;
        mFusionCandidates.clear();
        mDiscordantFragments.clear();
        mSingleGeneFragments.clear();
        mReadIds.clear();
    }

}

