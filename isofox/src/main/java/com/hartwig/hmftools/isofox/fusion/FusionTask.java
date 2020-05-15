package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_NOVEL_SJ_DISTANCE;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.FS_DOWNSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.FS_UPSTREAM;

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
    private final Map<String,List<FusionFragment>> mRealignCandidateFragments; // keyed by chromosome, since single-sided

    private final FusionWriter mFusionWriter;

    private final PerformanceCounter[] mPerfCounters;

    private static final int PC_DISCOVERY = 0;
    private static final int PC_ANNOTATE = 1;

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
        mRealignCandidateFragments = Maps.newHashMap();
        mReadIds = Sets.newHashSet();

        mFusionWriter = fusionWriter;

        mPerfCounters = new PerformanceCounter[PC_ANNOTATE + 1];
        mPerfCounters[PC_DISCOVERY] = new PerformanceCounter("Discovery");
        mPerfCounters[PC_ANNOTATE] = new PerformanceCounter("Annotate");
   }

    public final Map<String,List<FusionReadData>> getFusionCandidates() { return mFusionCandidates; }
    public final Map<String,List<FusionFragment>> getUnfusedFragments() { return mDiscordantFragments; }

    @Override
    public Long call()
    {
        ISF_LOGGER.info("{}: processing {} chimeric fragments", mTaskId, mAllFragments.size());

        mPerfCounters[PC_DISCOVERY].start();;

        formInitialFusions();

        // assign any discordant reads
        ISF_LOGGER.debug("{}: assigning unfused fragments", mTaskId);
        assignDiscordantFragments();
        createLocalFusions();
        assignRealignedFragments();

        mPerfCounters[PC_DISCOVERY].stop();

        mPerfCounters[PC_ANNOTATE].start();

        annotateFusions();

        mPerfCounters[PC_ANNOTATE].stop();

        writeData();

        if(mAllFragments.size() > 1000)
            ISF_LOGGER.info("{}: fusion task complete", mTaskId);

        return (long)1;
    }

    public final PerformanceCounter[] getPerfCounters() { return mPerfCounters; }

    private void formInitialFusions()
    {
        int junctioned = 0;
        for (FusionFragment fragment : mAllFragments)
        {
            mReadIds.add(fragment.readId());

            if (fragment.type() == MATCHED_JUNCTION)
            {
                fragment.setJunctionBases(mConfig.RefFastaSeqFile);

                if(fragment.hasSuppAlignment())
                {
                    ++junctioned;
                    createOrUpdateFusion(fragment);
                }
                else
                {
                    // may be a local fusion candidate, or just supporting another fusion but without supp alignment data
                    cacheDiscordantFragment(fragment);
                }
            }
            else if(fragment.type() == REALIGN_CANDIDATE)
            {
                cacheRealignCandidateFragment(fragment);
            }
            else
            {
                cacheDiscordantFragment(fragment);
            }
        }

        ISF_LOGGER.info("{}: chimeric fragments({} disc={} candRealgn={} junc={}) fusions(loc={} gene={} total={})",
                mTaskId, mAllFragments.size(), mDiscordantFragments.values().stream().mapToInt(x -> x.size()).sum(),
                mRealignCandidateFragments.values().stream().mapToInt(x -> x.size()).sum(), junctioned,
                mFusionCandidates.size(), mFusionsByGene.size(), mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum());
    }

    private void annotateFusions()
    {
        ISF_LOGGER.debug("{}: annotating fusions", mTaskId);

        markRelatedFusions();

        for (List<FusionReadData> fusions : mFusionCandidates.values())
        {
            for (FusionReadData fusionData : fusions)
            {
                fusionData.calcJunctionDepth(mChrGeneDepthMap);
                fusionData.calcMaxSplitMappedLength();
            }
        }
    }

    private void writeData()
    {
        // write results
        mFusionWriter.writeFusionData(mFusionCandidates);
        mFusionWriter.writeUnfusedFragments(mDiscordantFragments);
        mFusionWriter.writeUnfusedFragments(mRealignCandidateFragments);
    }

    private FusionReadData createOrUpdateFusion(final FusionFragment fragment)
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
        else
        {
            for (final FusionReadData fusionData : fusions)
            {
                if (fusionData.junctionMatch(fragment))
                {
                    fusionData.addFusionFragment(fragment);

                    // mark donor-acceptor types whether strands are known or not
                    fragment.junctionTypes()[SE_START] = fusionData.getInitialFragment().junctionTypes()[SE_START];
                    fragment.junctionTypes()[SE_END] = fusionData.getInitialFragment().junctionTypes()[SE_END];
                    return fusionData;
                }
            }
        }

        int fusionId = mTaskId * 1000000 + mNextFusionId++; // keep unique across threaded tasks
        final FusionReadData fusionData = new FusionReadData(fusionId, fragment);

        setGeneData(fusionData);

        fusions.add(fusionData);

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && fragment.locationIds()[SE_START].equals(fragment.locationIds()[SE_END]))
                break;

            List<FusionReadData> fusionsByGene = mFusionsByGene.get(fragment.locationIds()[se]);

            if(fusionsByGene == null)
                mFusionsByGene.put(fragment.locationIds()[se], Lists.newArrayList(fusionData));
            else
                fusionsByGene.add(fusionData);
        }

        return fusionData;
    }

    private void cacheDiscordantFragment(final FusionFragment fragment)
    {
        List<FusionFragment> fragments = mDiscordantFragments.get(fragment.locationPair());

        if(fragments == null)
            mDiscordantFragments.put(fragment.locationPair(), Lists.newArrayList(fragment));
        else
            fragments.add(fragment);
    }

    private void cacheRealignCandidateFragment(final FusionFragment fragment)
    {
        List<FusionFragment> fragments = mRealignCandidateFragments.get(fragment.locationIds()[SE_START]);

        if(fragments == null)
            mRealignCandidateFragments.put(fragment.locationIds()[SE_START], Lists.newArrayList(fragment));
        else
            fragments.add(fragment);
    }

    private void setGeneData(final FusionReadData fusionData)
    {
        // get the genes supporting the splice junction in the terms of an SV (ie lower chromosome and lower position first)
        final List<List<EnsemblGeneData>> genesByPosition = Lists.newArrayList(Lists.newArrayList(), Lists.newArrayList());

        FusionFragment initialFragment = fusionData.getInitialFragment();

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<String> spliceGeneIds = fusionData.getInitialFragment().getGeneIds(se);

            if(!spliceGeneIds.isEmpty())
            {
                genesByPosition.set(se, spliceGeneIds.stream()
                        .map(x -> mGeneTransCache.getGeneDataById(x)).collect(Collectors.toList()));
            }

            // purge any invalid transcript exons when the splice junction is known
            final List<TranscriptData> transDataList = Lists.newArrayList();
            spliceGeneIds.forEach(x -> transDataList.addAll(mGeneTransCache.getTranscripts(x)));

            // purge any invalid transcript-exons and mark the junction as known if applicable
            initialFragment.validateTranscriptExons(transDataList, se);
        }

        // organise genes by strand based on the orientations around the splice junction
        // a positive orientation implies either an upstream +ve strand gene or a downstream -ve strand gene
        final byte[] sjOrientations = fusionData.junctionOrientations();

        boolean foundBothStreams = false;
        boolean foundOneStream = false;

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
                if(foundBothStreams)
                {
                    // both combinations have possible gene-pairings
                    ISF_LOGGER.debug("fusion({}) has multiple gene pairings by strand and orientation", fusionData.toString());
                    fusionData.setIncompleteData();
                    break;
                }

                foundBothStreams = true;
                fusionData.setStreamData(upstreamGenes, downstreamGenes, se == SE_START);
            }
            else if(!foundBothStreams && !foundOneStream && (!upstreamGenes.isEmpty() || !downstreamGenes.isEmpty()))
            {
                // take just one stream if clear
                foundOneStream = true;
                fusionData.setStreamData(upstreamGenes, downstreamGenes, se == SE_START);
            }
        }

        fusionData.cacheTranscriptData();

        initialFragment.setJunctionTypes(mConfig.RefFastaSeqFile, fusionData.getGeneStrands());

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
        // attempt to allocate discordant fragments to fusions
        for(Map.Entry<String,List<FusionFragment>> entry : mDiscordantFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionCandidates.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final Set<FusionFragment> allocatedFragments = Sets.newHashSet();

            for (FusionFragment fragment : fragments)
            {
                // only allocate discordant fragments if they support a gene & transcript at both ends
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                for(FusionReadData fusionData : fusions)
                {
                    if(fusionData.canAddDiscordantFragment(fragment, mConfig.MaxFragmentLength))
                    {
                        // check if one of the split read ends can be realigned to support the fusion junction
                        if(fusionData.canRelignFragmentToJunction(fragment))
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

    private void createLocalFusions()
    {
        // create fusions from fragments with 1 or both junctions matching known splice sites between genes without supp alignment
        // and then reassign any other fragments to these new fusions

        final Set<FusionReadData> newFusions = Sets.newHashSet();

        for(Map.Entry<String,List<FusionFragment>> entry : mDiscordantFragments.entrySet())
        {
            final List<FusionFragment> fragments = entry.getValue();

            final Set<FusionFragment> allocatedFragments = Sets.newHashSet();

            for (FusionFragment fragment : fragments)
            {
                if(fragment.type() == MATCHED_JUNCTION)
                {
                    FusionReadData fusionData = createOrUpdateFusion(fragment);

                    newFusions.add(fusionData);
                    allocatedFragments.add(fragment);
                }
            }

            allocatedFragments.forEach(x -> fragments.remove(x));
        }

        // now re-check remaining non-split-junction against the new fusions only
        for(Map.Entry<String,List<FusionFragment>> entry : mDiscordantFragments.entrySet())
        {
            final List<FusionReadData> fusions = newFusions.stream()
                    .filter(x -> x.locationId().equals(entry.getKey())).collect(Collectors.toList());

            if(fusions.isEmpty())
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final Set<FusionFragment> allocatedFragments = Sets.newHashSet();

            for (FusionFragment fragment : fragments)
            {
                // only allocate discordant fragments if they support a gene & transcript at both ends
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                for(FusionReadData fusionData : fusions)
                {
                    if(fusionData.canAddDiscordantFragment(fragment, mConfig.MaxFragmentLength))
                    {
                        // check if one of the split read ends can be realigned to support the fusion junction
                        if(fusionData.canRelignFragmentToJunction(fragment))
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
        for(Map.Entry<String,List<FusionFragment>> entry : mRealignCandidateFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionsByGene.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final Set<FusionFragment> allocatedFragments = Sets.newHashSet();

            for (FusionFragment fragment : fragments)
            {
                for(FusionReadData fusionData : fusions)
                {
                    if(fusionData.canRelignFragmentToJunction(fragment))
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
        mRealignCandidateFragments.clear();
        mReadIds.clear();
    }

}

