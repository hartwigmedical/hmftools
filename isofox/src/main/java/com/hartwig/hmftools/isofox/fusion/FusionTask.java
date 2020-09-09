package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;

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

    private int mNextFusionId;
    private final Map<String,List<FusionReadData>> mFusionCandidates; // keyed by the chromosome pair
    private Map<String,Map<String,FusionReadData>> mFusionsByLocation; // keyed by the chromosome pair, then precise position (hashed)
    private final Map<String,List<FusionReadData>> mFusionsByGene; // keyed by the locationId
    private final Map<String,List<FusionFragment>> mDiscordantFragments; // keyed by the chromosome pair
    private final Map<String,List<FusionFragment>> mRealignCandidateFragments; // keyed by chromosome, since single-sided

    private final FusionWriter mFusionWriter;

    private final PerformanceCounter mPerfCounter;

    public FusionTask(
            int taskId, final IsofoxConfig config, final EnsemblDataCache geneTransCache,
            final Map<String,Map<Integer,BaseDepth>> geneDepthMap, final FusionWriter fusionWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mGeneTransCache = geneTransCache;
        mChrGeneDepthMap = geneDepthMap;

        mAllFragments = Lists.newArrayList();

        mNextFusionId = 0;
        mFusionCandidates = Maps.newHashMap();
        mFusionsByLocation = Maps.newHashMap();
        mFusionsByGene = Maps.newHashMap();
        mDiscordantFragments = Maps.newHashMap();
        mRealignCandidateFragments = Maps.newHashMap();

        mFusionWriter = fusionWriter;

        mPerfCounter = new PerformanceCounter("FusionTask");
   }

    public final Map<String,List<FusionReadData>> getFusionCandidates() { return mFusionCandidates; }
    public final Map<String,List<FusionFragment>> getUnfusedFragments() { return mDiscordantFragments; }
    public final List<FusionFragment> getFragments() { return mAllFragments; }

    @Override
    public Long call()
    {
        int initialFragmentCount = mAllFragments.size();

        ISF_LOGGER.info("{}: processing {} chimeric fragments", mTaskId, initialFragmentCount);

        mPerfCounter.start();;

        formInitialFusions();
        reconcileFusions();

        // assign any discordant reads
        ISF_LOGGER.debug("{}: assigning unfused fragments", mTaskId);
        assignDiscordantFragments();
        createLocalFusions();
        assignRealignedFragments();

        annotateFusions();

        mPerfCounter.stop();

        writeData();

        ISF_LOGGER.info("{}: fusion task complete, fusions({}) unassigned(disc={} realgn={})",
                mTaskId, mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum(),
                mDiscordantFragments.values().stream().mapToInt(x -> x.size()).sum(),
                mRealignCandidateFragments.values().stream().mapToInt(x -> x.size()).sum());

        return (long)1;
    }

    public final PerformanceCounter getPerfCounter() { return mPerfCounter; }

    private void formInitialFusions()
    {
        int junctioned = 0;
        for (FusionFragment fragment : mAllFragments)
        {
            if (fragment.type() == MATCHED_JUNCTION)
            {
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

        // free up the set of initial fragments now they've all been assigned
        mAllFragments.clear();
        mFusionsByLocation.clear();
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

    private FusionReadData findExistingFusion(final FusionFragment fragment)
    {
        final Map<String,FusionReadData> fusionsByPosition = mFusionsByLocation.get(formChromosomePair(fragment.chromosomes()));

        if(fusionsByPosition == null)
            return null;

        return fusionsByPosition.get(fragment.positionHash());
    }

    private FusionReadData createOrUpdateFusion(final FusionFragment fragment)
    {
        // scenarios:
        // 1. New fusion with correct splice-junction support - may or may not match a known transcript and exon
        // 2. Potential discordant or realigned fragment

        // fusions will be stored in a map keyed by their location pair (chromosome + geneCollectionId)
        // and in an additional map of precise positions to avoid mismatches on gene collections
        FusionReadData existingFusion = findExistingFusion(fragment);

        if(existingFusion != null)
        {
            existingFusion.addFusionFragment(fragment);

            // mark donor-acceptor types whether strands are known or not
            fragment.junctionTypes()[SE_START] = existingFusion.getInitialFragment().junctionTypes()[SE_START];
            fragment.junctionTypes()[SE_END] = existingFusion.getInitialFragment().junctionTypes()[SE_END];
            return existingFusion;
        }

        List<FusionReadData> fusions = mFusionCandidates.get(fragment.locationPair());

        if(fusions == null)
        {
            fusions = Lists.newArrayList();
            mFusionCandidates.put(fragment.locationPair(), fusions);
        }

        int fusionId = mTaskId * 1000000 + mNextFusionId++; // keep unique across threaded tasks
        final FusionReadData fusionData = new FusionReadData(fusionId, fragment);

        fusionData.setJunctionBases(mConfig.RefGenome);

        setGeneData(fusionData);

        fusions.add(fusionData);

        // add to precise-location store
        Map<String,FusionReadData> fusionsByPosition = mFusionsByLocation.get(formChromosomePair(fragment.chromosomes()));

        if(fusionsByPosition == null)
        {
            fusionsByPosition = Maps.newHashMap();
            mFusionsByLocation.put(formChromosomePair(fragment.chromosomes()), fusionsByPosition);
        }

        fusionsByPosition.put(fragment.positionHash(), fusionData);

        // add to the cache by gene for later assignment or realigned fragments
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
        final List<EnsemblGeneData>[] genesByPosition = new List[] { Lists.newArrayList(), Lists.newArrayList() };
        final List<TranscriptData>[] validTransDataList = new List[] { Lists.newArrayList(), Lists.newArrayList() };

        FusionFragment initialFragment = fusionData.getInitialFragment();

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<String> spliceGeneIds = initialFragment.getGeneIds(se);

            // purge any invalid transcript exons when the splice junction is known
            final List<TranscriptData> transDataList = Lists.newArrayList();
            spliceGeneIds.forEach(x -> transDataList.addAll(mGeneTransCache.getTranscripts(x)));

            // purge any invalid transcript-exons and mark the junction as known if applicable
            initialFragment.validateTranscriptExons(transDataList, se);

            // and collect again
            spliceGeneIds.clear();
            spliceGeneIds.addAll(initialFragment.getGeneIds(se));

            if(!spliceGeneIds.isEmpty())
            {
                genesByPosition[se] = spliceGeneIds.stream().map(x -> mGeneTransCache.getGeneDataById(x)).collect(Collectors.toList());

                final int seIndex = se;
                validTransDataList[se] = transDataList.stream()
                        .filter(x -> initialFragment.getTransExonRefs()[seIndex].stream().anyMatch(y -> x.TransId == y.TransId))
                        .collect(Collectors.toList());
            }
        }

        if(genesByPosition[SE_START].size() > 1 || genesByPosition[SE_END].size() > 1)
        {
            boolean matched = false;

            for(EnsemblGeneData gene1 : genesByPosition[SE_START])
            {
                for(EnsemblGeneData gene2 : genesByPosition[SE_END])
                {
                    if(mConfig.Fusions.KnownFusions.hasKnownFusion(gene1.GeneName, gene2.GeneName)
                    || mConfig.Fusions.KnownFusions.hasKnownFusion(gene2.GeneName, gene1.GeneName))
                    {
                        matched = true;
                        genesByPosition[SE_START].clear();
                        genesByPosition[SE_START].add(gene1);
                        genesByPosition[SE_END].clear();
                        genesByPosition[SE_END].add(gene2);
                        break;
                    }
                }

                if(matched)
                    break;
            }

            if(!matched)
            {
                for(int se = SE_START; se <= SE_END; ++se)
                {
                    if(genesByPosition[se].size() > 1)
                        prioritiseKnownFusionGene(genesByPosition[se]);
                }
            }
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

            final List<EnsemblGeneData> upstreamGenes = genesByPosition[upstreamIndex].stream()
                    .filter(x -> x.Strand == sjOrientations[upstreamIndex]).collect(Collectors.toList());

            final List<EnsemblGeneData> downstreamGenes = genesByPosition[downstreamIndex].stream()
                    .filter(x -> x.Strand == -sjOrientations[downstreamIndex]).collect(Collectors.toList());

            // where multiple (non-known) genes are still considered, take the longest coding one
            prioritiseLongestCodingFusionGene(upstreamGenes, validTransDataList[upstreamIndex]);
            prioritiseLongestCodingFusionGene(downstreamGenes, validTransDataList[downstreamIndex]);

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

        initialFragment.setJunctionTypes(mConfig.RefGenome, fusionData.getGeneStrands(), fusionData.junctionSpliceBases());
    }

    private void prioritiseLongestCodingFusionGene(final List<EnsemblGeneData> geneList, final List<TranscriptData> transDataList)
    {
        if(geneList.isEmpty())
            return;
        else if(geneList.size() == 1)
            return;

        // take the longest protein coding
        EnsemblGeneData longestGene = null;
        int longestLength = 0;

        for(EnsemblGeneData geneData : geneList)
        {
            int maxCodingBases = transDataList.stream()
                    .filter(x -> x.GeneId.equals(geneData.GeneId))
                    .filter(x -> x.CodingStart != null)
                    .mapToInt(x -> x.CodingEnd - x.CodingStart)
                    .max().orElse(0);

            if(maxCodingBases > longestLength)
            {
                longestLength = maxCodingBases;
                longestGene = geneData;
            }
        }

        geneList.clear();
        geneList.add(longestGene);
    }

    private void prioritiseKnownFusionGene(final List<EnsemblGeneData> geneList)
    {
        if(geneList.isEmpty())
            return;
        else if(geneList.size() == 1)
            return;

        // first look for a known pair gene
        List<EnsemblGeneData> culledList = geneList.stream()
                .filter(x -> mConfig.Fusions.KnownFusions.hasKnownPairGene(x.GeneName))
                .collect(Collectors.toList());

        if(culledList.isEmpty())
        {
            culledList = geneList.stream()
                    .filter(x -> mConfig.Fusions.KnownFusions.hasPromiscuousFiveGene(x.GeneName)
                            || mConfig.Fusions.KnownFusions.hasPromiscuousThreeGene(x.GeneName))
                    .collect(Collectors.toList());
        }

        if(culledList.size() == 1)
        {
            geneList.clear();
            geneList.add(culledList.get(0));
            return;
        }
    }

    private void reconcileFusions()
    {
        // merge fusions if they match on homology
        // otherwise if fusion junctions are close enough, take the one with largest amount of support
        for(List<FusionReadData> fusions : mFusionCandidates.values())
        {
            if(fusions.size() == 1)
                continue;

            // annotate similar fusions for post-run analysis
            for(int i = 0; i < fusions.size() - 1; ++i)
            {
                FusionReadData fusion1 = fusions.get(i);

                int j = i + 1;
                while(j < fusions.size())
                {
                    FusionReadData fusion2 = fusions.get(j);

                    if(fusion1.matchWithinHomology(fusion2))
                    {
                        ISF_LOGGER.debug("fusion({}) homology merge with fusion({})", fusion1.id(), fusion2.id());

                        ISF_LOGGER.debug("fusion1({}) homology({}/{}) start(junc={} adj={}) end(junc={} adj={})",
                                fusion1.toString(), fusion1.junctionHomology()[SE_START], fusion1.junctionHomology()[SE_END],
                                fusion1.junctionBases()[SE_START], fusion1.adjacentJunctionBases()[SE_START],
                                fusion1.junctionBases()[SE_END], fusion1.adjacentJunctionBases()[SE_END]);

                        ISF_LOGGER.debug("fusion2({}) homology({}/{}) start(junc={} adj={}) end(junc={} adj={})",
                                fusion2.toString(), fusion2.junctionHomology()[SE_START], fusion2.junctionHomology()[SE_END],
                                fusion2.junctionBases()[SE_START], fusion2.adjacentJunctionBases()[SE_START],
                                fusion2.junctionBases()[SE_END], fusion2.adjacentJunctionBases()[SE_END]);

                        final FusionReadData fusion1Const = fusion1;
                        fusion2.getFragments(MATCHED_JUNCTION).forEach(x -> fusion1Const.addFusionFragment(x));
                        fusions.remove(j);
                        continue;
                    }

                    if(fusion1.isCloseMatch(fusion2))
                    {
                        // keep the fusion with more support, discard the other
                        if(fusion1.getFragments(MATCHED_JUNCTION).size() < fusion2.getFragments(MATCHED_JUNCTION).size())
                        {
                            fusions.set(i, fusion2);
                            fusion1 = fusion2;
                        }

                        fusions.remove(j);
                    }
                    else
                    {
                        ++j;
                    }
                }
            }
        }
    }

    private static final int POSITION_REALIGN_DISTANCE = 20;

    private void markRelatedFusions()
    {
        for( Map.Entry<String,List<FusionReadData>> entry : mFusionCandidates.entrySet())
        {
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
                            fusion2.addRelatedFusion(fusion1.id(), true);
                            continue;
                        }
                    }
                    else if(isUnspliced && fusion2.isKnownSpliced())
                    {
                        if(hasTranscriptExonMatch(upRefs1, fusion2.getTransExonRefsByStream(FS_UPSTREAM))
                        && hasTranscriptExonMatch(fusion2.getTransExonRefsByStream(FS_DOWNSTREAM), downRefs1, -1))
                        {
                            fusion1.addRelatedFusion(fusion2.id(), true);
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
                        fusion1.addRelatedFusion(fusion2.id(), false);
                        fusion2.addRelatedFusion(fusion1.id(), false);
                    }
                }
            }
        }
    }

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
    }

}

