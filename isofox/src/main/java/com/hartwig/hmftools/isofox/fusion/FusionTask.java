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
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionTask implements Callable
{
    private final int mTaskId;
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    // private final Map<String,List<FusionFragment>> mChrPairFragments;
    private final List<FusionFragment> mAllFragments;
    private final Set<String> mReadIds;

    private int mNextFusionId;
    private final Map<String,List<FusionReadData>> mFusionCandidates; // keyed by the chromosome pair
    private final Map<String,List<FusionFragment>> mUnfusedFragments;

    private final FusionReadDepth mFusionReadDepth;
    private final FusionWriter mFusionWriter;

    private final PerformanceCounter[] mPerfCounters;

    private static final int PC_DISCOVERY = 0;
    private static final int PC_ANNOTATE = 1;
    private static final int PC_DEPTH = 2;

    public FusionTask(
            int taskId, final IsofoxConfig config, final EnsemblDataCache geneTransCache, final List<FusionFragment> fragments,
            final FusionWriter fusionWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mAllFragments = Lists.newArrayList(fragments);

        mNextFusionId = 0;
        mFusionCandidates = Maps.newHashMap();
        mUnfusedFragments = Maps.newHashMap();
        mReadIds = Sets.newHashSet();

        mFusionReadDepth = new FusionReadDepth(mConfig, mReadIds);
        mFusionWriter = fusionWriter;

        mPerfCounters = new PerformanceCounter[PC_DEPTH + 1];
        mPerfCounters[PC_DISCOVERY] = new PerformanceCounter("Discovery");
        mPerfCounters[PC_ANNOTATE] = new PerformanceCounter("Annotate");
        mPerfCounters[PC_DEPTH] = new PerformanceCounter("ReadDepth");
   }

    public final Map<String,List<FusionReadData>> getFusionCandidates() { return mFusionCandidates; }
    public final Map<String,List<FusionFragment>> getUnfusedFragments() { return mUnfusedFragments; }
    public final FusionReadDepth getFusionReadDepth() { return mFusionReadDepth; }

    @Override
    public Long call()
    {
        formInitialFusions();

        calculateReadDepth();

        annotateFusions();

        writeData();

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
            else
            {
                cacheUnfusedFragment(fragment);
            }
        }

        ISF_LOGGER.info("{}: chimeric fragments({} unfused={} junc={}) fusions(loc={} total={})",
                mTaskId, mAllFragments.size(), mUnfusedFragments.values().stream().mapToInt(x -> x.size()).sum(), junctioned,
                mFusionCandidates.size(), mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum());

        mPerfCounters[PC_DISCOVERY].stop();
    }

    private void calculateReadDepth()
    {
        mPerfCounters[PC_DEPTH].start();

        ISF_LOGGER.info("{}: calculating junction depth", mTaskId);

        // classify / analyse fusions
        int nextLog = 1000;
        int fusionCount = 0;
        for (List<FusionReadData> fusions : mFusionCandidates.values())
        {
            mFusionReadDepth.calcFusionReadDepth(fusions);

            if (++fusionCount >= nextLog)
            {
                ISF_LOGGER.info("processed {} fusions", fusionCount);
                nextLog += 1000;
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
        assignUnfusedFragments();

        mPerfCounters[PC_ANNOTATE].stop();
    }

    private void writeData()
    {
        // write results
        // ISF_LOGGER.debug("{}: writing results", mTaskId);
        mFusionWriter.writeFusionData(mFusionCandidates);
        mFusionWriter.writeUnfusedFragments(mUnfusedFragments);
    }

    private void cacheUnfusedFragment(final FusionFragment fragment)
    {
        List<FusionFragment> fragments = mUnfusedFragments.get(fragment.locationPair());

        if(fragments == null)
        {
            fragments = Lists.newArrayList();
            mUnfusedFragments.put(fragment.locationPair(), fragments);
        }

        fragments.add(fragment);
    }

    private void createOrUpdateFusion(final FusionFragment fragment)
    {
        // scenarios:
        // 1. New fusion with correct splice-junction support - may or may not match a known transcript and exon
        // 3. Potential discordant or realigned fragment

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
                    ISF_LOGGER.warn("fusion({}) has multiple gene pairings by strand and orientation", fusionData.toString());
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

    private void assignUnfusedFragments()
    {
        for(Map.Entry<String,List<FusionFragment>> entry : mUnfusedFragments.entrySet())
        {
            final List<FusionReadData> fusions = mFusionCandidates.get(entry.getKey());

            if(fusions == null)
                continue;

            final List<FusionFragment> fragments = entry.getValue();

            final List<FusionFragment> allocatedFragments = Lists.newArrayList();

            for (FusionFragment fragment : fragments)
            {
                if(fragment.getTransExonRefs()[SE_START].isEmpty() || fragment.getTransExonRefs()[SE_END].isEmpty())
                    continue;

                for(FusionReadData fusion : fusions)
                {
                    if(!fusion.isValid())
                        continue;

                    if(fusion.canAddUnfusedFragment(fragment, mConfig.MaxFragmentLength))
                    {
                        allocatedFragments.add(fragment);

                        if(fusion.isRelignedFragment(fragment)) // was also fragment.isSingleGene()
                        {
                            fragment.setType(REALIGNED);
                        }
                        else
                        {
                            fragment.setType(DISCORDANT);
                        }

                        fusion.addFusionFragment(fragment);
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
        mUnfusedFragments.clear();
        mReadIds.clear();
    }

}

