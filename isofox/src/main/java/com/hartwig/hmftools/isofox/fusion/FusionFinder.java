package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.ReadRecord.NO_GENE_ID;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionConfig.LOG_READ_ID;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.softClippedReadSupportsJunction;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.checkMissingGeneData;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.ReadGroup.mergeChimericReadMaps;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionFinder implements Callable
{
    private final String mTaskId;
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final List<FusionFragment> mAllFragments;

    private final List<ReadGroup> mSpanningReadGroups; // temporary caching for read groups spanning gene collections
    private final Map<String,ReadGroup> mChimericPartialReadGroups;

    private final Map<String,List<FusionReadData>> mFusionCandidates; // keyed by the chromosome pair
    private final Map<String,Map<String,FusionReadData>> mFusionsByLocation; // keyed by the chromosome pair, then precise position (hashed)
    private final Map<String,List<FusionReadData>> mFusionsByGene; // keyed by the chr + gene collectionId
    private final Map<String,List<FusionFragment>> mDiscordantFragments; // keyed by the chromosome pair
    private final Map<String,List<FusionFragment>> mRealignCandidateFragments; // keyed by chr + gene collectionId

    private final FusionWriter mFusionWriter;
    private final FusionGeneFilters mGeneFilters;

    private final PerformanceCounter mPerfCounter;

    public FusionFinder(
            final String taskId, final IsofoxConfig config, final EnsemblDataCache geneTransCache,
            final FusionGeneFilters fusionGeneFilters, final FusionWriter fusionWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mAllFragments = Lists.newArrayList();

        mChimericPartialReadGroups = Maps.newHashMap();
        mSpanningReadGroups = Lists.newArrayList();

        mFusionCandidates = Maps.newHashMap();
        mFusionsByLocation = Maps.newHashMap();
        mFusionsByGene = Maps.newHashMap();
        mDiscordantFragments = Maps.newHashMap();
        mRealignCandidateFragments = Maps.newHashMap();

        mFusionWriter = fusionWriter;
        mGeneFilters = fusionGeneFilters;

        mPerfCounter = new PerformanceCounter("FusionTask");
   }

    public final Map<String,List<FusionReadData>> getFusionCandidates() { return mFusionCandidates; }
    public final Map<String,List<FusionFragment>> getUnfusedFragments() { return mDiscordantFragments; }
    public final Map<String,ReadGroup> getChimericPartialReadGroups() { return mChimericPartialReadGroups; }
    public final Map<String,List<FusionFragment>> getRealignCandidateFragments() { return mRealignCandidateFragments; }

    public final List<FusionFragment> getFragments() { return mAllFragments; }
    public final List<ReadGroup> getSpanningReadGroups() { return mSpanningReadGroups; }

    public void clearState()
    {
        mAllFragments.clear();
        mFusionCandidates.clear();
        mFusionsByLocation.clear();
        mFusionsByGene.clear();
        mDiscordantFragments.clear();
    }

    public List<ReadGroup> processNewChimericReadGroups(
            final GeneCollection geneCollection, final BaseDepth baseDepth, final Map<String,ReadGroup> newReadGroups)
    {
        List<ReadGroup> completeReadGroups = Lists.newArrayList();

        mergeChimericReadMaps(mChimericPartialReadGroups, completeReadGroups, newReadGroups);

        // identify any read groups with reads spanning into a future gene collection
        // and fill in any missing gene info for reads (partial or complete) which link to this gene collections
        final List<ReadGroup> spanningGroups = newReadGroups.values().stream()
                .filter(x -> x.Reads.stream().anyMatch(y -> y.getGeneCollectons()[SE_END] == NO_GENE_ID))
                .collect(Collectors.toList());

        final List<ReadGroup> geneCompletedGroups = reconcileSpanningReadGroups(geneCollection, spanningGroups, baseDepth);

        spanningGroups.stream().forEach(x -> completeReadGroups.remove(x));
        geneCompletedGroups.stream().filter(x -> !completeReadGroups.contains(x)).forEach(x -> completeReadGroups.add(x));

        return completeReadGroups;
    }

    private List<ReadGroup> reconcileSpanningReadGroups(
            final GeneCollection geneCollection, final List<ReadGroup> spanningReadGroups, final BaseDepth baseDepth)
    {
        List<ReadGroup> completeGroups = Lists.newArrayList();

        // check pending groups which needed their upper gene and depth info populated by this gene collection
        int index = 0;
        while(index < mSpanningReadGroups.size())
        {
            ReadGroup readGroup = mSpanningReadGroups.get(index);
            boolean missingGeneInfo = false;

            for(ReadRecord read : readGroup.Reads)
            {
                if(read.getGeneCollectons()[SE_END] != NO_GENE_ID)
                    continue;

                if(!positionWithin(read.getCoordsBoundary(SE_END),
                        geneCollection.getNonGenicPositions()[SE_START], geneCollection.getNonGenicPositions()[SE_END]))
                {
                    missingGeneInfo = true;
                    continue;
                }

                if(positionWithin(read.getCoordsBoundary(SE_END),
                        geneCollection.regionBounds()[SE_START], geneCollection.regionBounds()[SE_END]))
                {
                    read.setGeneCollection(SE_END, geneCollection.id(), true);
                    checkMissingGeneData(read, geneCollection.getTranscripts());
                }
                else
                {
                    read.setGeneCollection(SE_END, geneCollection.id(), false);
                }

                // fill in junction positions depth from reads which spanned into this GC (eg from long N-split reads)
                read.setReadJunctionDepth(baseDepth);
            }

            if(!missingGeneInfo)
            {
                mSpanningReadGroups.remove(index);

                if(readGroup.isComplete())
                    completeGroups.add(readGroup);
            }
            else
            {
                ++index;
            }
        }

        mSpanningReadGroups.addAll(spanningReadGroups);
        return completeGroups;
    }

    public Map<String,Map<String,ReadGroup>> extractIncompleteReadGroups(final String chromosome)
    {
        Map<String,Map<String,ReadGroup>> chrIncompleteReadsGroups = Maps.newHashMap();
        for(ReadGroup readGroup : mChimericPartialReadGroups.values())
        {
            String otherChromosome = readGroup.findOtherChromosome(chromosome);

            if(otherChromosome != null)
            {
                if(!HumanChromosome.contains(otherChromosome))
                    continue;

                Map<String, ReadGroup> readGroupMap = chrIncompleteReadsGroups.get(otherChromosome);
                if(readGroupMap == null)
                {
                    readGroupMap = Maps.newHashMap();
                    chrIncompleteReadsGroups.put(otherChromosome, readGroupMap);
                }

                readGroupMap.put(readGroup.id(), readGroup);
            }
        }

        mChimericPartialReadGroups.clear();
        return chrIncompleteReadsGroups;
    }

    public final Map<String,List<FusionFragment>> extractRealignCandidateFragments(final Map<String,Map<String,ReadGroup>> chrReadGroups)
    {
        // extract any RAC fragment which supports an (incomplete) inter-chromosomal read group
        final Map<String,List<FusionFragment>> matchingRacFragments = Maps.newHashMap();

        for(Map<String,ReadGroup> readGroupMap : chrReadGroups.values())
        {
            for(ReadGroup readGroup : readGroupMap.values())
            {
                ReadRecord saRead = readGroup.Reads.stream().filter(x -> x.hasSuppAlignment()).findFirst().orElse(null);

                if(saRead == null)
                    continue;

                final int scIndex = saRead.longestSoftClippedEnd();
                int junctionPosition = saRead.getCoordsBoundary(scIndex);
                byte junctionOrientation = scIndex == SE_START ? NEG_ORIENT : POS_ORIENT;

                // now find any potentially supporting RAC fragments
                ChrGeneCollectionPair chrGenePair = new ChrGeneCollectionPair(saRead.Chromosome, saRead.getGeneCollectons()[SE_START]);

                final List<FusionFragment> fragments = mRealignCandidateFragments.get(chrGenePair.toString());

                if(fragments == null)
                    continue;

                List<FusionFragment> matchingFragments = null;

                for(FusionFragment fragment : fragments)
                {
                    if(fragment.reads().stream().
                            anyMatch(x -> softClippedReadSupportsJunction(x, scIndex, junctionPosition, junctionOrientation, null)))
                    {
                        if(matchingFragments == null)
                        {
                            matchingFragments = Lists.newArrayList();
                            matchingRacFragments.put(chrGenePair.toString(), matchingFragments);
                        }

                        matchingFragments.add(fragment);
                    }
                }
            }
        }

        return matchingRacFragments;
    }

    public void processLocalReadGroups(final List<ReadGroup> readGroups)
    {
        processReadGroups(readGroups, false);
    }

    public void processInterChromosomalReadGroups(final List<ReadGroup> readGroups)
    {
        processReadGroups(readGroups, true);
    }

    private static final int LOG_COUNT = 10000;

    private void processReadGroups(final List<ReadGroup> readGroups, boolean isInterChromosomal)
    {
        // read groups are guaranteed to be complete

        clearState();

        if(readGroups.isEmpty())
            return;

        mPerfCounter.start();

        // first turn them into fragments, then look for fusions
        ISF_LOGGER.debug("chr({}) processing {} {} chimeric read groups",
                mTaskId, readGroups.size(), isInterChromosomal ? "inter-chromosome" : "local");

        int readGroupCount = 0;
        int nextLog = LOG_COUNT;

        for(ReadGroup readGroup : readGroups)
        {
            ++readGroupCount;

            if(readGroupCount >= nextLog)
            {
                nextLog += LOG_COUNT;
                ISF_LOGGER.info("chr({}) processed {} {} chimeric read groups",
                        mTaskId, isInterChromosomal ? "inter-chromosome" : "local", readGroupCount);
            }

            final List<ReadRecord> reads = readGroup.Reads;

            if(reads.get(0).Id.equals(LOG_READ_ID))
            {
                ISF_LOGGER.debug("specific read: {}", reads.get(0));
            }

            if(reads.stream().anyMatch(x -> mGeneFilters.skipRead(x.mateChromosome(), x.mateStartPosition())))
                continue;

            FusionFragment fragment = new FusionFragment(readGroup);

            if(fragment.type() == FusionFragmentType.UNKNOWN)
            {
                mFusionWriter.writeReadData(reads, "INVALID_FRAG");
                continue;
            }

            mAllFragments.add(fragment);
        }

        processFragments();

        if(!isInterChromosomal) // otherwise need to wait for RAC fragment assignment
        {
            assignRealignCandidateFragments(mRealignCandidateFragments);
            writeFusionSummary();
        }

        if(readGroupCount > LOG_COUNT)
        {
            ISF_LOGGER.info("{}: fusion task complete, fusions({}) unassigned(disc={} realgn={})",
                    mTaskId, mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum(),
                    mDiscordantFragments.values().stream().mapToInt(x -> x.size()).sum(),
                    mRealignCandidateFragments.values().stream().mapToInt(x -> x.size()).sum());
        }

        mPerfCounter.stop();
    }

    @Override
    public Long call()
    {
        mPerfCounter.start();;
        processFragments();
        assignRealignCandidateFragments(mRealignCandidateFragments);
        writeFusionSummary();
        mPerfCounter.stop();

        return (long)1;
    }

    private void processFragments()
    {
        // ISF_LOGGER.debug("{}: processing {} chimeric fragments", mTaskId, initialFragmentCount);

        formInitialFusions();
        reconcileFusions();

        // assign any discordant reads
        // ISF_LOGGER.debug("{}: assigning unfused fragments", mTaskId);
        assignDiscordantFragments();
        createLocalFusions();
    }

    private void writeFusionSummary()
    {
        annotateFusions();
        writeData();
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

        ISF_LOGGER.debug("{}: chimeric fragments({} disc={} candRealgn={} junc={}) fusions(loc={} gene={} total={})",
                mTaskId, mAllFragments.size(), mDiscordantFragments.values().stream().mapToInt(x -> x.size()).sum(),
                mRealignCandidateFragments.values().stream().mapToInt(x -> x.size()).sum(), junctioned,
                mFusionCandidates.size(), mFusionsByGene.size(), mFusionCandidates.values().stream().mapToInt(x -> x.size()).sum());

        // free up the set of initial fragments now they've all been assigned
        mAllFragments.clear();
        mFusionsByLocation.clear();
    }

    private void annotateFusions()
    {
        markRelatedFusions();
    }

    private void writeData()
    {
        // write results
        mFusionWriter.writeFusionData(mFusionCandidates);
        mFusionWriter.writeUnfusedFragments(mDiscordantFragments);
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
            existingFusion.addFusionFragment(fragment, mConfig.Fusions.CacheFragments);

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

        // int fusionId = mTaskId * 1000000 + mNextFusionId++; // keep unique across threaded tasks
        int fusionId = mFusionWriter.getNextFusionId();
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
            if(se == SE_END && fragment.isSingleGeneCollection())
                break;

            final String chrGenePair = fragment.chrGeneCollection(se).toString();
            List<FusionReadData> fusionsByGene = mFusionsByGene.get(chrGenePair);

            if(fusionsByGene == null)
                mFusionsByGene.put(chrGenePair, Lists.newArrayList(fusionData));
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

    private void setGeneData(final FusionReadData fusionData)
    {
        // get the genes supporting the splice junction in the terms of an SV (ie lower chromosome and lower position first)
        final List<GeneData>[] genesByPosition = new List[] { Lists.newArrayList(), Lists.newArrayList() };
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

            for(GeneData gene1 : genesByPosition[SE_START])
            {
                for(GeneData gene2 : genesByPosition[SE_END])
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

            final List<GeneData> upstreamGenes = genesByPosition[upstreamIndex].stream()
                    .filter(x -> x.Strand == sjOrientations[upstreamIndex]).collect(Collectors.toList());

            final List<GeneData> downstreamGenes = genesByPosition[downstreamIndex].stream()
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

    private void prioritiseLongestCodingFusionGene(final List<GeneData> geneList, final List<TranscriptData> transDataList)
    {
        if(geneList.isEmpty())
            return;
        else if(geneList.size() == 1)
            return;

        // take the longest protein coding
        GeneData longestGene = null;
        int longestLength = 0;

        for(GeneData geneData : geneList)
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

    private void prioritiseKnownFusionGene(final List<GeneData> geneList)
    {
        if(geneList.isEmpty())
            return;
        else if(geneList.size() == 1)
            return;

        // first look for a known pair gene
        List<GeneData> culledList = geneList.stream()
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

            for(int i = 0; i < fusions.size() - 1; ++i)
            {
                FusionReadData fusion1 = fusions.get(i);

                boolean knownSpliceJunction = fusion1.isKnownSpliced();

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

                        if(mConfig.Fusions.CacheFragments)
                            fusion2.getFragments(MATCHED_JUNCTION).forEach(x -> fusion1Const.addFusionFragment(x, mConfig.Fusions.CacheFragments));
                        else
                            fusion1Const.addFragmentTypeCount(MATCHED_JUNCTION, fusion2.getFragmentTypeCount(MATCHED_JUNCTION));

                        fusions.remove(j);
                        continue;
                    }

                    if(!(knownSpliceJunction && fusion2.isKnownSpliced()) && fusion1.isCloseMatch(fusion2))
                    {
                        // keep the fusion with more support, discard the other
                        if(fusion1.getFragmentTypeCount(MATCHED_JUNCTION) < fusion2.getFragmentTypeCount(MATCHED_JUNCTION))
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

                final List<TransExonRef> upRefs1 = fusion1.getTransExonRefsByStream(FS_UP);
                final List<TransExonRef> downRefs1 = fusion1.getTransExonRefsByStream(FS_DOWN);

                for(int j = i + 1; j < fusions.size() - 1; ++j)
                {
                    FusionReadData fusion2 = fusions.get(j);

                    if(isSpliced && fusion2.isUnspliced())
                    {
                        if(hasTranscriptExonMatch(upRefs1, fusion2.getTransExonRefsByStream(FS_UP))
                        && hasTranscriptExonMatch(downRefs1, fusion2.getTransExonRefsByStream(FS_DOWN), -1))
                        {
                            fusion2.addRelatedFusion(fusion1.id(), true);
                            continue;
                        }
                    }
                    else if(isUnspliced && fusion2.isKnownSpliced())
                    {
                        if(hasTranscriptExonMatch(upRefs1, fusion2.getTransExonRefsByStream(FS_UP))
                        && hasTranscriptExonMatch(fusion2.getTransExonRefsByStream(FS_DOWN), downRefs1, -1))
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

                        fusionData.addFusionFragment(fragment, mConfig.Fusions.CacheFragments);
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

                        fusionData.addFusionFragment(fragment, mConfig.Fusions.CacheFragments);
                        allocatedFragments.add(fragment);
                    }
                }
            }

            allocatedFragments.forEach(x -> fragments.remove(x));
        }
    }

    private void cacheRealignCandidateFragment(final FusionFragment fragment)
    {
        final String chrPair = fragment.chrGeneCollection(SE_START).toString();
        List<FusionFragment> fragments = mRealignCandidateFragments.get(chrPair);

        if(fragments == null)
            mRealignCandidateFragments.put(chrPair, Lists.newArrayList(fragment));
        else
            fragments.add(fragment);
    }

    private void assignRealignCandidateFragments(final Map<String,List<FusionFragment>> racFragments)
    {
        if(racFragments.isEmpty())
            return;

        for(Map.Entry<String,List<FusionReadData>> entry : mFusionsByGene.entrySet())
        {
            final ChrGeneCollectionPair chrGenePair = ChrGeneCollectionPair.from(entry.getKey());
            final List<FusionReadData> fusions = entry.getValue();

            final List<FusionFragment> realignCandidates = racFragments.get(chrGenePair.toString());

            if(realignCandidates == null || realignCandidates.isEmpty())
                continue;

            for(final FusionReadData fusionData : fusions)
            {
                for (FusionFragment fragment : realignCandidates)
                {
                    if(fusionData.canRelignFragmentToJunction(fragment))
                    {
                        fragment.setType(REALIGNED);
                        fusionData.addFusionFragment(fragment, mConfig.Fusions.CacheFragments);
                    }
                }
            }
        }
    }

    public void assignInterChromosomalRacFragments(final Map<String,List<FusionFragment>> racFragments)
    {
        mRealignCandidateFragments.clear();
        assignRealignCandidateFragments(racFragments);
        writeFusionSummary();
    }
}

