package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_NOVEL_SJ_DISTANCE;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.addChimericReads;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitRead;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitReadJunction;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.hasRealignableSoftClip;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.isInversion;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.setHasMultipleKnownSpliceGenes;
import static com.hartwig.hmftools.isofox.fusion.LocalJunctionData.setMaxSplitMappedLength;
import static com.hartwig.hmftools.isofox.fusion.ReadGroup.hasSuppAlignment;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;

public class ChimericReadTracker
{
    private final IsofoxConfig mConfig;
    private final boolean mRunFusions;
    private final boolean mEnabled; // required for alt splice junctions

    private final List<String[]> mKnownPairGeneIds;

    private GeneCollection mGeneCollection; // the current collection being processed
    private final Map<String,ReadGroup> mChimericReadMap;
    private final List<ReadGroup> mLocalCompleteGroups;

    // junction position from fusion junction candidate reads are cached to identify candidate realignable reads
    private final Set<Integer> mJunctionPositions;

    private final List<List<ReadRecord>> mLocalChimericReads; // fragments to re-evaluate as alternate splice sites
    private final Map<String,ReadGroup> mCandidateRealignedReadMap;

    // to avoid double-processing reads falling after a gene collection
    private final Map<String,List<ReadRecord>> mPostGeneReadMap;
    private final Map<String,List<ReadRecord>> mPreviousPostGeneReadMap;
    private final ChimericStats mChimericStats;

    public ChimericReadTracker(final IsofoxConfig config)
    {
        mConfig = config;
        mRunFusions = mConfig.Functions.contains(FUSIONS);
        mEnabled = mRunFusions || mConfig.Functions.contains(ALT_SPLICE_JUNCTIONS);

        mKnownPairGeneIds = Lists.newArrayList();
        mChimericStats = new ChimericStats();
        mChimericReadMap = Maps.newHashMap();
        mJunctionPositions = Sets.newHashSet();
        mLocalChimericReads = Lists.newArrayList();
        mLocalCompleteGroups = Lists.newArrayList();
        mCandidateRealignedReadMap = Maps.newHashMap();
        mPostGeneReadMap = Maps.newHashMap();
        mPreviousPostGeneReadMap = Maps.newHashMap();
        mGeneCollection = null;
    }

    public boolean enabled() { return mEnabled; }

    public final Map<String,ReadGroup> getReadMap() { return mChimericReadMap; }
    public final Set<Integer> getJunctionPositions() { return mJunctionPositions; }
    public final List<List<ReadRecord>> getLocalChimericReads() { return mLocalChimericReads; }
    public ChimericStats getStats() { return mChimericStats; }

    public boolean isChimeric(final ReadRecord read1, final ReadRecord read2, boolean isDuplicate, boolean isMultiMapped)
    {
        if(read1.isChimeric() || read2.isChimeric() || !read1.withinGeneCollection() || !read2.withinGeneCollection())
            return true;

        if(!isDuplicate && !isMultiMapped && enabled() && (read1.containsSplit() || read2.containsSplit()))
        {
            return setHasMultipleKnownSpliceGenes(Lists.newArrayList(read1, read2), mKnownPairGeneIds);
        }

        return false;
    }

    public void registerKnownFusionPairs(final EnsemblDataCache geneTransCache)
    {
        for(final KnownFusionData knownPair : mConfig.Fusions.KnownFusions.getDataByType(KNOWN_PAIR))
        {
            final GeneData upGene = geneTransCache.getGeneDataByName(knownPair.FiveGene);
            final GeneData downGene = geneTransCache.getGeneDataByName(knownPair.ThreeGene);

            if(upGene != null && downGene != null)
                mKnownPairGeneIds.add(new String[] { upGene.GeneId, downGene.GeneId });
        }
    }

    public void initialise(final GeneCollection geneCollection)
    {
        mGeneCollection = geneCollection;

        mPreviousPostGeneReadMap.clear();
        mPreviousPostGeneReadMap.putAll(mPostGeneReadMap);
        mPostGeneReadMap.clear();

        // only purge junction positions which are now outside the regions to be processed
        Set<Integer> pastJuncPositions = mJunctionPositions.stream()
                .filter(x -> x < geneCollection.getNonGenicPositions()[SE_START]).collect(Collectors.toSet());

        pastJuncPositions.forEach(x -> mJunctionPositions.remove(x));
    }

    public void clear() { clear(false); }
    public void clearAll() { clear(true); }

    private void clear(boolean full)
    {
        mChimericReadMap.clear();
        mLocalCompleteGroups.clear();
        mCandidateRealignedReadMap.clear();
        mChimericStats.clear();
        mLocalChimericReads.clear();

        if(full)
        {
            mPreviousPostGeneReadMap.clear();
            mPostGeneReadMap.clear();
            mJunctionPositions.clear();
        }
    }

    public void addRealignmentCandidates(final ReadRecord read1, final ReadRecord read2)
    {
        if(read1.isDuplicate() || read2.isDuplicate()) // group complete so drop these
            return;

        mCandidateRealignedReadMap.put(read1.Id, new ReadGroup(read1, read2));
    }

    public void addChimericReadPair(final ReadRecord read1, final ReadRecord read2)
    {
        if(inExcludedRegion(read1, false) || inExcludedRegion(read2, false))
            return;

        if(!read1.isDuplicate() && !read2.isDuplicate())
        {
            // populate transcript info for intronic reads since it will be used in fusion matching
            addIntronicTranscriptData(read1);
            addIntronicTranscriptData(read2);
        }

        // add the pair when it's clear there aren't others with the same ID in the map
        if(mConfig.RunValidations && mChimericReadMap.containsKey(read1.Id))
        {
            // shouldn't occur
            ISF_LOGGER.error("overriding chimeric read({})", read1.Id);

            final ReadGroup existingGroup = mChimericReadMap.get(read1.Id);

            for(ReadRecord read : existingGroup.Reads)
            {
                ISF_LOGGER.error("existing read: {}", read);
            }

            ISF_LOGGER.error("new read: {}", read1);
            ISF_LOGGER.error("new read: {}", read2);

            existingGroup.Reads.add(read1);
            existingGroup.Reads.add(read2);
        }
        else
        {
            // cache info about any local dual-junction group
            ReadGroup group = new ReadGroup(read1, read2);

            if(group.size() == 2 && group.isComplete() && !group.hasDuplicateRead())
            {
                final int[] junctPositions = findCandidateJunctions(group.Reads);
                if(junctPositions[SE_START] > 0 && junctPositions[SE_END] > 0)
                {
                    final byte[] junctOrientations = {1, -1};
                    matchOrAddLocalJunctionGroup(group, junctPositions, junctOrientations);
                    return;
                }
            }

            mChimericReadMap.put(read1.Id, group);
        }
    }

    private void matchOrAddLocalJunctionGroup(final ReadGroup group, final int[] junctPositions, final byte[] junctOrientations)
    {
        LocalJunctionData matchData = null;
        for(ReadGroup readGroup : mLocalCompleteGroups)
        {
            if(!Arrays.equals(readGroup.localJunctionData().JunctionPositions, junctPositions))
                continue;

            if(!Arrays.equals(readGroup.localJunctionData().JunctionOrientations, junctOrientations))
                continue;

            matchData = readGroup.localJunctionData();
            ++matchData.MatchedGroupCount;
            ++mChimericStats.MatchedJunctions;
            break;
        }

        if(matchData == null)
        {
            matchData = new LocalJunctionData(junctPositions, junctOrientations);
            group.setLocalJunctionData(matchData);
            mLocalCompleteGroups.add(group);
            mChimericReadMap.put(group.id(), group);
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            setMaxSplitMappedLength(
                    se, group.Reads, junctPositions, junctOrientations, matchData.MaxSplitLengths);
        }
    }

    private void addIntronicTranscriptData(final ReadRecord read)
    {
        if(read.overlapsGeneCollection() && read.getMappedRegions().isEmpty())
            read.addIntronicTranscriptRefs(mGeneCollection.getTranscripts());
    }

    public void postProcessChimericReads(final BaseDepth baseDepth, final FragmentTracker fragmentTracker)
    {
        // check any lone reads - this cannot be one of a pair of non-genic reads since they will have already been dismissed
        // so will either be a supplementary or a read linked to another gene collection
        for(Object object : fragmentTracker.getValues())
        {
            final ReadRecord read = (ReadRecord)object;

            if(read.isMateUnmapped() || inExcludedRegion(read, true) || read.isSecondaryAlignment())
                continue;

            if(!read.isDuplicate())
            {
                baseDepth.processRead(read.getMappedRegionCoords());
                addIntronicTranscriptData(read);
            }

            addChimericReads(mChimericReadMap, read);
        }

        // migrate any local chimeric fragments for analysis as alternate splice junctions
        final List<String> fragsToRemove = Lists.newArrayList();

        for(final ReadGroup readGroup : mChimericReadMap.values())
        {
            // skip reads if all will be processed later or have been already
            final List<ReadRecord> reads = readGroup.Reads;
            final String readId = reads.get(0).Id;

            /*
            if(readId.equals(""))
            {
                ISF_LOGGER.debug("specific read: {}", readId);
            }
            */

            int readCount = reads.size();
            boolean readGroupComplete = readGroup.isComplete();

            // duplicates are kept until a group is complete, since not all reads are marked as duplicates and those would otherwise
            // be orphaned later on

            if(reads.stream().anyMatch(x -> x.isDuplicate()) && reads.size() >= 2)
            {
                // chimeric read groups with duplicates will be dropped later on once the group is complete
                mGeneCollection.addCount(DUPLICATE, 1);
            }

            if(mRunFusions && skipNonGenicReads(reads))
            {
                fragsToRemove.add(readId);
                continue;
            }

            boolean readsRemoved = reads.size() < readCount;

            if(!readsRemoved && !keepChimericGroup(reads, readGroupComplete))
            {
                fragsToRemove.add(readId);
                continue;
            }

            if(mRunFusions)
                collectCandidateJunctions(readGroup);
        }

        if(!fragsToRemove.isEmpty())
            fragsToRemove.forEach(x -> mChimericReadMap.remove(x));

        mChimericStats.ChimericJunctions += mJunctionPositions.size();

        int chimericCount = mChimericReadMap.size();
        mGeneCollection.addCount(TOTAL, chimericCount);
        mGeneCollection.addCount(CHIMERIC, chimericCount);

        if(mRunFusions)
        {
            addRealignCandidates();

            // chimeric reads will be processed by the fusion-finding routine, so need to capture transcript and exon data
            // and free up other gene & region read data (to avoid retaining large numbers of references/memory)
            for(final ReadGroup readGroup : mChimericReadMap.values())
            {
                readGroup.Reads.forEach(x -> x.captureGeneInfo(true));
                readGroup.Reads.forEach(x -> x.setReadJunctionDepth(baseDepth));
            }
        }
        else
        {
            // clear other chimeric state except for local junction information
            mChimericReadMap.clear();
            mLocalCompleteGroups.clear();
            mCandidateRealignedReadMap.clear();
            mChimericStats.clear();
        }
    }

    private boolean inExcludedRegion(final ReadRecord read, boolean checkMate)
    {
        // check the read and its supplementary data if present
        if(mConfig.Filters.skipRead(read.Chromosome, read.PosStart))
            return true;

        if(checkMate && mConfig.Filters.skipRead(read.mateChromosome(), read.mateStartPosition()))
            return true;

        // only skip fragments in immune regions if both junction positions are in one
        boolean inImmuneRegion = mConfig.Filters.ImmuneGeneRegions.stream()
                .anyMatch(x -> x.Chromosome.equals(read.Chromosome) && positionsOverlap(read.PosStart, read.PosEnd, x.start(), x.end()));

        if(inImmuneRegion
        && mConfig.Filters.ImmuneGeneRegions.stream().anyMatch(x -> x.containsPosition(read.mateChromosome(), read.mateStartPosition())))
        {
            return true;
        }

        if(read.hasSuppAlignment())
        {
            SupplementaryReadData suppData = SupplementaryReadData.from(read.getSuppAlignment());

            if(suppData != null && mConfig.Filters.skipRead(suppData.Chromosome, suppData.Position))
                return true;

            if(inImmuneRegion && suppData != null
            && mConfig.Filters.ImmuneGeneRegions.stream().anyMatch(x -> x.containsPosition(suppData.Chromosome, suppData.Position)))
            {
                return true;
            }
        }

        return false;
    }

    private boolean keepChimericGroup(final List<ReadRecord> reads, boolean readGroupComplete)
    {
        if(reads.stream().anyMatch(x -> x.isTranslocation()))
        {
            ++mChimericStats.Translocations;
            return true;
        }

        if(readGroupComplete && isInversion(reads))
        {
            ++mChimericStats.Inversions;
            return true;
        }

        boolean spanGeneSpliceSites = reads.stream().anyMatch(x -> x.hasInterGeneSplit()) ?
                true : setHasMultipleKnownSpliceGenes(reads, mKnownPairGeneIds);

        if(reads.stream().anyMatch(x -> x.spansGeneCollections()))
        {
            // may turn out to just end in the next pre-gene section but cannot say at this time
            ++mChimericStats.LocalInterGeneFrags;
            return true;
        }

        if(!readGroupComplete)
            return true;

        if(reads.stream().anyMatch(x -> !x.withinGeneCollection()))
        {
            // some reads are non-genic in full or part
            if(reads.stream().filter(x -> x.fullyNonGenic()).count() == reads.size())
            {
                // all reads non-genic - drop these entirely
                return false;
            }

            int minPosition = reads.stream().mapToInt(x -> x.getCoordsBoundary(SE_START)).min().orElse(0);
            int maxPosition = reads.stream().mapToInt(x -> x.getCoordsBoundary(SE_END)).max().orElse(0);

            if(mGeneCollection.regionBounds()[SE_START] - minPosition > MAX_NOVEL_SJ_DISTANCE
            || maxPosition - mGeneCollection.regionBounds()[SE_END] > MAX_NOVEL_SJ_DISTANCE)
            {
                // too far from the gene boundaries so consider these chimeric
                return true;
            }
        }

        // check whether 2 genes must be involved, or whether just one gene can explain the junction
        // NOTE: since not all chimeric reads may be available at this point, this test is repeated in the fusion routine
        if(spanGeneSpliceSites)
        {
            ++mChimericStats.LocalInterGeneFrags;
            return true;
        }

        // all reads within the gene - treat as alternative SJ candidates
        mLocalChimericReads.add(reads);
        return false;
    }

    private boolean skipNonGenicReads(final List<ReadRecord> reads)
    {
        // any set of entirely post-gene read(s) will be skipped and then picked up by the next gene collection's processing
        // otherwise record that they were processed to avoid double-processing them in the next gene collection
        List<ReadRecord> postGeneReads = !mGeneCollection.isEndOfChromosome() ? reads.stream()
                .filter(x -> x.PosStart > mGeneCollection.regionBounds()[SE_END])
                .collect(Collectors.toList()) : Lists.newArrayList();

        if(postGeneReads.size() == reads.size())
            return true;

        List<ReadRecord> preGeneReads = reads.stream()
                .filter(x -> x.PosStart < mGeneCollection.regionBounds()[SE_START])
                .collect(Collectors.toList());

        if(!preGeneReads.isEmpty())
        {
            // remove any previously processed reads
            final String readId = preGeneReads.get(0).Id;
            List<ReadRecord> prevPostGeneReads = mPreviousPostGeneReadMap.get(readId);

            if(prevPostGeneReads != null)
            {
                preGeneReads.stream().filter(x -> prevPostGeneReads.stream().anyMatch(y -> y.matches(x))).forEach(x -> reads.remove(x));

                if(reads.isEmpty())
                    return true;
            }
        }

        // cache and stop processing this group
        if(!postGeneReads.isEmpty())
            mPostGeneReadMap.put(reads.get(0).Id, postGeneReads);

        return false;
    }

    private void addRealignCandidates()
    {
        // in addition to the group having a least one read with the requires soft-clipping, the other read cannot extend past this
        // possible point of junction support
        Set<Integer>[] supportedJunctions = new Set[SE_PAIR];
        supportedJunctions[SE_START] = Sets.newHashSetWithExpectedSize(2); // from start boundaries, orientation -1
        supportedJunctions[SE_END] = Sets.newHashSetWithExpectedSize(2); // from end boundaries, orientation +1

        for(final ReadGroup readGroup : mCandidateRealignedReadMap.values())
        {
            supportedJunctions[SE_START].clear();
            supportedJunctions[SE_END].clear();

            for(ReadRecord read : readGroup.Reads)
            {
                for(int se = SE_START; se <= SE_END; ++se)
                {
                    final int seIndex = se;
                    if(!read.isSoftClipped(se))
                        continue;

                    int readBoundary = read.getCoordsBoundary(seIndex);

                    if(mJunctionPositions.stream().anyMatch(x -> positionWithin(readBoundary,
                            x - SOFT_CLIP_JUNC_BUFFER, x + SOFT_CLIP_JUNC_BUFFER)))
                    {
                        supportedJunctions[se].add(readBoundary);
                    }
                }
            }

            if(supportedJunctions[SE_START].isEmpty() && supportedJunctions[SE_END].isEmpty())
                continue;

            boolean validGroup = true;

            if(!supportedJunctions[SE_START].isEmpty()
            && readGroup.Reads.stream().anyMatch(x -> supportedJunctions[SE_START].stream().anyMatch(y -> x.PosStart < y - SOFT_CLIP_JUNC_BUFFER)))
            {
                validGroup = false;
            }
            else if(!supportedJunctions[SE_END].isEmpty()
            && readGroup.Reads.stream().anyMatch(x -> supportedJunctions[SE_END].stream().anyMatch(y -> x.PosEnd > y + SOFT_CLIP_JUNC_BUFFER)))
            {
                validGroup = false;
            }

            if(validGroup)
            {
                mChimericReadMap.put(readGroup.id(), readGroup);
                ++mChimericStats.CandidateRealignFrags;
            }
        }
    }

    private void collectCandidateJunctions(final ReadGroup readGroup)
    {
        if(readGroup.localJunctionData() != null)
        {
            mJunctionPositions.add(readGroup.localJunctionData().JunctionPositions[SE_START]);
            mJunctionPositions.add(readGroup.localJunctionData().JunctionPositions[SE_END]);
            return;
        }

        final int[] junctionPositions = findCandidateJunctions(readGroup.Reads);

        if(junctionPositions[SE_START] > 0)
            mJunctionPositions.add(junctionPositions[SE_START]);

        if(junctionPositions[SE_END] > 0)
            mJunctionPositions.add(junctionPositions[SE_END]);
    }

    private int[] findCandidateJunctions(final List<ReadRecord> reads)
    {
        final int[] junctionPositions = new int[SE_PAIR];

        final ReadRecord splitRead = findSplitRead(reads);

        if(splitRead != null)
        {
            return findSplitReadJunction(splitRead);
        }

        if(hasSuppAlignment(reads))
        {
            for(ReadRecord read : reads)
            {
                if(read.hasSuppAlignment())
                {
                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        if(hasRealignableSoftClip(read, se, false))
                            junctionPositions[se] = read.getCoordsBoundary(se);
                    }
                }
            }

            return junctionPositions;
        }

        // otherwise must either have a junction supported by 2 facing soft-clipped reads or a supplementary read
        // logic needs to match the type and junction assignment in FusionFragmentBuilder
        if(reads.size() == 1)
        {
            final ReadRecord read = reads.get(0);

            if(hasRealignableSoftClip(read, SE_START, false))
                junctionPositions[SE_START] = read.getCoordsBoundary(SE_START);

            if(hasRealignableSoftClip(read, SE_END, false))
                junctionPositions[SE_END] = read.getCoordsBoundary(SE_END);
        }
        else
        {
            int[] scPositions = {-1, -1};

            for(ReadRecord read : reads)
            {
                if(hasRealignableSoftClip(read, SE_START, false))
                    scPositions[SE_END] = read.getCoordsBoundary(SE_START);
                else if(hasRealignableSoftClip(read, SE_END, false))
                    scPositions[SE_START] = read.getCoordsBoundary(SE_END);
            }

            if(scPositions[SE_START] > 0 && scPositions[SE_END] > 0 && scPositions[SE_START] < scPositions[SE_END])
            {
                return scPositions;
            }
        }

        return junctionPositions;
    }
}
