package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_3;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.PROMISCUOUS_5;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_NOVEL_SJ_DISTANCE;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.fusion.ChimericUtils.findSplitReadJunction;
import static com.hartwig.hmftools.isofox.fusion.ChimericUtils.isInversion;
import static com.hartwig.hmftools.isofox.fusion.ChimericUtils.setHasMultipleKnownSpliceGenes;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionRead.convertReads;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
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

    private GeneCollection mGeneCollection; // the current collection being processed
    private final Map<String,ChimericReadGroup> mChimericReadMap;
    private final Map<String, FusionReadGroup> mFusionReadGroupMap;
    private final List<FusionReadGroup> mLocalCompleteGroups; // 2-read same-gene-collection groups with a split junction

    // junction position from fusion junction candidate reads are cached to identify candidate realignable reads
    private JunctionRacFragments mJunctionRacGroups;

    private final List<List<ReadRecord>> mLocalChimericReads; // fragments to re-evaluate as alternate splice sites
    private final List<ChimericReadGroup> mCandidateRealignedGroups;

    // map of candidate junction groups with supp data
    private final Set<String> mReferenceKnownGeneIds; // all as defined by the cache
    private final Set<String> mKnownGeneIds; // in the current gene collection
    private final List<String[]> mKnownPairGeneIds;

    private final Map<Integer,List<SupplementaryJunctionData>> mSupplementaryJunctions;
    private final Map<String,Set<String>> mHardFilteredReadIds;
    private Map<String,Set<Integer>> mKnownSpliteSites;

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
        mReferenceKnownGeneIds = Sets.newHashSet();
        mKnownGeneIds = Sets.newHashSet();
        mChimericStats = new ChimericStats();
        mChimericReadMap = Maps.newHashMap();
        mFusionReadGroupMap = Maps.newHashMap();
        mJunctionRacGroups = null;
        mLocalChimericReads = Lists.newArrayList();
        mLocalCompleteGroups = Lists.newArrayList();
        mCandidateRealignedGroups = Lists.newArrayList();
        mPostGeneReadMap = Maps.newHashMap();
        mPreviousPostGeneReadMap = Maps.newHashMap();
        mSupplementaryJunctions = Maps.newHashMap();
        mHardFilteredReadIds = Maps.newHashMap();
        mGeneCollection = null;
        mKnownSpliteSites = null;
    }

    public boolean enabled() { return mEnabled; }

    public Map<String, FusionReadGroup> getReadMap() { return mFusionReadGroupMap; }

    public JunctionRacFragments extractJunctionRacFragments()
    {
        mJunctionRacGroups.purgeEmptyJunctions();
        mChimericStats.ChimericJunctions += mJunctionRacGroups.junctionCount();
        mChimericStats.CandidateRealignFrags += mJunctionRacGroups.fragmentCount();
        return mJunctionRacGroups;
    }

    public void setKnownSpliteSites(final Map<String,Set<Integer>> knownSpliceites) { mKnownSpliteSites = knownSpliceites; }

    public JunctionRacFragments getJunctionRacGroups() { return mJunctionRacGroups; } // testing only
    public List<List<ReadRecord>> getLocalChimericReads() { return mLocalChimericReads; }
    public ChimericStats getStats() { return mChimericStats; }
    public Map<String,Set<String>> getHardFilteredReadIds() { return mHardFilteredReadIds; }

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
        List<KnownFusionType> knownTypes = Lists.newArrayList(KNOWN_PAIR, PROMISCUOUS_5, PROMISCUOUS_3);

        for(KnownFusionType type : knownTypes)
        {
            List<KnownFusionData> typeData = mConfig.Fusions.KnownFusions.getDataByType(type);

            if(typeData == null)
                continue;

            for(final KnownFusionData knownData : typeData)
            {
                final GeneData upGene = !knownData.FiveGene.isEmpty() ? geneTransCache.getGeneDataByName(knownData.FiveGene) : null;
                final GeneData downGene = !knownData.ThreeGene.isEmpty() ? geneTransCache.getGeneDataByName(knownData.ThreeGene) : null;

                if(type == KNOWN_PAIR && upGene != null && downGene != null)
                {
                    mKnownPairGeneIds.add(new String[] { upGene.GeneId, downGene.GeneId });
                }

                if(upGene != null)
                    mReferenceKnownGeneIds.add(upGene.GeneId);

                if(downGene != null)
                    mReferenceKnownGeneIds.add(downGene.GeneId);
            }
        }
    }

    public void initialise(final GeneCollection geneCollection)
    {
        mGeneCollection = geneCollection;
        mKnownGeneIds.clear();
        mGeneCollection.geneIds().stream().filter(x -> mReferenceKnownGeneIds.contains(x)).forEach(x -> mKnownGeneIds.add(x));
        mJunctionRacGroups = new JunctionRacFragments();

        mPreviousPostGeneReadMap.clear();
        mPreviousPostGeneReadMap.putAll(mPostGeneReadMap);
        mPostGeneReadMap.clear();
    }

    public void clear() { clear(false); }
    public void clearAll() { clear(true); }

    private void clear(boolean full)
    {
        mChimericReadMap.clear();
        mFusionReadGroupMap.clear();
        mLocalCompleteGroups.clear();
        mCandidateRealignedGroups.clear();
        mChimericStats.clear();
        mLocalChimericReads.clear();
        mSupplementaryJunctions.clear();

        if(full)
        {
            mPreviousPostGeneReadMap.clear();
            mPostGeneReadMap.clear();
            mJunctionRacGroups = null;
            mHardFilteredReadIds.clear();
        }
    }

    public void addRealignmentCandidates(final ReadRecord read1, final ReadRecord read2)
    {
        if(read1.isDuplicate() || read2.isDuplicate()) // group complete so drop these
            return;

        mCandidateRealignedGroups.add(new ChimericReadGroup(read1, read2));
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

            final ChimericReadGroup existingGroup = mChimericReadMap.get(read1.Id);

            for(ReadRecord read : existingGroup.reads())
            {
                ISF_LOGGER.error("existing read: {}", read);
            }

            ISF_LOGGER.error("new read: {}", read1);
            ISF_LOGGER.error("new read: {}", read2);

            existingGroup.addRead(read1);
            existingGroup.addRead(read2);
        }
        else
        {
            mChimericReadMap.put(read1.Id, new ChimericReadGroup(read1, read2));
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

            ChimericReadGroup chimericReads = mChimericReadMap.get(read.Id);
            if (chimericReads == null)
                mChimericReadMap.put(read.Id, new ChimericReadGroup(read));
            else
                chimericReads.addRead(read);
        }

        // migrate any local chimeric fragments for analysis as alternate splice junctions
        final List<String> fragsToRemove = Lists.newArrayList();

        for(final ChimericReadGroup readGroup : mChimericReadMap.values())
        {
            // skip reads if all will be processed later or have been already
            final List<ReadRecord> reads = readGroup.reads();
            final String readId = reads.get(0).Id;

            int readCount = reads.size();
            boolean readGroupComplete = readGroup.isComplete();

            // duplicates are kept until a group is complete, since not all reads are marked as duplicates and those would otherwise
            // be orphaned later on

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
            {
                collectCandidateJunctions(readGroup);
                cacheSupplementaryJunctionCandidate(readGroup);
            }
        }

        if(!fragsToRemove.isEmpty())
            fragsToRemove.forEach(x -> mChimericReadMap.remove(x));

        int chimericCount = mChimericReadMap.size();
        mGeneCollection.addCount(CHIMERIC, chimericCount);

        applyHardFilter();

        if(mRunFusions)
        {
            addRealignCandidates();

            // chimeric reads will be processed by the fusion-finding routine, so need to capture transcript and exon data
            // and free up other gene & region read data (to avoid retaining large numbers of references/memory)
            for(final ChimericReadGroup readGroup : mChimericReadMap.values())
            {
                List<FusionRead> reads = convertReads(readGroup.reads());
                reads.forEach(x -> x.setReadJunctionDepth(baseDepth));
                mFusionReadGroupMap.put(readGroup.id(), new FusionReadGroup(readGroup.id(), reads));
            }

            for(FusionFragment fragment : mJunctionRacGroups.getRacFragments())
            {
                fragment.reads().forEach(x -> x.setReadJunctionDepth(baseDepth));
            }
        }
        else
        {
            // clear other chimeric state except for local junction information
            mChimericReadMap.clear();
            mLocalCompleteGroups.clear();
            mCandidateRealignedGroups.clear();
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
            SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read.getSuppAlignment());

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
        // cache each RAC against the junction(s) it may support, and discard any which support none
        // RACs are no longer added to the chimeric read map since they now only exist in here
        mCandidateRealignedGroups.forEach(x -> mJunctionRacGroups.checkAddCandidateGroup(x));
    }

    private void collectCandidateJunctions(final ChimericReadGroup readGroup)
    {
        // type 1: split reads
        final ReadRecord splitRead = readGroup.reads().stream()
            .filter(x -> x.containsSplit())
            .filter(x -> x.spansGeneCollections() || x.hasInterGeneSplit())
            .findFirst().orElse(null);

        if(splitRead != null)
        {
            int[] splitJunction = findSplitReadJunction(splitRead);

            addJunction(splitRead, SE_START, splitJunction[SE_START], POS_ORIENT);
            addJunction(splitRead, SE_END, splitJunction[SE_END], NEG_ORIENT);

            return;
        }

        final int[] junctionPositions = new int[SE_PAIR];

        // type 2: supplementary with clipping
        final ReadRecord suppRead = readGroup.reads().stream().filter(x -> x.hasSuppAlignment()).findFirst().orElse(null);
        if(suppRead != null)
        {
            ClippedSide scSide = ReadRecord.clippedSide(suppRead);

            if(scSide != null && scSide.Length >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH)
            {
                junctionPositions[scSide.Side] = suppRead.getCoordsBoundary(scSide.Side);

                addJunction(suppRead, scSide.Side, junctionPositions[scSide.Side], scSide.Side == SE_START ? NEG_ORIENT : POS_ORIENT);
            }

            return;
        }

        // type 3: discordant with at least 1 known-pair gene, must be a single-read group since the other will be distant
        if(readGroup.reads().size() > 1 || mKnownGeneIds.isEmpty() || !readInKnownGene(readGroup))
            return;

        // select the side with the longest soft-clipping
        ReadRecord read = readGroup.reads().get(0);

        ClippedSide scSide = ReadRecord.clippedSide(read);

        if(scSide != null)
        {
            if(scSide.isLeft())
            {
                junctionPositions[SE_START] = read.getCoordsBoundary(SE_START);
                addJunction(read, SE_START, junctionPositions[SE_START], NEG_ORIENT);
            }
            else
            {
                junctionPositions[SE_END] = read.getCoordsBoundary(SE_END);
                addJunction(read, SE_END, junctionPositions[SE_END], POS_ORIENT);
            }
        }
    }

    private void addJunction(final ReadRecord read, int seIndex, int juncPosition, byte juncOrientation)
    {
        if(juncPosition > 0)
        {
            read.setJunctionPosition(seIndex, juncPosition);
            mJunctionRacGroups.addJunction(juncPosition, juncOrientation);
        }
    }

    private boolean readInKnownGene(final ChimericReadGroup readGroup)
    {
        if(mKnownGeneIds.isEmpty())
            return false;

        for(ReadRecord read : readGroup.reads())
        {
            if(read.getMappedRegions().keySet().stream()
                    .anyMatch(x -> x.getTransExonRefs().stream().anyMatch(y -> mKnownGeneIds.contains(y.GeneId))))
            {
                return true;
            }
        }

        return false;
    }

    private void applyHardFilter()
    {
        if(mConfig.Fusions.MinHardFilterFrags <= 1 || !mKnownGeneIds.isEmpty())
            return;

        HardFilteredCache.applyHardFilter(
                mSupplementaryJunctions, mHardFilteredReadIds, mChimericReadMap, mGeneCollection.chromosome(),
                mConfig.Fusions.MinHardFilterFrags, mKnownSpliteSites);

        mSupplementaryJunctions.clear();
    }

    private void cacheSupplementaryJunctionCandidate(final ChimericReadGroup readGroup)
    {
        if(mConfig.Fusions.MinHardFilterFrags <= 1 || readInKnownGene(readGroup))
            return;

        SupplementaryJunctionData.cacheSupplementaryJunctionCandidate(readGroup, mSupplementaryJunctions);
    }
}
