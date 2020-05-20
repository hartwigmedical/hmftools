package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MAX_NOVEL_SJ_DISTANCE;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionFinder.addChimericReads;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentBuilder.hasSuppAlignment;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.findSplitReadJunction;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.hasRealignableSoftClip;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.isInversion;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;

public class ChimericReadTracker
{
    private final IsofoxConfig mConfig;

    private GeneCollection mGeneCollection; // the current collection being processed
    private final Map<String, List<ReadRecord>> mChimericReadMap;
    private final Set<Integer> mJunctionPositions; // from amongst the chimeric fragments with evidence of a fusion junction
    private final List<List<ReadRecord>> mLocalChimericReads; // fragments to re-evaluate as alternate splice sites
    private final Map<String,List<ReadRecord>> mCandidateRealignedReadMap;
    private final Map<String,Integer> mSecondaryReadMap;
    private final Set<String> mDuplicateReadIds; // currently only used to store chimeric duplicates

    // to avoid double-processing reads falling after a gene collection
    private final Map<String,List<ReadRecord>> mPostGeneReadMap;
    private final Map<String,List<ReadRecord>> mPreviousPostGeneReadMap;
    private final ChimericStats mChimericStats;

    public ChimericReadTracker(final IsofoxConfig config)
    {
        mConfig = config;
        mChimericStats = new ChimericStats();
        mChimericReadMap = Maps.newHashMap();
        mJunctionPositions = Sets.newHashSet();
        mDuplicateReadIds = Sets.newHashSet();
        mLocalChimericReads = Lists.newArrayList();
        mCandidateRealignedReadMap = Maps.newHashMap();
        mSecondaryReadMap = Maps.newHashMap();
        mPostGeneReadMap = Maps.newHashMap();
        mPreviousPostGeneReadMap = Maps.newHashMap();
        mGeneCollection = null;
    }

    public final Map<String,List<ReadRecord>> getReadMap() { return mChimericReadMap; }
    public final Set<Integer> getJunctionPositions() { return mJunctionPositions; }
    public final List<List<ReadRecord>> getLocalChimericReads() { return mLocalChimericReads; }
    public Set<String> getReadIds() { return mDuplicateReadIds; }
    public ChimericStats getStats() { return mChimericStats; }

    public void initialise(final GeneCollection geneCollection)
    {
        mGeneCollection = geneCollection;

        mPreviousPostGeneReadMap.clear();
        mPreviousPostGeneReadMap.putAll(mPostGeneReadMap);
        mPostGeneReadMap.clear();
    }

    public void clear()
    {
        mChimericReadMap.clear();
        mCandidateRealignedReadMap.clear();
        mChimericStats.clear();
        mSecondaryReadMap.clear();
        mLocalChimericReads.clear();
        mJunctionPositions.clear();
        mDuplicateReadIds.clear();
    }

    public void addRealignmentCandidates(final ReadRecord read1, final ReadRecord read2)
    {
        if(read1.isDuplicate() || read2.isDuplicate()) // group complete so drop these
            return;

        mCandidateRealignedReadMap.put(read1.Id, Lists.newArrayList(read1, read2));
    }

    public void registerSecondaryRead(final String readId)
    {
        // used for filtering fusions
        Integer count = mSecondaryReadMap.get(readId);
        if(count == null)
            mSecondaryReadMap.put(readId, 1);
        else
            mSecondaryReadMap.put(readId, count + 1);
    }

    public void addChimericReadPair(final ReadRecord read1, final ReadRecord read2)
    {
        if(mGeneCollection.inEnrichedRegion(read1.PosStart, read1.PosEnd) || mGeneCollection.inEnrichedRegion(read2.PosStart, read2.PosEnd))
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

            final List<ReadRecord> existingReads = mChimericReadMap.get(read1.Id);

            for(ReadRecord read : existingReads)
            {
                ISF_LOGGER.error("existing read: {}", read);
            }

            ISF_LOGGER.error("new read: {}", read1);
            ISF_LOGGER.error("new read: {}", read2);

            existingReads.add(read1);
            existingReads.add(read2);
        }
        else
        {
            mChimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));
        }
    }

    private void addIntronicTranscriptData(final ReadRecord read)
    {
        if(read.overlapsGeneCollection() && read.getMappedRegions().isEmpty())
            read.addIntronicTranscriptRefs(mGeneCollection.getTranscripts());
    }

    private static final String LOG_READ_ID = "";
    // private static final String LOG_READ_ID = "NB500901:18:HTYNHBGX2:3:21506:2131:4925";

    public void postProcessChimericReads(final BaseDepth baseDepth, final FragmentTracker fragmentTracker)
    {
        // check any lone reads - this cannot be one of a pair of non-genic reads since they will have already been dismissed
        // so will either be a supplementary or a read linked to another gene collection
        for(Object object : fragmentTracker.getValues())
        {
            final ReadRecord read = (ReadRecord)object;

            if(read.isMateUnmapped() || mGeneCollection.inEnrichedRegion(read.PosStart, read.PosEnd))
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

        for(List<ReadRecord> reads : mChimericReadMap.values())
        {
            // skip reads if all will be processed later or have been already
            final String readId = reads.get(0).Id;

            if(readId.equals(LOG_READ_ID))
            {
                ISF_LOGGER.debug("specific read: {}", readId);
            }

            int readCount = reads.size();
            boolean readGroupComplete = (readCount == 3 || (!hasSuppAlignment(reads) && readCount == 2));

            // if an entire group is a duplicate it can be dropped, otherwise record it's ID for the reads in other gene collections then drop it
            if(reads.stream().anyMatch(x -> x.isDuplicate()))
            {
                mDuplicateReadIds.add(readId);
                fragsToRemove.add(readId);
                continue;
            }

            if(skipNonGenicReads(reads))
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

            collectCandidateJunctions(reads, mJunctionPositions);
        }

        if(!fragsToRemove.isEmpty())
            fragsToRemove.forEach(x -> mChimericReadMap.remove(x));

        mChimericStats.ChimericJunctions += mJunctionPositions.size();

        int chimericCount = mChimericReadMap.size();
        mGeneCollection.addCount(TOTAL, chimericCount);
        mGeneCollection.addCount(CHIMERIC, chimericCount);

        for(List<ReadRecord> reads : mCandidateRealignedReadMap.values())
        {
            boolean addRead = false;

            for(ReadRecord read : reads)
            {
                for(int se = SE_START; se <= SE_END; ++se)
                {
                    final int seIndex = se;
                    if(read.isSoftClipped(se))
                    {
                        if(mJunctionPositions.stream().anyMatch(x -> positionWithin(read.getCoordsBoundary(seIndex),
                                x - SOFT_CLIP_JUNC_BUFFER, x + SOFT_CLIP_JUNC_BUFFER)))
                        {
                            addRead = true;
                            mChimericReadMap.put(read.Id, reads);
                            ++mChimericStats.CandidateRealignFrags;
                            break;
                        }
                    }
                }

                if(addRead)
                    break;
            }
        }

        // chimeric reads will be processed by the fusion-finding routine, so need to capture transcript and exon data
        // and free up other gene & region read data (to avoid retaining large numbers of references/memory)
        for(Map.Entry<String,List<ReadRecord>> entry : mChimericReadMap.entrySet())
        {
            final List<ReadRecord> reads = entry.getValue();

            Integer secondaryCount = mSecondaryReadMap.get(entry.getKey());

            reads.forEach(x -> x.captureGeneInfo(true));

            if(secondaryCount != null)
                reads.forEach(x -> x.setSecondaryReadCount(secondaryCount));
        }
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

        boolean spanGeneSpliceSites = hasMultipleKnownSpliceGenes(reads);

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

    private boolean hasMultipleKnownSpliceGenes(final List<ReadRecord> reads)
    {
        Set<Integer> junctionPositions = Sets.newHashSet();

        ReadRecord splitRead = null;

        for(ReadRecord read : reads)
        {
            if(read.containsSplit())
            {
                // find the largest N-split to mark the junction
                final int[] splitJunction = findSplitReadJunction(read);

                if(splitJunction != null)
                {
                    junctionPositions.add(splitJunction[SE_START]);
                    junctionPositions.add(splitJunction[SE_END]);
                    splitRead = read;
                }

                break;
            }

            if(hasRealignableSoftClip(read, SE_START, false))
                junctionPositions.add(read.getCoordsBoundary(SE_START));

            if(hasRealignableSoftClip(read, SE_END, false))
                junctionPositions.add(read.getCoordsBoundary(SE_END));
        }

        if(junctionPositions.size() < 2)
            return false;

        List<Set<String>> positionGeneLists = Lists.newArrayList();

        for(int juncPosition : junctionPositions)
        {
            Set<String> geneIds = Sets.newHashSet();

            if(splitRead != null)
            {
                splitRead.getMappedRegions().entrySet().stream()
                        .filter(x -> x.getKey().PosStart == juncPosition || x.getKey().PosEnd == juncPosition)
                        .forEach(x -> x.getKey().getTransExonRefs().stream().forEach(y -> geneIds.add(y.GeneId)));
            }
            else
            {
                reads.stream().forEach(x -> x.getMappedRegions().entrySet().stream()
                        .filter(y -> y.getKey().PosStart == juncPosition || y.getKey().PosEnd == juncPosition)
                        .forEach(y -> y.getKey().getTransExonRefs().stream().forEach(z -> geneIds.add(z.GeneId))));
            }

            positionGeneLists.add(geneIds);
        }

        for(int i = 0; i < positionGeneLists.size() - 1; ++i)
        {
            Set<String> geneIds1 = positionGeneLists.get(i);

            if(geneIds1.isEmpty())
                continue;

            for(int j = i + 1; j < positionGeneLists.size(); ++j)
            {
                Set<String> geneIds2 = positionGeneLists.get(j);

                if(geneIds2.isEmpty())
                    continue;

                if(!geneIds1.stream().anyMatch(x -> geneIds2.contains(x)))
                {
                    if(splitRead != null)
                        splitRead.setHasInterGeneSplit(); // will make use of this when handling fusions

                    return true;
                }
            }
        }

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

    private void collectCandidateJunctions(final List<ReadRecord> reads, final Set<Integer> junctionPositions)
    {
        for(ReadRecord read : reads)
        {
            if(read.hasInterGeneSplit() || (read.containsSplit() && read.spansGeneCollections()))
            {
                final int[] splitJunction = findSplitReadJunction(read);
                junctionPositions.add(splitJunction[SE_START]);
                junctionPositions.add(splitJunction[SE_END]);
                break;
            }

            if(hasRealignableSoftClip(read, SE_START, false))
                junctionPositions.add(read.getCoordsBoundary(SE_START));

            if(hasRealignableSoftClip(read, SE_END, false))
                junctionPositions.add(read.getCoordsBoundary(SE_END));
        }
    }



}
