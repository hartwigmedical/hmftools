package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionFinder.addChimericReads;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentBuilder.hasSuppAlignment;

import java.util.List;
import java.util.Map;
import java.util.Set;

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

    private final Map<String, List<ReadRecord>> mChimericReadMap;
    private final Map<String,List<ReadRecord>> mCandidateRealignedReadMap;
    private final Map<String,Integer> mSecondaryReadMap;
    private GeneCollection mGeneCollection;
    private final ChimericStats mChimericStats;

    public ChimericReadTracker(final IsofoxConfig config)
    {
        mConfig = config;
        mChimericStats = new ChimericStats();
        mChimericReadMap = Maps.newHashMap();
        mCandidateRealignedReadMap = Maps.newHashMap();
        mSecondaryReadMap = Maps.newHashMap();
        mGeneCollection = null;
    }

    public final Map<String,List<ReadRecord>> getReadMap() { return mChimericReadMap; }
    public ChimericStats getStats() { return mChimericStats; }

    public void initialise(final GeneCollection geneCollection)
    {
        mGeneCollection = geneCollection;
    }

    public void clear()
    {
        mChimericReadMap.clear();
        mCandidateRealignedReadMap.clear();
        mChimericStats.clear();
        mSecondaryReadMap.clear();
    }

    public static boolean isRealignedFragmentCandidate(final ReadRecord read)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if (read.isSoftClipped(se))
            {
                int scLength =
                        se == SE_START ? read.Cigar.getFirstCigarElement().getLength() : read.Cigar.getLastCigarElement().getLength();

                if (scLength >= REALIGN_MIN_SOFT_CLIP_BASE_LENGTH && scLength <= REALIGN_MAX_SOFT_CLIP_BASE_LENGTH)
                    return true;
            }
        }

        return false;
    }

    public void addRealignmentCandidates(final ReadRecord read1, final ReadRecord read2)
    {
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

        // populate transcript info for intronic reads since it will be used in fusion matching
        addIntronicTranscriptData(read1);
        addIntronicTranscriptData(read2);

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

    public void postProcessChimericReads(final BaseDepth baseDepth, final FragmentTracker fragmentTracker)
    {
        // check any lone reads - this cannot be one of a pair of non-genic reads since they will have already been dismissed
        // so will either be a supplementary or a read linked to another gene collection
        for(Object object : fragmentTracker.getValues())
        {
            final ReadRecord read = (ReadRecord)object;

            if(read.isMateUnmapped() || mGeneCollection.inEnrichedRegion(read.PosStart, read.PosEnd))
                continue;

            baseDepth.processRead(read.getMappedRegionCoords());
            addIntronicTranscriptData(read);
            addChimericReads(mChimericReadMap, read);
        }

        // find any split chimeric reads and use this to select from the candidates for realignment
        final Set<Integer> chimericJunctions = Sets.newHashSet();
        for(List<ReadRecord> reads : mChimericReadMap.values())
        {
            for(ReadRecord read : reads)
            {
                if(read.getSuppAlignment() == null)
                    continue;

                int scLeft = read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0;
                int scRight = read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0;

                int junctionSeIndex = scLeft > scRight ? SE_START : SE_END;
                chimericJunctions.add(read.getCoordsBoundary(junctionSeIndex));
            }
        }

        mChimericStats.ChimericJunctions += chimericJunctions.size();

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
                        if(chimericJunctions.stream().anyMatch(x -> positionWithin(read.getCoordsBoundary(seIndex),
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
        List<String> nonGenicFrags = Lists.newArrayList();
        for(Map.Entry<String,List<ReadRecord>> entry : mChimericReadMap.entrySet())
        {
            final List<ReadRecord> reads = entry.getValue();

            Integer secondaryCount = mSecondaryReadMap.get(entry.getKey());

            reads.forEach(x -> x.captureGeneInfo());

            if(secondaryCount != null)
                reads.forEach(x -> x.setSecondaryReadCount(secondaryCount));

            boolean hasSuppData = hasSuppAlignment(reads);
            boolean isTranslocation = reads.stream().anyMatch(x -> x.isTranslocation());
            boolean spanGenes = reads.stream().anyMatch(x -> x.spansGeneCollections());
            boolean isInversion = reads.stream().anyMatch(x -> x.isInversion());

            if(!isTranslocation && !spanGenes && reads.stream().anyMatch(x -> x.preGeneCollection()))
            {
                int preGeneReads = (int)reads.stream().filter(x -> x.preGeneCollection()).count();

                if(preGeneReads == reads.size() && (reads.size() == 3 || reads.size() == 2 && !hasSuppData))
                {
                    nonGenicFrags.add(entry.getKey());
                    continue;
                }

                ++mChimericStats.LocalPreGeneFrags;
            }

            if(hasSuppData)
                ++mChimericStats.SupplementaryFrags;

            if(isTranslocation)
                ++mChimericStats.Translocations;
            else if(isInversion)
                ++mChimericStats.Inversions;
            else if(spanGenes)
                ++mChimericStats.LocalInterGeneFrags;

        }

        if(!nonGenicFrags.isEmpty())
            nonGenicFrags.forEach(x -> mChimericReadMap.remove(x));
    }


}
