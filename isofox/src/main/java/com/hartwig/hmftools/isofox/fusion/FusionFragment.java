package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.ReadRecord.NO_GENE_ID;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.common.RnaUtils.canonicalAcceptor;
import static com.hartwig.hmftools.isofox.common.RnaUtils.canonicalDonor;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionFragment
{
    private final ReadGroup mReadGroup;
    private final boolean mHasSupplementaryAlignment;

    private final int[] mGeneCollections;
    private final String[] mChromosomes;
    private final byte[] mOrientations;
    private final int[] mJunctionPositions; // fusion junction if exists
    private final byte[] mJunctionOrientations; // orientation at junction if exists
    private final FusionJunctionType[] mJunctionTypes;
    private FusionFragmentType mType;
    private final String[] mLocationIds; // used to group proximate fragments and fusions

    private final RegionMatchType[] mRegionMatchTypes; // top-ranking region match type from the reads
    private final List<TransExonRef>[] mTransExonRefs;

    private final Set<FusionReadData> mFusions;

    public FusionFragment(final ReadGroup readGroup)
    {
        mReadGroup = readGroup;
        mHasSupplementaryAlignment = mReadGroup.hasSuppAlignment();

        mGeneCollections = new int[] {NO_GENE_ID, NO_GENE_ID};
        mJunctionPositions = new int[] {-1, -1};
        mChromosomes = new String[] {"", ""};
        mJunctionOrientations = new byte[] {0, 0};
        mOrientations = new byte[] {0, 0};
        mRegionMatchTypes = new RegionMatchType[] { RegionMatchType.NONE, RegionMatchType.NONE };
        mJunctionTypes = new FusionJunctionType[] { FusionJunctionType.UNKNOWN, FusionJunctionType.UNKNOWN };
        mLocationIds = new String[] {"", ""};

        mTransExonRefs = new List[SE_PAIR];
        mTransExonRefs[SE_START] = Lists.newArrayList();
        mTransExonRefs[SE_END] = Lists.newArrayList();

        mFusions = Sets.newHashSet();

        mType = UNKNOWN;

        FusionFragmentBuilder.setFragmentProperties(this);

        if(mType != UNKNOWN)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                mLocationIds[se] = formLocation(mChromosomes[se], mGeneCollections[se], true);
            }
        }

        extractTranscriptExonData();
    }

    public static FusionFragment from(final List<ReadRecord> reads) { return new FusionFragment(new ReadGroup(reads)); }

    public String readId() { return mReadGroup.id(); }

    public final List<ReadRecord> reads() { return mReadGroup.Reads; }
    public final ReadGroup readGroup() { return mReadGroup; }

    public FusionFragmentType type() { return mType; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final byte[] orientations() { return mOrientations; }
    public ChrGeneCollectionPair chrGeneCollection(int se) { return new ChrGeneCollectionPair(mChromosomes[se], mGeneCollections[se]); }
    public final Set<FusionReadData> assignedFusions() { return mFusions; }

    public boolean hasSuppAlignment() { return mHasSupplementaryAlignment; }

    public boolean isSingleChromosome() { return mChromosomes[SE_START].equals(mChromosomes[SE_END]); }

    public boolean isSingleGeneCollection()
    {
        return mChromosomes[SE_START].equals(mChromosomes[SE_END]) && mGeneCollections[SE_START] == mGeneCollections[SE_END];
    }

    public final int[] junctionPositions() { return mJunctionPositions; }
    public final byte[] junctionOrientations() { return mJunctionOrientations; }

    public final RegionMatchType[] regionMatchTypes() { return mRegionMatchTypes; }
    public final FusionJunctionType[] junctionTypes() { return mJunctionTypes; }

    public boolean isUnspliced() { return mRegionMatchTypes[SE_START] == INTRON && mRegionMatchTypes[SE_END] == INTRON; }
    public boolean isSpliced() { return exonBoundary(mRegionMatchTypes[SE_START]) && exonBoundary(mRegionMatchTypes[SE_END]); }

    public String locationPair()
    {
        if(isSingleGeneCollection())
            return mLocationIds[SE_START];
        else
            return String.format("%s_%s", mLocationIds[SE_START], mLocationIds[SE_END]);
    }

    public String positionHash()
    {
        return String.format("%d_%d_%d_%d",
                mJunctionPositions[SE_START], mJunctionOrientations[SE_START],
                mJunctionPositions[SE_END], mJunctionOrientations[SE_END]);
    }

    public void setType(FusionFragmentType type) { mType = type; }

    public StructuralVariantType getImpliedSvType()
    {
        return impliedSvType(mChromosomes, mJunctionOrientations);
    }

    public final List<TransExonRef>[] getTransExonRefs() { return mTransExonRefs; }

    public List<String> getGeneIds(int seIndex)
    {
        final List<String> geneIds = Lists.newArrayList();

        for(TransExonRef transExonRef : mTransExonRefs[seIndex])
        {
            if(!geneIds.contains(transExonRef.GeneId))
                geneIds.add(transExonRef.GeneId);
        }

        return geneIds;
    }

    public final List<ReadRecord> readsByLocation(final int se)
    {
        if(mType == UNKNOWN)
            return Lists.newArrayList();

        // doesn't take junctions into consideration
        if(isSingleGeneCollection())
            return mReadGroup.Reads;

        return mReadGroup.Reads.stream()
                .filter(x -> x.Chromosome.equals(mChromosomes[se]))
                .filter(x -> x.getGeneCollectons()[SE_START] == mGeneCollections[se] || x.getGeneCollectons()[SE_END] == mGeneCollections[se])
                .collect(Collectors.toList());
    }

    private void extractTranscriptExonData()
    {
        // set transcript & exon info for each junction from each applicable read
        // only take the highest matches
        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<ReadRecord> reads = readsByLocation(se);

            for (final ReadRecord read : reads)
            {
                if(mJunctionPositions[se] > 0 && !positionWithin(mJunctionPositions[se], read.PosStart, read.PosEnd))
                    continue;

                final Map<RegionMatchType,List<TransExonRef>> transExonRefMap = read.getTransExonRefs(se);

                for (Map.Entry<RegionMatchType, List<TransExonRef>> entry : transExonRefMap.entrySet())
                {
                    RegionMatchType matchType = entry.getKey();

                    if(matchRank(matchType) >= matchRank(mRegionMatchTypes[se]))
                        mRegionMatchTypes[se] = matchType;
                    else
                        continue;

                    for(TransExonRef readTransExonRef : entry.getValue())
                    {
                        boolean found = false;

                        for(TransExonRef transExonRef : mTransExonRefs[se])
                        {
                            if (transExonRef.TransId == readTransExonRef.TransId)
                            {
                                found = true;

                                if (transExonRef.ExonRank != readTransExonRef.ExonRank)
                                {
                                    ISF_LOGGER.trace("multi-exon: read({} cigar={}) ref1({}) ref2({})",
                                            read.Id, read.Cigar.toString(), transExonRef.toString(), readTransExonRef.toString());

                                    // will be handled later on
                                    mTransExonRefs[se].add(readTransExonRef);
                                }

                                break;
                            }
                        }

                        if(!found)
                        {
                            mTransExonRefs[se].add(readTransExonRef);
                        }
                    }
                }
            }
        }
    }

    public void validateTranscriptExons(final List<TranscriptData> transDataList, int seIndex)
    {
        if(mJunctionPositions[seIndex] < 0 || !exonBoundary(mRegionMatchTypes[seIndex]))
            return;

        int junctionPosition = mJunctionPositions[seIndex];
        byte junctionOrientation = mJunctionOrientations[seIndex];

        int index = 0;
        while(index < mTransExonRefs[seIndex].size())
        {
            TransExonRef transExonRef = mTransExonRefs[seIndex].get(index);
            TranscriptData transData = transDataList.stream().filter(x -> x.TransId == transExonRef.TransId).findFirst().orElse(null);

            if(transData == null)
            {
                mTransExonRefs[seIndex].remove(index);
                continue;
            }

            ExonData exon = transData.exons().stream()
                    .filter(x -> (junctionOrientation == 1 && x.End == junctionPosition)
                            || (junctionOrientation == -1 && x.Start == junctionPosition))
                    .findFirst().orElse(null);

            if(exon == null || transExonRef.ExonRank != exon.Rank)
            {
                mTransExonRefs[seIndex].remove(index);
                continue;
            }

            mJunctionTypes[seIndex] = KNOWN;
            ++index;
        }
    }

    public void setJunctionTypes(final RefGenomeInterface refGenome, final byte[] junctionStrands, final String[] junctionSpliceBases)
    {
        if(refGenome == null)
            return;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if (mJunctionTypes[se] == KNOWN)
            {
                continue;
            }
            else if(mJunctionPositions[se] > 0)
            {
                String daBases = junctionSpliceBases[se];

                if(junctionStrands != null)
                {
                    boolean isDonor = (mJunctionOrientations[se] == junctionStrands[se]);

                    if (isDonor && canonicalDonor(daBases, junctionStrands[se]))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                    else if (!isDonor && canonicalAcceptor(daBases, junctionStrands[se]))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                }
                else
                {
                    // try them both
                    byte asDonorStrand = mJunctionOrientations[se];
                    byte asAcceptorStrand = (byte)(-mJunctionOrientations[se]);

                    if (canonicalDonor(daBases, asDonorStrand) || canonicalAcceptor(daBases, asAcceptorStrand))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                }
            }
        }
    }

    public String toString()
    {
        return String.format("type(%s) chr(%s-%s) junc(%d-%d %d/%d %s)",
                mType, mChromosomes[SE_START], mChromosomes[SE_END], mJunctionPositions[SE_START], mJunctionPositions[SE_END],
                mJunctionOrientations[SE_START], mJunctionOrientations[SE_END], getImpliedSvType());
    }
}
