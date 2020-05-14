package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.GeneCollection.NON_GENIC_ID;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.common.RnaUtils.canonicalAcceptor;
import static com.hartwig.hmftools.isofox.common.RnaUtils.canonicalDonor;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.JUNCTION_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formChromosomePair;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocationPair;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class FusionFragment
{
    private final List<ReadRecord> mReads;
    private final List<ReadRecord>[] mReadsByLocation;
    private final boolean mHasSupplementaryAlignment;

    private final int[] mGeneCollections;
    private final boolean[] mInGenicRegions;
    private final String[] mChromosomes;
    private final byte[] mOrientations;
    private final int[] mJunctionPositions; // fusion junction if exists
    private final byte[] mJunctionOrientations; // orientation at junction if exists
    private final FusionJunctionType[] mJunctionTypes;
    private final String[] mJunctionBases; // the 10 bases leading up to the junction if it exists
    private final String[] mJunctionSpliceBases; // the 2 donor/acceptor bases
    private FusionFragmentType mType;
    private final String[] mLocationIds; // used to group proximate fragments and fusions

    private final RegionMatchType[] mRegionMatchTypes; // top-ranking region match type from the reads
    private final List<TransExonRef>[] mTransExonRefs;

    public FusionFragment(final List<ReadRecord> reads)
    {
        mReads = reads;
        mHasSupplementaryAlignment = FusionFragmentBuilder.hasSuppAlignment(mReads);

        mReadsByLocation = new List[SE_PAIR];
        mReadsByLocation[SE_START] = Lists.newArrayList();
        mReadsByLocation[SE_END] = Lists.newArrayList();

        mGeneCollections = new int[] {NON_GENIC_ID, NON_GENIC_ID};
        mInGenicRegions = new boolean[] {true, true};
        mJunctionPositions = new int[] {-1, -1};
        mChromosomes = new String[] {"", ""};
        mJunctionOrientations = new byte[] {0, 0};
        mOrientations = new byte[] {0, 0};
        mRegionMatchTypes = new RegionMatchType[] { RegionMatchType.NONE, RegionMatchType.NONE };
        mJunctionTypes = new FusionJunctionType[] { FusionJunctionType.UNKNOWN, FusionJunctionType.UNKNOWN };
        mJunctionBases = new String[] {"", ""};
        mJunctionSpliceBases = new String[] {"", ""};

        mTransExonRefs = new List[SE_PAIR];
        mTransExonRefs[SE_START] = Lists.newArrayList();
        mTransExonRefs[SE_END] = Lists.newArrayList();

        mType = UNKNOWN;

        FusionFragmentBuilder.setLocationData(this);
        FusionFragmentBuilder.setJunctionData(this);
        mLocationIds = FusionFragmentBuilder.setLocationIds(this);

        extractTranscriptExonData();
    }

    public String readId() { return mReads.get(0).Id; }

    public final List<ReadRecord> readsByLocation(int se) { return mReadsByLocation[se]; }
    public final List<ReadRecord> reads() { return mReads; }

    public FusionFragmentType type() { return mType; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final boolean[] inGenicRegions() { return mInGenicRegions; }
    public final byte[] orientations() { return mOrientations; }
    public final String[] locationIds() { return mLocationIds; }

    public boolean hasSuppAlignment() { return mHasSupplementaryAlignment; }

    public boolean isSingleGene()
    {
        return mChromosomes[SE_START].equals(mChromosomes[SE_END]) && mGeneCollections[SE_START] == mGeneCollections[SE_END]
                && mInGenicRegions[SE_START] == mInGenicRegions[SE_END];
    }

    public final int[] junctionPositions() { return mJunctionPositions; }
    public final byte[] junctionOrientations() { return mJunctionOrientations; }
    public final String[] junctionBases() { return mJunctionBases; }

    public final RegionMatchType[] regionMatchTypes() { return mRegionMatchTypes; }
    public final FusionJunctionType[] junctionTypes() { return mJunctionTypes; }

    public boolean isUnspliced() { return mRegionMatchTypes[SE_START] == INTRON && mRegionMatchTypes[SE_END] == INTRON; }
    public boolean isSpliced() { return exonBoundary(mRegionMatchTypes[SE_START]) && exonBoundary(mRegionMatchTypes[SE_END]); }

    public String locationPair()
    {
        if(isSingleGene())
            return mLocationIds[SE_START];
        else
            return String.format("%s_%s", mLocationIds[SE_START], mLocationIds[SE_END]);
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

    private void extractTranscriptExonData()
    {
        // set transcript & exon info for each junction from each applicable read
        // only take the highest matches
        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<ReadRecord> reads = isSingleGene() ? mReads : readsByLocation(se);

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
                    .filter(x -> (junctionOrientation == 1 && x.ExonEnd == junctionPosition)
                            || (junctionOrientation == -1 && x.ExonStart == junctionPosition))
                    .findFirst().orElse(null);

            if(exon == null || transExonRef.ExonRank != exon.ExonRank)
            {
                mTransExonRefs[seIndex].remove(index);
                continue;
            }

            mJunctionTypes[seIndex] = KNOWN;
            ++index;
        }
    }

    public void setJunctionTypes(final IndexedFastaSequenceFile refGenome, final byte[] junctionStrands)
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
                String daBases = mJunctionSpliceBases[se];

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

    public void setJunctionBases(final IndexedFastaSequenceFile refGenome)
    {
        if(mType != MATCHED_JUNCTION)
            return;

        if(refGenome != null)
        {
            try
            {
                for (int se = SE_START; se <= SE_END; ++se)
                {
                    int junctionBase = mJunctionPositions[se];

                    if (junctionOrientations()[se] == 1)
                    {
                        String junctionBases = refGenome.getSubsequenceAt(
                                mChromosomes[se], junctionBase - JUNCTION_BASE_LENGTH + 1, junctionBase + 2).getBaseString();

                        mJunctionBases[se] = junctionBases.substring(0, JUNCTION_BASE_LENGTH);
                        mJunctionSpliceBases[se] = junctionBases.substring(JUNCTION_BASE_LENGTH);
                    }
                    else
                    {
                        String junctionBases = refGenome.getSubsequenceAt(
                                mChromosomes[se], junctionBase - 2, junctionBase + JUNCTION_BASE_LENGTH - 1).getBaseString();

                        mJunctionBases[se] = junctionBases.substring(2);
                        mJunctionSpliceBases[se] = junctionBases.substring(0, 2);
                    }
                }
            }
            catch(Exception e)
            {
                // junction may be in an invalid region, just ignore these
            }
        }
        else
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                final int seIndex = se;
                int junctionBase = mJunctionPositions[se];

                ReadRecord read = mReads.stream()
                        .filter(x -> x.Chromosome.equals(mChromosomes[seIndex]))
                        .filter(x -> x.getMappedRegionCoords().stream().anyMatch(y -> positionWithin(junctionBase, y[SE_START], y[SE_END])))
                        // .filter(x -> x.getCoordsBoundary(junctionOrientations()[seIndex] == 1 ? SE_END : SE_START) == junctionBase)
                        .findFirst().orElse(null);

                if(read.getCoordsBoundary(junctionOrientations()[seIndex] == 1 ? SE_END : SE_START) == junctionBase)
                {
                    if (junctionOrientations()[se] == 1)
                    {
                        mJunctionBases[se] = read.ReadBases.substring(read.Length - JUNCTION_BASE_LENGTH, read.Length);
                    }
                    else
                    {
                        mJunctionBases[se] = read.ReadBases.substring(0, JUNCTION_BASE_LENGTH);
                    }
                }
                else
                {
                    // find the read bases around the junction manually
                    int baseIndex = 0;
                    for(int[] coord : read.getMappedRegionCoords())
                    {
                        if(junctionBase > coord[SE_END])
                        {
                            baseIndex += coord[SE_END] - coord[SE_START] + 1;
                        }
                        else
                        {
                            baseIndex += junctionBase - coord[SE_START];
                            break;
                        }
                    }

                    if (junctionOrientations()[se] == 1)
                    {
                        ++baseIndex;
                        mJunctionBases[se] = read.ReadBases.substring(baseIndex - JUNCTION_BASE_LENGTH, baseIndex);
                    }
                    else
                    {
                        mJunctionBases[se] = read.ReadBases.substring(baseIndex, baseIndex + JUNCTION_BASE_LENGTH);
                    }
                }
            }
        }
    }

}
