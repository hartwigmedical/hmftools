package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.common.ReadRecord.NO_GENE_ID;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.exonBoundary;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.matchRank;
import static com.hartwig.hmftools.isofox.common.CommonUtils.canonicalAcceptor;
import static com.hartwig.hmftools.isofox.common.CommonUtils.canonicalDonor;
import static com.hartwig.hmftools.isofox.common.CommonUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.UNKNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionTransExon.mergeUnique;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocation;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.RegionMatchType;

public class FusionFragment
{
    private final FusionReadGroup mReadGroup;
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
    private final List<FusionTransExon>[] mTransExonRefs;

    public FusionFragment(final FusionReadGroup readGroup)
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

    public String readId() { return mReadGroup.ReadId; }
    public List<FusionRead> reads() { return mReadGroup.Reads; }
    public FusionReadGroup readGroup() { return mReadGroup; }

    public FusionFragmentType type() { return mType; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final byte[] orientations() { return mOrientations; }

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

    public final List<FusionTransExon>[] getTransExonRefs() { return mTransExonRefs; }

    /*
    public List<String> getGeneIds(int seIndex)
    {
        final List<String> geneIds = Lists.newArrayList();

        for(FusionTransExonRef transExonRef : mTransExonRefs[seIndex])
        {
            if(!geneIds.contains(transExonRef.GeneId))
                geneIds.add(transExonRef.GeneId);
        }

        return geneIds;
    }
    */

    public final List<FusionRead> readsByLocation(final int se)
    {
        if(mType == UNKNOWN)
            return Lists.newArrayList();

        // doesn't take junctions into consideration
        if(isSingleGeneCollection())
            return mReadGroup.Reads;

        return mReadGroup.Reads.stream()
                .filter(x -> x.Chromosome.equals(mChromosomes[se]))
                .filter(x -> x.GeneCollections[SE_START] == mGeneCollections[se] || x.GeneCollections[SE_END] == mGeneCollections[se])
                .collect(Collectors.toList());
    }

    private void extractTranscriptExonData()
    {
        // set transcript & exon info for each junction from each applicable read, taking only the highest matches
        for(int se = SE_START; se <= SE_END; ++se)
        {
            for(final FusionRead read : mReadGroup.Reads)
            {
                if(!isSingleGeneCollection())
                {
                    if(!read.Chromosome.equals(mChromosomes[se]))
                        continue;

                    if(read.GeneCollections[SE_START] != mGeneCollections[se] && read.GeneCollections[SE_END] != mGeneCollections[se])
                        continue;
                }

                if(mJunctionPositions[se] > 0 && !positionWithin(mJunctionPositions[se], read.Positions[SE_START], read.Positions[SE_END]))
                    continue;

                if(se == SE_END && read.getTransExonRefs(SE_END) == null)
                    continue;

                RegionMatchType readMatchType = read.getRegionMatchType(se);

                if(matchRank(readMatchType) < matchRank(mRegionMatchTypes[se]))
                    continue;

                if(matchRank(readMatchType) > matchRank(mRegionMatchTypes[se]))
                {
                    mRegionMatchTypes[se] = readMatchType;
                    mTransExonRefs[se].clear();
                }

                List<FusionTransExon> readTransExonRefs = read.getTransExonRefs(se);
                mergeUnique(mTransExonRefs[se], readTransExonRefs);
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
            FusionTransExon transExonRef = mTransExonRefs[seIndex].get(index);
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
            if(mJunctionTypes[se] == KNOWN)
            {
                continue;
            }
            else if(mJunctionPositions[se] > 0)
            {
                String daBases = junctionSpliceBases[se];

                if(junctionStrands != null)
                {
                    boolean isDonor = (mJunctionOrientations[se] == junctionStrands[se]);

                    if(isDonor && canonicalDonor(daBases, junctionStrands[se]))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                    else if(!isDonor && canonicalAcceptor(daBases, junctionStrands[se]))
                        mJunctionTypes[se] = FusionJunctionType.CANONICAL;
                }
                else
                {
                    // try them both
                    byte asDonorStrand = mJunctionOrientations[se];
                    byte asAcceptorStrand = (byte)(-mJunctionOrientations[se]);

                    if(canonicalDonor(daBases, asDonorStrand) || canonicalAcceptor(daBases, asAcceptorStrand))
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
