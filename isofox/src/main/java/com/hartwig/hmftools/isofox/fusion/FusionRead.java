package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionRead
{
    // skip read ID since always obtainable from the read-group

    public final String Chromosome;
    public final int[] Positions;
    public final byte Orientation;
    public final String Cigar;
    public final String MateChromosome;
    public int MatePosStart;

    public final int[] SoftClipLengths;
    public final String[] BoundaryBases;
    public final int ReadBaseLength;

    public final int[] GeneCollections;
    public final boolean[] IsGenicRegion;
    public boolean HasInterGeneSplit;
    public boolean HasSuppAlignment;
    public final SupplementaryReadData SuppData;
    public boolean IsDuplicate;
    public boolean ContainsSplit;
    public int Flags;

    public final List<int[]> MappedCoords;

    // directly related to fusion junctions, may be set on one, both or no sides
    private final int[] mJunctionPositions;

    private Map<Integer,Integer> mBoundaryDepth; // depth at mapped coords boundaries
    private int[] mJunctionDepth; // depth at chimeric junctions

    private final Map<RegionMatchType,List<TransExonRef>> mTransExonRefs;
    private final Map<RegionMatchType,List<TransExonRef>> mUpperTransExonRefs; // TE refs for upper coords if a spanning read

    public FusionRead(final ReadRecord read)
    {
        Chromosome = read.Chromosome;
        Positions = new int[] { read.PosStart, read.PosEnd};
        Orientation = read.orientation();
        MateChromosome = read.mateChromosome();
        MatePosStart = read.mateStartPosition();
        MappedCoords = read.getMappedRegionCoords(false);
        Cigar = read.Cigar.toString();
        GeneCollections = read.getGeneCollectons();
        IsGenicRegion = read.getIsGenicRegion();
        HasInterGeneSplit = read.hasInterGeneSplit();
        HasSuppAlignment = read.hasSuppAlignment();
        IsDuplicate = read.isDuplicate();
        ContainsSplit = read.containsSplit();
        Flags = read.flags();

        SuppData = read.hasSuppAlignment() ? SupplementaryReadData.from(read.getSuppAlignment()) : null;

        SoftClipLengths = new int[]
                { read.isSoftClipped(SE_START) ? read.Cigar.getFirstCigarElement().getLength() : 0,
                  read.isSoftClipped(SE_END) ? read.Cigar.getLastCigarElement().getLength() : 0 };

        ReadBaseLength = read.Length;
        int extraBasesBuffer = 5;
        int startBases = min(SoftClipLengths[SE_START] + extraBasesBuffer, read.Length);
        int endBases = min(SoftClipLengths[SE_END] + extraBasesBuffer, read.Length);

        BoundaryBases = new String[] { read.ReadBases.substring(0, startBases), read.ReadBases.substring(read.Length - endBases) };

        mTransExonRefs = Maps.newHashMap();
        mUpperTransExonRefs = Maps.newHashMap();

        if(!read.getMappedRegions().isEmpty())
        {
            for(Map.Entry<RegionReadData, RegionMatchType> entry : read.getMappedRegions().entrySet())
            {
                List<TransExonRef> transRefList = mTransExonRefs.get(entry.getValue());

                if(transRefList == null)
                {
                    mTransExonRefs.put(entry.getValue(), Lists.newArrayList(entry.getKey().getTransExonRefs()));
                }
                else
                {
                    transRefList.addAll(entry.getKey().getTransExonRefs());
                }
            }
        }
        else
        {
            mTransExonRefs.putAll(read.getTransExonRefs());
        }

        // depth will only be set for known junctions
        mBoundaryDepth = null;
        mJunctionDepth = null;

        mJunctionPositions = new int[SE_PAIR];

        if(read.junctionPositions() != null)
        {
            mJunctionPositions[SE_START] = read.junctionPositions()[SE_START];
            mJunctionPositions[SE_END] = read.junctionPositions()[SE_END];
            mJunctionDepth = new int[SE_PAIR];
        }
    }

    public int getCoordsBoundary(int se)
    {
        return se == SE_START ? MappedCoords.get(0)[SE_START] : MappedCoords.get(MappedCoords.size() - 1)[SE_END];
    }

    public int posStart() { return Positions[SE_START]; }
    public int posEnd() { return Positions[SE_END]; }
    public boolean spansGeneCollections()
    {
        return GeneCollections[SE_START] != GeneCollections[SE_END];
    }
    public boolean isSoftClipped(int se) { return SoftClipLengths[se] > 0; }
    public boolean isLongestSoftClip(int se) { return SoftClipLengths[se] > SoftClipLengths[switchIndex(se)]; }

    public final int[] junctionPositions() { return mJunctionPositions; }

    public final Map<RegionMatchType,List<TransExonRef>> getTransExonRefs() { return mTransExonRefs; }
    public final Map<RegionMatchType,List<TransExonRef>> getTransExonRefs(int se)
    {
        if(spansGeneCollections())
            return se == SE_START ? mTransExonRefs : mUpperTransExonRefs;
        else
            return mTransExonRefs;
    }

    public void setReadJunctionDepth(final BaseDepth baseDepth)
    {
        if(mJunctionPositions[SE_START] > 0 || mJunctionPositions[SE_END] > 0)
        {
            if(mJunctionDepth == null)
                mJunctionDepth = new int[SE_PAIR];

            if(mJunctionPositions[SE_START] > 0)
                mJunctionDepth[SE_START] = baseDepth.depthAtBase(mJunctionPositions[SE_START]);

            if(mJunctionPositions[SE_END] > 0)
                mJunctionDepth[SE_END] = baseDepth.depthAtBase(mJunctionPositions[SE_END]);

            return;
        }

        if(mBoundaryDepth == null)
            mBoundaryDepth = Maps.newHashMap();

        for(final int[] mappedCoords : MappedCoords)
        {
            mBoundaryDepth.put(mappedCoords[SE_START], baseDepth.depthAtBase(mappedCoords[SE_START]));
            mBoundaryDepth.put(mappedCoords[SE_END], baseDepth.depthAtBase(mappedCoords[SE_END]));
        }
    }

    public int junctionDepth(int se, int junctionPosition)
    {
        if(mJunctionPositions != null && mJunctionDepth != null)
            return mJunctionPositions[se] == junctionPosition ? mJunctionDepth[se] : 0;

        if(mBoundaryDepth != null)
            return mBoundaryDepth.containsKey(junctionPosition) ? mBoundaryDepth.get(junctionPosition) : 0;

        return 0;
    }

    public static List<FusionRead> convertReads(final List<ReadRecord> reads)
    {
        List<FusionRead> fusionReads = reads.stream().map(x -> new FusionRead(x)).collect(Collectors.toList());
        return fusionReads;
    }

    public String toString()
    {
        return String.format("range(%s: %d -> %d) cigar(%s) junc(%d - %d) gc(%d - %d) sup=%s igs=%s",
                Chromosome, Positions[SE_START], Positions[SE_END], Cigar,
                mJunctionPositions != null ? mJunctionPositions[SE_START] : 0, mJunctionPositions != null ? mJunctionPositions[SE_END] : 0,
                GeneCollections[SE_START], GeneCollections[SE_END], HasSuppAlignment, HasInterGeneSplit);
    }

}
