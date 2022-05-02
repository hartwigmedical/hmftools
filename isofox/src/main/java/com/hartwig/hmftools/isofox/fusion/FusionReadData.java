package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.utils.Strings.reverseString;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.CommonUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.JUNCTION_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.setMaxSplitMappedLength;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.GeneData;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionReadData
{
    private final int mId;
    private final String mLocationId;
    private final FusionFragment mFragment; // the fragment which establishes this fusion and whose junction data is used

    private Map<FusionFragmentType,List<FusionFragment>> mFragments; // only cached if all chimeric fragments are written to file
    private final Map<FusionFragmentType,Integer> mFragmentCounts;

    private boolean mIncompleteData;

    private final Set<Integer> mRelatedProximateFusions;
    private final Set<Integer> mRelatedSplicedFusions;

    // not stored by stream
    private final String[] mChromosomes;
    private final int[] mJunctionPositions;
    private final byte[] mJunctionOrientations;
    private final String[] mJunctionBases; // the 10 bases leading up to the junction
    private final String[] mAdjacentJunctionBases; // the 10 bases continuing on from the junction
    private final String[] mJunctionSpliceBases; // the 2 donor/acceptor bases

    private final List<TransExonRef>[] mTransExonRefs;
    private final int[] mReadDepth;

    // the following data is stored by stream, not start/end
    private final GeneData[] mFusionGenes; // up and downstream genes if identified
    private final int[] mStreamIndices; // mapping of up & down stream to position data which is in SV terms
    private final int[] mMaxSplitLengths;
    private final int[] mJunctionHomology;

    public FusionReadData(int id, final FusionFragment fragment)
    {
        mId = id;
        mFragment = fragment;

        mChromosomes = new String[] { fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END] };
        mJunctionPositions = new int[] { fragment.junctionPositions()[SE_START], fragment.junctionPositions()[SE_END] };
        mJunctionOrientations = new byte[]{ fragment.junctionOrientations()[SE_START], fragment.junctionOrientations()[SE_END] };

        mJunctionBases = new String[] {"", ""};
        mJunctionSpliceBases = new String[] {"", ""};
        mAdjacentJunctionBases = new String[] {"", ""};

        mFragments = null;
        mFragmentCounts = Maps.newHashMap();

        mLocationId = fragment.locationPair();

        mRelatedSplicedFusions = Sets.newHashSet();
        mRelatedProximateFusions = Sets.newHashSet();

        mFusionGenes = new GeneData[] { null, null};
        mStreamIndices = new int[] { SE_START, SE_END };

        mReadDepth = new int[] {0, 0};

        // extract depth from reads
        for(int se = SE_START; se <= SE_END; ++se)
        {
            for(FusionRead read : mFragment.reads())
            {
                int depth = read.junctionDepth(se, mJunctionPositions[se]);
                if(depth > 0)
                {
                    mReadDepth[se] = depth;
                    break;
                }
            }
        }

        mMaxSplitLengths = new int[] {0, 0};
        mJunctionHomology = new int[] {0, 0};

        mIncompleteData = false;

        mTransExonRefs = new List[FS_PAIR];
        mTransExonRefs[SE_START] = Lists.newArrayList();
        mTransExonRefs[SE_END] = Lists.newArrayList();

        addFusionFragment(fragment, true);
    }

    public int id() { return mId; }
    public String locationId() { return mLocationId; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] junctionPositions() { return mJunctionPositions; }
    public final byte[] junctionOrientations() { return mJunctionOrientations; }
    public final String[] junctionBases() { return mJunctionBases; }
    public final String[] adjacentJunctionBases() { return mAdjacentJunctionBases; }
    public final String[] junctionSpliceBases() { return mJunctionSpliceBases; }
    public final int[] junctionHomology() { return mJunctionHomology; }
    public final int[] streamIndices() { return mStreamIndices; }

    public boolean hasIncompleteData() { return mIncompleteData; }
    public void setIncompleteData() { mIncompleteData = true; }

    public List<TransExonRef> getTransExonRefsByPos(int se) { return mTransExonRefs[se]; }

    public List<TransExonRef> getTransExonRefsByStream(int fs)
    {
        if(hasViableGenes())
            return mTransExonRefs[mStreamIndices[fs]];

        return mTransExonRefs[fs];
    }

    public final Map<FusionFragmentType,List<FusionFragment>> getFragments() { return mFragments; }

    public final List<FusionFragment> getFragments(FusionFragmentType type)
    {
        return mFragments.containsKey(type) ? mFragments.get(type) : Lists.newArrayList();
    }

    public int getFragmentTypeCount(FusionFragmentType type)
    {
        return mFragmentCounts.containsKey(type) ? mFragmentCounts.get(type) : 0;
    }

    public void addFragmentTypeCount(FusionFragmentType type, int count)
    {
        Integer existingCount = mFragmentCounts.get(type);
        if(existingCount == null)
            mFragmentCounts.put(type, count);
        else
            mFragmentCounts.put(type, existingCount + count);
    }

    public void addFusionFragment(final FusionFragment fragment, boolean cacheFragment)
    {
        addFragmentTypeCount(fragment.type(), 1);

        if(fragment.type().isJunctionType() || fragment == mFragment)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                setMaxSplitMappedLength(se, fragment.readsByLocation(se), mJunctionPositions, mJunctionOrientations, mMaxSplitLengths);
            }
        }

        if(cacheFragment)
        {
            if(mFragments == null)
                mFragments = Maps.newHashMap();

            List<FusionFragment> fragments = mFragments.get(fragment.type());

            if(fragments == null)
                mFragments.put(fragment.type(), Lists.newArrayList(fragment));
            else
                fragments.add(fragment);
        }
    }

    public boolean isKnownSpliced() { return getInitialFragment().isSpliced(); }
    public boolean isUnspliced() { return getInitialFragment().isUnspliced() && getInitialFragment().type().isJunctionType(); }

    public boolean hasViableGenes() { return mFusionGenes[FS_UP] != null && mFusionGenes[FS_DOWN] != null; }

    public void setJunctionBases(final RefGenomeInterface refGenome)
    {
        if(!mFragment.type().isJunctionType())
            return;

        try
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                int junctionBase = mJunctionPositions[se];

                if (junctionOrientations()[se] == POS_ORIENT)
                {
                    String junctionBases = refGenome.getBaseString(
                            mChromosomes[se], junctionBase - JUNCTION_BASE_LENGTH + 1, junctionBase + JUNCTION_BASE_LENGTH);

                    mJunctionBases[se] = junctionBases.substring(0, JUNCTION_BASE_LENGTH);
                    mAdjacentJunctionBases[se] = junctionBases.substring(JUNCTION_BASE_LENGTH);
                    mJunctionSpliceBases[se] = junctionBases.substring(JUNCTION_BASE_LENGTH, JUNCTION_BASE_LENGTH + 2);
                }
                else
                {
                    String junctionBases = refGenome.getBaseString(
                            mChromosomes[se], junctionBase - JUNCTION_BASE_LENGTH, junctionBase + JUNCTION_BASE_LENGTH - 1);

                    mJunctionBases[se] = junctionBases.substring(JUNCTION_BASE_LENGTH);
                    mJunctionSpliceBases[se] = junctionBases.substring(JUNCTION_BASE_LENGTH - 2, JUNCTION_BASE_LENGTH);
                    mAdjacentJunctionBases[se] = junctionBases.substring(0, JUNCTION_BASE_LENGTH);
                }
            }
        }
        catch(Exception e)
        {
            // junction may be in an invalid region, just ignore these
        }

        setHomologyOffsets();
    }

    private void setHomologyOffsets()
    {
        // test moving the junction point back and forth from the current positions

        // in the int-pair array the first element is the number of bases the junction could move
        boolean startHasPosOrient = mJunctionOrientations[SE_START] == POS_ORIENT;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int psIndex = switchIndex(se);

            // first of all the junction start bases will be compared to the adjacent end bases to see if they were transferred from
            // the junction bases to the adjacent bases, would there be a match

            // for the start, compare the start bases with the other junction's post-junction bases
            for(int i = 1; i < JUNCTION_BASE_LENGTH; ++i)
            {
                if(mJunctionBases[se].length() < i || mAdjacentJunctionBases[psIndex].length() < i)
                    break;

                String junctionStr = mJunctionOrientations[se] == POS_ORIENT ?
                        mJunctionBases[se].substring(mJunctionBases[se].length() - i) :
                        mJunctionBases[se].substring(0, i);

                String adjacentStr = mJunctionOrientations[psIndex] == POS_ORIENT ?
                        mAdjacentJunctionBases[psIndex].substring(0, i) :
                        mAdjacentJunctionBases[psIndex].substring(mAdjacentJunctionBases[psIndex].length() - i);

                if(i > 1 && mJunctionOrientations[SE_START] == mJunctionOrientations[SE_END])
                {
                    reverseString(adjacentStr);
                }

                if(!junctionStr.equals(adjacentStr))
                    break;

                boolean testingPosRetreating = (se == SE_START && startHasPosOrient) || (se == SE_END && !startHasPosOrient);

                if(testingPosRetreating)
                    mJunctionHomology[SE_START] = -i;
                else
                    mJunctionHomology[SE_END] = i;
            }
        }
    }

    public boolean matchWithinHomology(final FusionReadData other)
    {
        // must be offset by the same number of bases
        if(mJunctionPositions[SE_START] - other.junctionPositions()[SE_START] != mJunctionPositions[SE_END] - other.junctionPositions()[SE_END])
            return false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(mJunctionOrientations[se] != other.junctionOrientations()[se])
                return false;

            if(mJunctionPositions[se] == other.junctionPositions()[se])
                continue;

            if(mJunctionPositions[se] < other.junctionPositions()[se])
            {
                if(mJunctionPositions[se] + mJunctionHomology[SE_END] < other.junctionPositions()[se])
                    return false;
            }
            else
            {
                if(mJunctionPositions[se] + mJunctionHomology[SE_START] > other.junctionPositions()[se])
                    return false;
            }
        }

        return true;
    }

    private static final int CLOSE_MATCH_BASES = 10;
    private static final double CLOSE_MATCH_FRAG_RATIO = 0.2;

    public boolean isCloseMatch(final FusionReadData other)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(abs(mJunctionPositions[se] - other.junctionPositions()[se]) > CLOSE_MATCH_BASES)
               return false;
        }

        double fragmentsRatio = getFragmentTypeCount(MATCHED_JUNCTION) / (double)other.getFragmentTypeCount(MATCHED_JUNCTION);
        return fragmentsRatio >= 1/CLOSE_MATCH_FRAG_RATIO || fragmentsRatio <= CLOSE_MATCH_FRAG_RATIO;
    }

    public void setStreamData(final List<GeneData> upstreamGenes, final List<GeneData> downstreamGenes, boolean startIsUpstream)
    {
        if(!upstreamGenes.isEmpty())
        {
            mStreamIndices[FS_UP] = startIsUpstream ? SE_START : SE_END;

            if(downstreamGenes.isEmpty())
                mStreamIndices[FS_DOWN] = switchStream(mStreamIndices[FS_UP]);

            mFusionGenes[FS_UP] = upstreamGenes.get(0);
        }

        if(!downstreamGenes.isEmpty())
        {
            mStreamIndices[FS_DOWN] = startIsUpstream ? SE_END : SE_START;

            if(upstreamGenes.isEmpty())
                mStreamIndices[FS_UP] = switchStream(mStreamIndices[FS_DOWN]);

            mFusionGenes[FS_DOWN] = downstreamGenes.get(0);
        }
    }

    public byte[] getGeneStrands()
    {
        if(!hasViableGenes())
            return null;

        if(mStreamIndices[FS_UP] == SE_START)
            return new byte[] { mFusionGenes[SE_START].Strand, mFusionGenes[SE_END].Strand };
        else
            return new byte[] { mFusionGenes[SE_END].Strand, mFusionGenes[SE_START].Strand };
    }

    public final Set<Integer> getRelatedFusions() { return mRelatedSplicedFusions; }

    public void addRelatedFusion(int id, boolean isSpliced)
    {
        if(isSpliced)
            mRelatedSplicedFusions.add(id);
        else
            mRelatedProximateFusions.add(id);
    }

    public StructuralVariantType getImpliedSvType()
    {
        return impliedSvType(mChromosomes, mJunctionOrientations);
    }

    public FusionFragment getInitialFragment() { return mFragment; }

    public void cacheTranscriptData()
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            mTransExonRefs[se].addAll(mFragment.getTransExonRefs()[se]);
        }
    }

    public boolean canAddDiscordantFragment(final FusionFragment fragment, int maxFragmentDistance)
    {
        // a discordant read spans both genes and cannot be outside the standard long fragment length either in intronic terms
        // or by exons if exonic

        // a realigned fragment must touch one of the fusion junctions with soft-clipping

        // the 2 reads' bounds need to fall within 2 or less exons away
        // apply max fragment distance criteria

        int impliedFragmentLength = fragment.reads().get(0).ReadBaseLength * 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<TransExonRef> fragmentRefs = fragment.getTransExonRefs()[se];
            final List<TransExonRef> fusionRefs = getTransExonRefsByPos(se);

            // must match the orientations of the fusion junction
            if(fragment.orientations()[se] != mJunctionOrientations[se])
                return false;

            boolean isUpstream = (mStreamIndices[FS_UP] == se);

            int permittedExonDiff;

            if(isUpstream)
                permittedExonDiff = -2;
            else if(fragment.regionMatchTypes()[se] == INTRON)
                permittedExonDiff = 1;
            else
                permittedExonDiff = 2;

            if(!hasTranscriptExonMatch(fusionRefs, fragmentRefs, permittedExonDiff))
                return false;

            final List<FusionRead> reads = fragment.readsByLocation(se);

            if(reads.isEmpty())
                return false;

            final FusionRead read = reads.get(0);

            int fragmentPosition = read.getCoordsBoundary(switchIndex(se));

            // cannot be on the wrong side of the junction
            if((mJunctionOrientations[se] == 1) == (fragmentPosition > mJunctionPositions[se]))
            {
                // check if a mis-mapping explains the over-hang
                if(!softClippedReadSupportsJunction(read, se))
                    return false;
            }

            // of the fusion is unspliced or the fragment is intronic, then measure the genomic distance vs permitted fragment length
            if(fragment.regionMatchTypes()[se] == INTRON || mFragment.regionMatchTypes()[se] == INTRON)
            {
                if(fragmentPosition < mJunctionPositions[se])
                    impliedFragmentLength += mJunctionPositions[se] - fragmentPosition;
                else
                    impliedFragmentLength += fragmentPosition - mJunctionPositions[se];
            }
        }

        if(impliedFragmentLength > maxFragmentDistance)
            return false;

        return true;
    }

    public boolean canRelignFragmentToJunction(final FusionFragment fragment)
    {
        boolean hasSupportingRead = false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final int seIndex = se;

            for(FusionRead read : fragment.reads())
            {
                if(!read.Chromosome.equals(mChromosomes[seIndex]) || mFragment.geneCollections()[seIndex] != read.GeneCollections[seIndex])
                    continue;

                // check that none of the other reads are on the incorrect side of this fusion junction
                if(mJunctionOrientations[se] == 1 && read.getCoordsBoundary(SE_END) > mJunctionPositions[se] + SOFT_CLIP_JUNC_BUFFER)
                    return false;
                else if(mJunctionOrientations[se] == -1 && read.getCoordsBoundary(SE_START) < mJunctionPositions[se] - SOFT_CLIP_JUNC_BUFFER)
                    return false;

                if(!hasSupportingRead && softClippedReadSupportsJunction(read, se))
                {
                    hasSupportingRead = true;
                }
            }
        }

        return hasSupportingRead;
    }

    private boolean softClippedReadSupportsJunction(final FusionRead read, int juncSeIndex)
    {
        return softClippedReadSupportsJunction(
                read, juncSeIndex, mJunctionPositions[juncSeIndex], mJunctionOrientations[juncSeIndex], mJunctionBases);
    }

    public static boolean softClippedReadSupportsJunction(
            final FusionRead read, int juncSeIndex, int junctionPosition, byte junctionOrientation, final String[] junctionBases)
    {
        // compare a minimum number of soft-clipped bases to the other side of the exon junction
        // if the read extends past break junction, include these bases in what is compared against the next junction to account for homology
        if(junctionOrientation == POS_ORIENT)
        {
            if(!read.isSoftClipped(SE_END))
                return false;

            int readBoundary = read.getCoordsBoundary(SE_END);

            // the fragment is limited to how far past the junction (into the other fused gene) it can overhang by a mis-map
            if(!positionWithin(readBoundary, junctionPosition, junctionPosition + SOFT_CLIP_JUNC_BUFFER))
                return false;

            // test that soft-clipped bases match the other junction's bases
            int scLength = read.SoftClipLengths[SE_END];

            if(scLength < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH || scLength > REALIGN_MAX_SOFT_CLIP_BASE_LENGTH)
                return false;

            if(junctionBases == null)
                return true;

            // if the junction is 1 base higher, then take 1 base off the soft-clipped bases
            int posAdjust = readBoundary > junctionPosition ? readBoundary - junctionPosition : 0;

            int boundaryBaseLength = read.BoundaryBases[SE_END].length();
            String extraBases = read.BoundaryBases[SE_END].substring(boundaryBaseLength - scLength - posAdjust, boundaryBaseLength);
            // String extraBases = read.ReadBases.substring(read.Length - scLength - posAdjust, read.Length);

            if(extraBases.length() > JUNCTION_BASE_LENGTH)
                extraBases = extraBases.substring(0, JUNCTION_BASE_LENGTH);

            return junctionBases[switchIndex(juncSeIndex)].startsWith(extraBases);
        }
        else
        {
            if(!read.isSoftClipped(SE_START))
                return false;

            int readBoundary = read.getCoordsBoundary(SE_START);

            if(!positionWithin(readBoundary, junctionPosition - SOFT_CLIP_JUNC_BUFFER, junctionPosition))
                return false;

            int scLength = read.SoftClipLengths[SE_START];

            if(scLength < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH || scLength > REALIGN_MAX_SOFT_CLIP_BASE_LENGTH)
                return false;

            if(junctionBases == null)
                return true;

            int posAdjust = readBoundary < junctionPosition ? junctionPosition - readBoundary : 0;

            String extraBases = read.BoundaryBases[SE_START].substring(0, scLength + posAdjust);
            // String extraBases = read.ReadBases.substring(0, scLength + posAdjust);

            if(extraBases.length() > JUNCTION_BASE_LENGTH)
                extraBases = extraBases.substring(extraBases.length() - JUNCTION_BASE_LENGTH, extraBases.length());

            return junctionBases[switchIndex(juncSeIndex)].endsWith(extraBases);
        }
    }

    public int[] getReadDepth() { return mReadDepth; }
    public int[] getMaxSplitLengths() { return mMaxSplitLengths; }

    public String getGeneName(int stream)
    {
        return mFusionGenes[stream] != null ? mFusionGenes[stream].GeneName : "";
    }

    public String toString()
    {
        return String.format("%d: chr(%s-%s) junc(%d-%d %d/%d %s) genes(%s-%s) frags(%d)",
                mId, mChromosomes[SE_START], mChromosomes[SE_END], mJunctionPositions[SE_START], mJunctionPositions[SE_END],
                mJunctionOrientations[SE_START], mJunctionOrientations[SE_END], getImpliedSvType(),
                getGeneName(FS_UP), getGeneName(FS_DOWN), mFragmentCounts.values().stream().mapToInt(x -> x).sum());
    }

    public static final String FUSION_ID_PREFIX = "Id_";
    public static final String FUSION_NONE = "NONE";

    public static String fusionId(int id) { return String.format("%s%d", FUSION_ID_PREFIX, id); }

    public FusionData toFusionData()
    {
        boolean isValid = hasViableGenes() && !hasIncompleteData();

        String[] chromosomes = new String[FS_PAIR];
        int[] junctionPositions = new int[FS_PAIR];
        byte[] junctionOrientations = new byte[FS_PAIR];

        final FusionFragment sampleFragment = getInitialFragment();

        FusionJunctionType[] junctionTypes = new FusionJunctionType[FS_PAIR];

        byte[] strands = new byte[FS_PAIR];
        String[] geneIds = new String[] {"", ""};
        String[] geneNames = new String[] {"", "" };

        int[] coverage = new int[FS_PAIR];
        int[] anchorDistance = new int[FS_PAIR];
        String[] transData = new String[] {"", ""};

        int splitFragments = 0;
        int realignedFragments = 0;
        int discordantFragments = 0;

        for(Map.Entry<FusionFragmentType,Integer> entry : mFragmentCounts.entrySet())
        {
            if(entry.getKey() == MATCHED_JUNCTION)
                splitFragments = entry.getValue();
            else if(entry.getKey() == DISCORDANT || entry.getKey() == DISCORDANT_JUNCTION)
                discordantFragments = entry.getValue();
            else if(entry.getKey() == REALIGNED)
                realignedFragments = entry.getValue();
        }

        int totalFragments = splitFragments + realignedFragments + discordantFragments;

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            chromosomes[fs] = mChromosomes[mStreamIndices[fs]];

            junctionPositions[fs] = mJunctionPositions[mStreamIndices[fs]];
            junctionOrientations[fs] = mJunctionOrientations[mStreamIndices[fs]];
            junctionTypes[fs] = sampleFragment.junctionTypes()[mStreamIndices[fs]];

            final GeneData geneData = mFusionGenes[fs];

            if(geneData != null)
            {
                geneIds[fs] = geneData.GeneId;
                geneNames[fs] = geneData.GeneName;
                strands[fs] = geneData.Strand;
            }

            // since depth of 1 may have been discarded from the BaseDepth, correct for this
            coverage[fs] = max(mReadDepth[mStreamIndices[fs]], splitFragments);

            anchorDistance[fs] = mMaxSplitLengths[mStreamIndices[fs]];

            final List<TransExonRef> transExonRefs = getTransExonRefsByStream(fs);
            if(transExonRefs.isEmpty())
            {
                transData[fs] = FUSION_NONE;
            }
            else
            {
                StringJoiner td = new StringJoiner(ITEM_DELIM);
                for(final TransExonRef transExonRef : transExonRefs)
                {
                    td.add(String.format("%s-%d", transExonRef.TransName, transExonRef.ExonRank));
                }

                transData[fs] = td.toString();
            }
        }

        String relatedFusionIds = "";

        if(!mRelatedSplicedFusions.isEmpty())
        {
            StringJoiner rf = new StringJoiner(ITEM_DELIM);
            mRelatedSplicedFusions.stream().map(x -> fusionId(x)).forEach(x -> rf.add(x));
            relatedFusionIds = rf.toString();
        }
        else
        {
            relatedFusionIds = FUSION_NONE;
        }

        String readType = sampleFragment.hasSuppAlignment() ? "SuppAlign" :
                (sampleFragment.reads().stream().anyMatch(x -> x.ContainsSplit) ? "Split" : (
                        sampleFragment.type() == DISCORDANT_JUNCTION ? "Discordant" : "Other"));

        return new FusionData(
                mId, isValid, chromosomes, junctionPositions, junctionOrientations, junctionTypes, getImpliedSvType().toString(),
                readType, geneIds, geneNames, strands, totalFragments, splitFragments, realignedFragments, discordantFragments, coverage,
                anchorDistance, transData, relatedFusionIds, getInitialFragment().readId());
    }
}
