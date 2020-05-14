package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsWithin;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.JUNCTION_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MAX_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.REALIGN_MIN_SOFT_CLIP_BASE_LENGTH;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.SOFT_CLIP_JUNC_BUFFER;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.FS_DOWNSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.FS_PAIR;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.FS_UPSTREAM;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.formLocationPair;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.switchStream;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;

import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionReadData
{
    private final int mId;
    private final String mLocationId;
    private final FusionFragment mFragment; // the one which establishes this fusion and whose junction data is used

    private final Map<FusionFragmentType,List<FusionFragment>> mFragments;

    private boolean mIncompleteData;

    private final List<Integer> mRelatedFusions;

    // not stored by stream
    private final String[] mChromosomes;
    private final int[] mGeneCollections;
    private final int[] mJunctionPositions;
    private final byte[] mJunctionOrientations;
    private final List<TransExonRef>[] mTransExonRefs;
    private final int[] mReadDepth;

    // the following data is stored by stream, not start/end
    private final List<EnsemblGeneData>[] mCandidateGenes; // up and downstream genes
    private final String[] mFusionGeneIds;
    private final int[] mStreamIndices; // mapping of up & down stream to position data which is in SV terms

    public FusionReadData(int id, final FusionFragment fragment)
    {
        mId = id;
        mFragment = fragment;

        mChromosomes = new String[] { fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END] };
        mGeneCollections = new int[] { fragment.geneCollections()[SE_START], fragment.geneCollections()[SE_END] };
        mJunctionPositions = new int[] { fragment.junctionPositions()[SE_START], fragment.junctionPositions()[SE_END] };
        mJunctionOrientations = new byte[]{ fragment.junctionOrientations()[SE_START], fragment.junctionOrientations()[SE_END] };

        mFragments = Maps.newHashMap();
        addFusionFragment(fragment);

        mLocationId = formLocationPair(mChromosomes, mGeneCollections, fragment.inGenicRegions());

        mRelatedFusions = Lists.newArrayList();
        mFusionGeneIds = new String[] {"", ""};
        mStreamIndices = new int[] { SE_START, SE_END };
        mReadDepth = new int[] {0, 0};

        mCandidateGenes = new List[FS_PAIR];
        mCandidateGenes[SE_START] = Lists.newArrayList();
        mCandidateGenes[SE_END] = Lists.newArrayList();
        mIncompleteData = false;

        mTransExonRefs = new List[FS_PAIR];
        mTransExonRefs[SE_START] = Lists.newArrayList();
        mTransExonRefs[SE_END] = Lists.newArrayList();
    }

    public int id() { return mId; }
    public String locationId() { return mLocationId; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final int[] junctionPositions() { return mJunctionPositions; }
    public final byte[] junctionOrientations() { return mJunctionOrientations; }
    private final String[] junctionBases() { return mFragment.junctionBases(); }

    public boolean hasIncompleteData() { return mIncompleteData; }
    public void setIncompleteData() { mIncompleteData = true; }

    public List<TransExonRef> getTransExonRefsByPos(int se) { return mTransExonRefs[se]; }

    public List<TransExonRef> getTransExonRefsByStream(int fs)
    {
        if(hasViableGenes())
            return mTransExonRefs[mStreamIndices[fs]];

        return mTransExonRefs[fs];
    }

    public final List<FusionFragment> getAllFragments()
    {
        if(mFragments.size() == 1)
            return mFragments.values().iterator().next();

        final List<FusionFragment> fragments = Lists.newArrayList();
        mFragments.values().forEach(x -> fragments.addAll(x));
        return fragments;
    }

    public int fragmentCount() { return mFragments.values().stream().mapToInt(x -> x.size()).sum(); }

    public final Map<FusionFragmentType,List<FusionFragment>> getFragments() { return mFragments; }
    public final List<FusionFragment> getFragments(FusionFragmentType type)
    {
        return mFragments.containsKey(type) ? mFragments.get(type) : Lists.newArrayList();
    }

    public void addFusionFragment(final FusionFragment fragment)
    {
        List<FusionFragment> fragments = mFragments.get(fragment.type());

        if (fragments == null)
        {
            mFragments.put(fragment.type(), Lists.newArrayList(fragment));
            return;
        }

        fragments.add(fragment);
    }

    public boolean hasJunctionFragments() { return mFragments.containsKey(MATCHED_JUNCTION); }

    public boolean isKnownSpliced() { return getInitialFragment().isSpliced(); }
    public boolean isUnspliced() { return getInitialFragment().isUnspliced() && getInitialFragment().type() == MATCHED_JUNCTION; }

    public boolean hasViableGenes() { return !mCandidateGenes[FS_UPSTREAM].isEmpty() && !mCandidateGenes[FS_DOWNSTREAM].isEmpty(); }

    public boolean isValid() { return hasViableGenes() && !hasIncompleteData(); }

    public void setStreamData(final List<EnsemblGeneData> upstreamGenes, final List<EnsemblGeneData> downstreamGenes, boolean startIsUpstream)
    {
        if(!upstreamGenes.isEmpty())
        {
            mStreamIndices[FS_UPSTREAM] = startIsUpstream ? SE_START : SE_END;

            if(downstreamGenes.isEmpty())
                mStreamIndices[FS_DOWNSTREAM] = switchStream(mStreamIndices[FS_UPSTREAM]);

            mCandidateGenes[FS_UPSTREAM] = upstreamGenes;

            // until a more informed decision can be made
            mFusionGeneIds[FS_UPSTREAM] = upstreamGenes.get(0).GeneId;
        }

        if(!downstreamGenes.isEmpty())
        {
            mStreamIndices[FS_DOWNSTREAM] = startIsUpstream ? SE_END : SE_START;

            if(upstreamGenes.isEmpty())
                mStreamIndices[FS_UPSTREAM] = switchStream(mStreamIndices[FS_DOWNSTREAM]);

            mCandidateGenes[FS_DOWNSTREAM] = downstreamGenes;
            mFusionGeneIds[FS_DOWNSTREAM] = downstreamGenes.get(0).GeneId;
        }
    }

    public byte[] getGeneStrands()
    {
        if(!hasViableGenes())
            return null;

        if(mStreamIndices[FS_UPSTREAM] == SE_START)
            return new byte[] { mCandidateGenes[SE_START].get(0).Strand, mCandidateGenes[SE_END].get(0).Strand };
        else
            return new byte[] { mCandidateGenes[SE_END].get(0).Strand, mCandidateGenes[SE_START].get(0).Strand };
    }

    public final List<Integer> getRelatedFusions() { return mRelatedFusions; }

    public void addRelatedFusion(int id)
    {
        if(!mRelatedFusions.contains(id))
            mRelatedFusions.add(id);
    }

    public StructuralVariantType getImpliedSvType()
    {
        return impliedSvType(mChromosomes, mJunctionOrientations);
    }

    public boolean junctionMatch(final FusionFragment fragment)
    {
        return fragment.type() == MATCHED_JUNCTION
                && mJunctionPositions[SE_START] == fragment.junctionPositions()[SE_START]
                && mJunctionPositions[SE_END] == fragment.junctionPositions()[SE_END]
                && mJunctionOrientations[SE_START] == fragment.junctionOrientations()[SE_START]
                && mJunctionOrientations[SE_END] == fragment.junctionOrientations()[SE_END];
    }

    public FusionFragment getInitialFragment() { return mFragment; }

    public void cacheTranscriptData()
    {
        for (int se = SE_START; se <= SE_END; ++se)
        {
            mTransExonRefs[se].addAll(mFragment.getTransExonRefs()[se]);
        }
    }

    public boolean canAddUnfusedFragment(final FusionFragment fragment, int maxFragmentDistance)
    {
        // a discordant read spans both genes and cannot be outside the standard long fragment length either in intronic terms
        // or by exons if exonic

        // a realigned fragment must touch one of the fusion junctions with soft-clipping

        // the 2 reads' bounds need to fall within 2 or less exons away
        // apply max fragment distance criteria

        int impliedFragmentLength = fragment.reads().get(0).Length * 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<TransExonRef> fragmentRefs = fragment.getTransExonRefs()[se];
            final List<TransExonRef> fusionRefs = getTransExonRefsByPos(se);

            // must match the orientations of the fusion junction
            if(fragment.orientations()[se] != mJunctionOrientations[se])
                return false;

            boolean isUpstream = (mStreamIndices[FS_UPSTREAM] == se);

            int permittedExonDiff;

            if(isUpstream)
                permittedExonDiff = -2;
            else if(fragment.regionMatchTypes()[se] == INTRON)
                permittedExonDiff = 1;
            else
                permittedExonDiff = 2;

            if(!hasTranscriptExonMatch(fusionRefs, fragmentRefs, permittedExonDiff))
                return false;

            final List<ReadRecord> reads = fragment.readsByLocation(se);

            if(reads.isEmpty())
                return false;

            final ReadRecord read = reads.get(0);

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
                if (fragmentPosition < mJunctionPositions[se])
                    impliedFragmentLength += mJunctionPositions[se] - fragmentPosition;
                else
                    impliedFragmentLength += fragmentPosition - mJunctionPositions[se];
            }
        }

        if(impliedFragmentLength > maxFragmentDistance)
            return false;

        return true;
    }

    public boolean isRelignedFragment(final FusionFragment fragment)
    {
        boolean hasSupportingRead = false;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            final int seIndex = se;
            List<ReadRecord> reads = fragment.reads().stream()
                    .filter(x -> mChromosomes[seIndex].equals(x.Chromosome))
                    .filter(x -> mGeneCollections[seIndex] == x.getGeneCollectons()[seIndex])
                    .filter(x -> mFragment.inGenicRegions()[seIndex] == x.getIsGenicRegion()[seIndex])
                    .collect(Collectors.toList());

            for (ReadRecord read : reads)
            {
                if (softClippedReadSupportsJunction(read, se))
                {
                    hasSupportingRead = true;
                }

                // check that none of the other reads are on the incorrect side of this fusion junction
                if(mJunctionOrientations[se] == 1 && read.getCoordsBoundary(SE_END) > mJunctionPositions[se] + SOFT_CLIP_JUNC_BUFFER)
                    return false;
                else if(mJunctionOrientations[se] == -1 && read.getCoordsBoundary(SE_START) < mJunctionPositions[se] - SOFT_CLIP_JUNC_BUFFER)
                    return false;
            }
        }

        return hasSupportingRead;
    }

    private boolean softClippedReadSupportsJunction(final ReadRecord read, int juncSeIndex)
    {
        // compare a minimum number of soft-clipped bases to the other side of the exon junction
        // if the read extends past break junction, include these bases in what is compared against the next junction to account for homology
        if(mJunctionOrientations[juncSeIndex] == 1)
        {
            if(!read.isSoftClipped(SE_END))
                return false;

            int readBoundary = read.getCoordsBoundary(SE_END);

            // the fragment is limited to how far past the junction (into the other fused gene) it can overhang by a mis-map
            if(!positionWithin(readBoundary, mJunctionPositions[juncSeIndex], mJunctionPositions[juncSeIndex] + SOFT_CLIP_JUNC_BUFFER))
                return false;

            // test that soft-clipped bases match the other junction's bases
            int scLength = read.Cigar.getLastCigarElement().getLength();

            if(scLength < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH || scLength > REALIGN_MAX_SOFT_CLIP_BASE_LENGTH)
                return false;

            // if the junction is 1 base higher, then take 1 base off the soft-clipped bases
            int posAdjust = readBoundary > mJunctionPositions[juncSeIndex] ? readBoundary - mJunctionPositions[juncSeIndex] : 0;

            String extraBases = read.ReadBases.substring(read.Length - scLength - posAdjust, read.Length);

            if(extraBases.length() > JUNCTION_BASE_LENGTH)
                extraBases = extraBases.substring(0, JUNCTION_BASE_LENGTH);

            return junctionBases()[switchIndex(juncSeIndex)].startsWith(extraBases);
        }
        else
        {
            if(!read.isSoftClipped(SE_START))
                return false;

            int readBoundary = read.getCoordsBoundary(SE_START);

            if(!positionWithin(readBoundary, mJunctionPositions[juncSeIndex] - SOFT_CLIP_JUNC_BUFFER, mJunctionPositions[juncSeIndex]))
                return false;

            int scLength = read.Cigar.getFirstCigarElement().getLength();

            if(scLength < REALIGN_MIN_SOFT_CLIP_BASE_LENGTH || scLength > REALIGN_MAX_SOFT_CLIP_BASE_LENGTH)
                return false;

            int posAdjust = readBoundary < mJunctionPositions[juncSeIndex] ? mJunctionPositions[juncSeIndex] - readBoundary : 0;

            String extraBases = read.ReadBases.substring(0, scLength + posAdjust);

            if(extraBases.length() > JUNCTION_BASE_LENGTH)
                extraBases = extraBases.substring(extraBases.length() - JUNCTION_BASE_LENGTH, extraBases.length());

            return junctionBases()[switchIndex(juncSeIndex)].endsWith(extraBases);
        }
    }

    private int[] maxSplitMappedLength()
    {
        // find the longest section mapped across the junction
        final List<FusionFragment> fragments = mFragments.get(MATCHED_JUNCTION);

        int[] maxSplitLengths = new int[SE_PAIR];

        if(fragments == null)
            return maxSplitLengths;

        for (int se = SE_START; se <= SE_END; ++se)
        {
            final int seIndex = se;
            for(final FusionFragment fragment : fragments)
            {
                final List<ReadRecord> reads = !fragment.isSingleGene() ? fragment.readsByLocation(se) :
                        fragment.reads().stream()
                                .filter(x -> positionWithin(mJunctionPositions[seIndex], x.PosStart, x.PosEnd)).collect(Collectors.toList());

                if(reads.isEmpty())
                {
                    ISF_LOGGER.error("fusion({}) fragment({}) has no reads on side", toString(), fragment.readId());
                    continue;
                }

                int mappedBases;
                if(reads.size() == 1)
                {
                    mappedBases = reads.get(0).getMappedRegionCoords(false).stream().mapToInt(x -> x[SE_END] - x[SE_START]).sum();
                }
                else
                {
                    final List<int[]> sharedRegions = deriveCommonRegions(
                            reads.get(0).getMappedRegionCoords(false), reads.get(1).getMappedRegionCoords(false));
                    mappedBases = sharedRegions.stream().mapToInt(x -> x[SE_END] - x[SE_START]).sum();
                }

                maxSplitLengths[se] = max(mappedBases, maxSplitLengths[se]);
            }
        }

        return maxSplitLengths;
    }

    public int[] getReadDepth() { return mReadDepth; }

    public String getGeneName(int stream)
    {
        if(mCandidateGenes[stream].isEmpty())
            return "";

        if(mFusionGeneIds[stream].isEmpty())
            return mCandidateGenes[stream].get(0).GeneName;

        return mCandidateGenes[stream].stream()
                .filter(x -> x.GeneId.equals(mFusionGeneIds[stream])).findFirst().map(x -> x.GeneName).orElse("");
    }

    public String toString()
    {
        return String.format("%d: chr(%s-%s) junc(%d-%d %d/%d %s) genes(%s-%s) frags(%d)",
                mId, mChromosomes[SE_START], mChromosomes[SE_END], mJunctionPositions[SE_START], mJunctionPositions[SE_END],
                mJunctionOrientations[SE_START], mJunctionOrientations[SE_END], getImpliedSvType(),
                getGeneName(FS_UPSTREAM), getGeneName(FS_DOWNSTREAM), mFragments.size());
    }

    public static String csvHeader()
    {
        return "FusionId,Valid,GeneIdUp,GeneNameUp,ChrUp,PosUp,OrientUp,StrandUp,JuncTypeUp"
                + ",GeneIdDown,GeneNameDown,ChrDown,PosDown,OrientDown,StrandDown,JuncTypeDown"
                + ",SVType,NonSupp,TotalFragments,SplitFrags,RealignedFrags,DiscordantFrags,MultiMapFrags,CoverageUp,CoverageDown"
                + ",MaxAnchorLengthUp,MaxAnchorLengthDown,TransDataUp,TransDataDown,OtherGenesUp,OtherGenesDown,RelatedFusions,RefReadId";
    }

    public static String fusionId(int id) { return String.format("Id_%d", id); }

    public String toCsv()
    {
        StringJoiner csvData = new StringJoiner(DELIMITER);

        csvData.add(fusionId(mId));
        csvData.add(String.valueOf(hasViableGenes() && !hasIncompleteData()));

        final FusionFragment sampleFragment = getInitialFragment();

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            final String geneId = mFusionGeneIds[fs];
            final List<EnsemblGeneData> genes = mCandidateGenes[fs];

            csvData.add(geneId);

            final EnsemblGeneData geneData = genes.stream()
                    .filter(x -> x.GeneId.equals(geneId)).findFirst().map(x -> x).orElse(null);

            csvData.add(geneData != null ? geneData.GeneName : "");

            csvData.add(mChromosomes[mStreamIndices[fs]]);
            csvData.add(String.valueOf(mJunctionPositions[mStreamIndices[fs]]));
            csvData.add(String.valueOf(mJunctionOrientations[mStreamIndices[fs]]));

            csvData.add(geneData != null ? String.valueOf(geneData.Strand) : "0");
            csvData.add(sampleFragment.junctionTypes()[mStreamIndices[fs]].toString());
        }

        csvData.add(getImpliedSvType().toString());
        csvData.add(String.valueOf(!sampleFragment.hasSuppAlignment()));

        int splitFragments = 0;
        int realignedFragments = 0;
        int discordantFragments = 0;
        int readsWithSecondaries = 0;

        for(Map.Entry<FusionFragmentType,List<FusionFragment>> entry : mFragments.entrySet())
        {
            if(entry.getKey() == MATCHED_JUNCTION)
                splitFragments = entry.getValue().size();
            else if(entry.getKey() == DISCORDANT)
                discordantFragments = entry.getValue().size();
            else if(entry.getKey() == REALIGNED)
                realignedFragments = entry.getValue().size();

            readsWithSecondaries += (int)entry.getValue().stream()
                    .filter(x -> x.reads().stream().anyMatch(y -> y.getSecondaryReadCount() > 0)).count();
        }

        int totalFragments = splitFragments + realignedFragments + discordantFragments;

        csvData.add(String.valueOf(totalFragments));
        csvData.add(String.valueOf(splitFragments));
        csvData.add(String.valueOf(realignedFragments));
        csvData.add(String.valueOf(discordantFragments));
        csvData.add(String.valueOf(readsWithSecondaries));

        // since depth of 1 may have been discarded from the BaseDepth, correct for this
        csvData.add(String.valueOf(max(mReadDepth[mStreamIndices[FS_UPSTREAM]], splitFragments)));
        csvData.add(String.valueOf(max(mReadDepth[mStreamIndices[FS_DOWNSTREAM]], splitFragments)));

        final int[] maxSplitLengths = maxSplitMappedLength();
        csvData.add(String.valueOf(maxSplitLengths[mStreamIndices[FS_UPSTREAM]]));
        csvData.add(String.valueOf(maxSplitLengths[mStreamIndices[FS_DOWNSTREAM]]));

        for (int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            final List<TransExonRef> transExonRefs = getTransExonRefsByStream(fs);
            if(transExonRefs.isEmpty())
            {
                csvData.add("NONE");
                continue;
            }

            String transData = "";
            for(final TransExonRef transExonRef : transExonRefs)
            {
                transData = appendStr(transData, String.format("%s-%d", transExonRef.TransName, transExonRef.ExonRank), ';');
            }

            csvData.add(transData);
        }

        if(hasViableGenes())
        {
            String[] otherGenes = new String[] {"", ""};

            for (int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
            {
                for (final EnsemblGeneData geneData : mCandidateGenes[fs])
                {
                    if (!geneData.GeneId.equals(mFusionGeneIds[fs]))
                    {
                        otherGenes[fs] = appendStr(otherGenes[fs], geneData.GeneName, ';');
                    }
                }

                csvData.add(!otherGenes[fs].isEmpty() ? otherGenes[fs] : "NONE");
            }
        }
        else
        {
            csvData.add("NONE");
            csvData.add("NONE");
        }

        if(!mRelatedFusions.isEmpty())
        {
            List<String> relatedFusions = mRelatedFusions.stream().map(x -> fusionId(x)).collect(Collectors.toList());
            csvData.add(appendStrList(relatedFusions, ';'));
        }
        else
        {
            csvData.add("NONE");
        }

        csvData.add(ISF_LOGGER.isDebugEnabled() ? getInitialFragment().readId() : "-");

        return csvData.toString();
    }

}
