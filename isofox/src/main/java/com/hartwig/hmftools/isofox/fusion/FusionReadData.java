package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.BOTH_JUNCTIONS;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.ONE_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionReadData
{
    private final int mId;
    private final String mLocationId;

    private final Map<FusionFragmentType,List<FusionFragment>> mFragments;

    private boolean mIncompleteData;

    private final List<Integer> mRelatedFusions;

    private final String[] mChromosomes;
    private final int[] mGeneCollections;
    private final int[] mJunctionPositions;
    private final byte[] mJunctionOrientations;
    private final List<TransExonRef>[] mTransExonRefs; // not stored by stream
    private final int[] mReadDepth;

    // the following data is stored by stream, not start/end
    private final List<EnsemblGeneData>[] mCandidateGenes; // up and downstream genes
    private final String[] mFusionGeneIds;
    private final int[] mStreamIndices; // mapping of up & down stream to position data which is in SV terms

    public static final int FS_UPSTREAM = 0;
    public static final int FS_DOWNSTREAM = 1;
    public static final int FS_PAIR = 2;

    public FusionReadData(int id, final FusionFragment fragment)
    {
        mId = id;

        mChromosomes = new String[] { fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END] };
        mGeneCollections = new int[] { fragment.geneCollections()[SE_START], fragment.geneCollections()[SE_END] };
        mJunctionPositions = new int[] { fragment.junctionPositions()[SE_START], fragment.junctionPositions()[SE_END] };
        mJunctionOrientations = new byte[]{ fragment.junctionOrientations()[SE_START], fragment.junctionOrientations()[SE_END] };

        mFragments = Maps.newHashMap();
        addFusionFragment(fragment);

        mLocationId = formLocationPair(mChromosomes, mGeneCollections);

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

        // record depth by any overlap in the fusion's junctions
        if(fragment.type() == BOTH_JUNCTIONS)
        {
            ++mReadDepth[SE_START];
            ++mReadDepth[SE_END];
        }
        else
        {
            for (int se = SE_START; se <= SE_END; ++se)
            {
                for (ReadRecord read : fragment.getReads())
                {
                    if (!read.Chromosome.equals(mChromosomes[se]))
                        continue;

                    final int seIndex = se;
                    if (read.getMappedRegionCoords().stream().anyMatch(x -> positionWithin(mJunctionPositions[seIndex], x[SE_START], x[SE_END])))
                    {
                        ++mReadDepth[se];
                        break;
                    }
                }
            }
        }
    }

    public boolean hasJunctionFragments() { return mFragments.containsKey(BOTH_JUNCTIONS); }

    public boolean isKnownSpliced() { return getSampleFragment().isSpliced(); }
    public boolean isUnspliced() { return getSampleFragment().isUnspliced() && getSampleFragment().hasBothJunctions(); }

    public List<EnsemblGeneData>[] getCandidateGenes() { return mCandidateGenes; }

    public boolean hasViableGenes() { return !mCandidateGenes[FS_UPSTREAM].isEmpty() && !mCandidateGenes[FS_DOWNSTREAM].isEmpty(); }

    public boolean isValid() { return hasViableGenes() && !hasIncompleteData(); }

    public void setStreamData(final List<EnsemblGeneData> upstreamGenes, final List<EnsemblGeneData> downstreamGenes, boolean startIsUpstream)
    {
        mStreamIndices[FS_UPSTREAM] = startIsUpstream ? SE_START : SE_END;
        mStreamIndices[FS_DOWNSTREAM] = startIsUpstream ? SE_END : SE_START;
        mCandidateGenes[FS_UPSTREAM] = upstreamGenes;
        mCandidateGenes[FS_DOWNSTREAM] = downstreamGenes;

        // until a more informed decision can be made
        mFusionGeneIds[FS_UPSTREAM] = upstreamGenes.get(0).GeneId;
        mFusionGeneIds[FS_DOWNSTREAM] = downstreamGenes.get(0).GeneId;
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

    public String chrPair() { return formChromosomePair(mChromosomes[SE_START], mChromosomes[SE_END]); }

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
        return fragment.hasBothJunctions()
                && mJunctionPositions[SE_START] == fragment.junctionPositions()[SE_START] && mJunctionPositions[SE_END] == fragment.junctionPositions()[SE_END]
                && mJunctionOrientations[SE_START] == fragment.junctionOrientations()[SE_START] && mJunctionOrientations[SE_END] == fragment.junctionOrientations()[SE_END];
    }

    public FusionFragment getSampleFragment()
    {
        if(mFragments.containsKey(BOTH_JUNCTIONS))
        {
            return mFragments.get(BOTH_JUNCTIONS).get(0);
        }
        else if(mFragments.containsKey(ONE_JUNCTION))
        {
            return mFragments.get(ONE_JUNCTION).get(0);
        }
        else
        {
            return mFragments.values().iterator().next().get(0);
        }
    }

    public void cacheTranscriptData()
    {
        // select a sample fragment from which to extract transcript-exon data
        final FusionFragment sampleFragment = getSampleFragment();

        for (int se = SE_START; se <= SE_END; ++se)
        {
            mTransExonRefs[se].addAll(sampleFragment.getTransExonRefs()[se]);
        }
    }

    public boolean canAddDiscordantFragment(final FusionFragment fragment, int maxFragmentDistance)
    {
        // a discordant read spans both reads and cannot be outside the standard long fragment length either in intronic terms
        // or by exons if exonic

        // the 2 reads' bounds need to fall within 2 or less exons away
        // apply max fragment distance criteria

        int impliedFragmentLength = fragment.getReads().get(0).Length * 2;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final List<TransExonRef> fragmentRefs = fragment.getTransExonRefs()[se];
            final List<TransExonRef> fusionRefs = getTransExonRefsByPos(se);

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

            int fragmentPosition = fragment.getReadBoundary(mChromosomes[se], mGeneCollections[se], switchIndex(se));

            // cannot be on the wrong side of the junction
            if((mJunctionOrientations[se] == 1) == (fragmentPosition > mJunctionPositions[se]))
                return false;

            if(fragment.regionMatchTypes()[se] == INTRON)
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
                + ",SVType,TotalFragments,SplitFrags,RealignedFrags,DiscordantFrags,SingleFrags,JuncDepthUp,JuncDepthDown"
                + ",TransDataUp,TransDataDown,OtherGenesUp,OtherGenesDown,RelatedFusions";
    }

    public static String fusionId(int id) { return String.format("Id_%d", id); }

    public String toCsv()
    {
        StringJoiner csvData = new StringJoiner(DELIMITER);

        csvData.add(fusionId(mId));
        csvData.add(String.valueOf(hasViableGenes() && !hasIncompleteData()));

        final FusionFragment sampleFragment = getSampleFragment();

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            final String geneId = mFusionGeneIds[fs];
            final List<EnsemblGeneData> genes = mCandidateGenes[fs];

            csvData.add(geneId);

            final EnsemblGeneData geneData = genes.stream()
                    .filter(x -> x.GeneId.equals(geneId)).findFirst().map(x -> x).orElse(null);

            if(geneData != null)
            {
                csvData.add(geneData.GeneName);

                csvData.add(mChromosomes[mStreamIndices[fs]]);
                csvData.add(String.valueOf(mJunctionPositions[mStreamIndices[fs]]));
                csvData.add(String.valueOf(mJunctionOrientations[mStreamIndices[fs]]));
                csvData.add(String.valueOf(geneData.Strand));
                csvData.add(sampleFragment.junctionTypes()[mStreamIndices[fs]].toString());
            }
            else
            {
                csvData.add("");
                csvData.add(mChromosomes[fs]);
                csvData.add(String.valueOf(mJunctionPositions[fs]));
                csvData.add(String.valueOf(mJunctionOrientations[fs]));
                csvData.add("0");
                csvData.add(FusionJunctionType.UNKNOWN.toString());
            }
        }

        csvData.add(getImpliedSvType().toString());

        int splitFragments = mFragments.containsKey(BOTH_JUNCTIONS) ? mFragments.get(BOTH_JUNCTIONS).size() : 0;
        int realignedFragments = mFragments.containsKey(REALIGNED) ? mFragments.get(REALIGNED).size() : 0;
        int discordantFragments = mFragments.containsKey(DISCORDANT) ? mFragments.get(DISCORDANT).size() : 0;
        int singleJuncFragments = mFragments.containsKey(ONE_JUNCTION) ? mFragments.get(ONE_JUNCTION).size() : 0;

        int totalFragments = splitFragments + realignedFragments + discordantFragments + singleJuncFragments;

        csvData.add(String.valueOf(totalFragments));
        csvData.add(String.valueOf(splitFragments));
        csvData.add(String.valueOf(realignedFragments));
        csvData.add(String.valueOf(discordantFragments));
        csvData.add(String.valueOf(singleJuncFragments));
        csvData.add(String.valueOf(mReadDepth[mStreamIndices[FS_UPSTREAM]]));
        csvData.add(String.valueOf(mReadDepth[mStreamIndices[FS_DOWNSTREAM]]));

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

        return csvData.toString();
    }

    public static boolean lowerChromosome(final String chr, final String otherChr)
    {
        return chromosomeRank(chr) < chromosomeRank(otherChr);
    }

    public static int chromosomeRank(final String chromosome)
    {
        if(!HumanChromosome.contains(chromosome))
            return -1;

        if(chromosome.equals("X"))
            return 23;
        else if(chromosome.equals("Y"))
            return 24;
        else if(chromosome.equals("MT"))
            return 25;
        else
            return Integer.parseInt(chromosome);
    }

    public static String formLocationPair(final String[] chromosomes, final int[] geneCollectionIds)
    {
        return String.format("%s:%d_%s:%d",
                chromosomes[SE_START], geneCollectionIds[SE_START], chromosomes[SE_END], geneCollectionIds[SE_END]);
    }

    public static String formChromosomePair(final String chr1, final String chr2) { return chr1 + "_" + chr2; }
    public static String[] getChromosomePair(final String chrPair) { return chrPair.split("_"); }


}
