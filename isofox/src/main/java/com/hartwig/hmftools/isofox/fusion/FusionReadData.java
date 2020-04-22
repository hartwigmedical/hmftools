package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.getHighestMatchType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragment.validPositions;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionReadData
{
    private final int mId;
    private final String mLocationId;

    private final List<FusionFragment> mFragments;

    private boolean mIncompleteData;

    private final List<Integer> mRelatedFusions;

    private final String[] mChromosomes;
    private final int[] mGeneCollections;
    private final long[] mSjPositions;
    private long[] mUnsplicedPositions; // used when this fusion has both spliced and unspliced fragments
    private final byte[] mSjOrientations;

    private final List<EnsemblGeneData>[] mCandidateGenes; // up and downstream genes
    private final String[] mFusionGeneIds; // stored by stream
    private final int[] mFusionIndices; // mapping of up & down stream to position data which is in SV terms
    private final List<List<TransExonRef>> mTransExonRefs;

    public static final int FS_UPSTREAM = 0;
    public static final int FS_DOWNSTREAM = 1;
    public static final int FS_PAIR = 2;

    public FusionReadData(int id, final FusionFragment fragment)
    {
        mId = id;

        mChromosomes = new String[] { fragment.chromosomes()[SE_START], fragment.chromosomes()[SE_END] };
        mGeneCollections = new int[] { fragment.geneCollections()[SE_START], fragment.geneCollections()[SE_END] };
        mSjPositions = new long[] { fragment.splicePositions()[SE_START], fragment.splicePositions()[SE_END] };
        mSjOrientations = new byte[]{ fragment.spliceOrientations()[SE_START], fragment.spliceOrientations()[SE_END] };

        mUnsplicedPositions = null;

        mFragments = Lists.newArrayList();
        addFusionFragment(fragment);

        mLocationId = formLocationPair(mChromosomes, mGeneCollections);

        mRelatedFusions = Lists.newArrayList();
        mFusionGeneIds = new String[] {"", ""};
        mFusionIndices = new int[] {-1, -1};

        mCandidateGenes = new List[FS_PAIR];
        mCandidateGenes[SE_START] = Lists.newArrayList();
        mCandidateGenes[SE_END] = Lists.newArrayList();
        mIncompleteData = false;

        mTransExonRefs = Lists.newArrayListWithCapacity(FS_PAIR);
        mTransExonRefs.add(Lists.newArrayList());
        mTransExonRefs.add(Lists.newArrayList());

        List<String>[] tmp = new List[SE_PAIR];
    }

    public int id() { return mId; }
    public String locationId() { return mLocationId; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final long[] splicePositions() { return mSjPositions; }
    public final byte[] spliceOrientations() { return mSjOrientations; }

    public boolean hasIncompleteData() { return mIncompleteData; }
    public void setIncompleteData() { mIncompleteData = true; }

    public List<TransExonRef> getTransExonRefsByPos(int se) { return mTransExonRefs.get(se); }
    public List<TransExonRef> getTransExonRefsByStream(int fs) { return mTransExonRefs.get(mFusionIndices[fs]); }

    public boolean hasSplicedFragments() { return mFragments.stream().anyMatch(x -> x.isSpliced()); }
    public boolean hasUnsplicedFragments() { return mFragments.stream().anyMatch(x -> x.isUnspliced()); }

    public void setStreamData(final List<EnsemblGeneData> upstreamGenes, final List<EnsemblGeneData> downstreamGenes, boolean startIsUpstream)
    {
        mFusionIndices[FS_UPSTREAM] = startIsUpstream ? SE_START : SE_END;
        mFusionIndices[FS_DOWNSTREAM] = startIsUpstream ? SE_END : SE_START;
        mCandidateGenes[FS_UPSTREAM] = upstreamGenes;
        mCandidateGenes[FS_DOWNSTREAM] = downstreamGenes;

        // until a more informed decision can be made
        mFusionGeneIds[FS_UPSTREAM] = upstreamGenes.get(0).GeneId;
        mFusionGeneIds[FS_DOWNSTREAM] = downstreamGenes.get(0).GeneId;
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
        return impliedSvType(mChromosomes, mSjOrientations);
    }

    public final List<FusionFragment> getFragments() { return mFragments; }
    public void addFusionFragment(final FusionFragment fragment)
    {
        mFragments.add(fragment);
    }

    public boolean spliceJunctionMatch(final FusionFragment fragment)
    {
        return validPositions(fragment.splicePositions()) && validPositions(mSjPositions)
                && mSjPositions[SE_START] == fragment.splicePositions()[SE_START] && mSjPositions[SE_END] == fragment.splicePositions()[SE_END]
                && mSjOrientations[SE_START] == fragment.spliceOrientations()[SE_START] && mSjOrientations[SE_END] == fragment.spliceOrientations()[SE_END];
    }

    public List<EnsemblGeneData>[] getCandidateGenes() { return mCandidateGenes; }

    public boolean hasViableGenes() { return !mCandidateGenes[FS_UPSTREAM].isEmpty() && !mCandidateGenes[FS_DOWNSTREAM].isEmpty(); }
    public boolean hasValidStreamData() { return mFusionIndices[FS_UPSTREAM] >= 0 && mFusionIndices[FS_DOWNSTREAM] >= 0; }
    public int streamStartEnd(int fs) { return mFusionIndices[fs] >= 0 ? mFusionIndices[fs] : SE_START; }

    public boolean isValid() { return hasViableGenes() && hasValidStreamData() && !hasIncompleteData(); }

    public static boolean hasTranscriptExonMatch(final List<TransExonRef> list1, final List<TransExonRef> list2)
    {
        return list1.stream().anyMatch(x -> list2.stream().anyMatch(y -> x.matches(y)));
    }

    public static boolean hasTranscriptNextExonMatch(final List<TransExonRef> list1, final List<TransExonRef> list2)
    {
        return list1.stream().anyMatch(x -> list2.stream().anyMatch(y -> x.matchesNext(y)));
    }

    public void mergeFusionData(final FusionReadData other)
    {
        other.getFragments().forEach(x -> mFragments.add(x));
        other.getRelatedFusions().forEach(x -> mRelatedFusions.add(x));

        if(other.hasUnsplicedFragments())
        {
            final FusionFragment fragment = other.getFragments().stream().filter(x -> x.isUnspliced()).findFirst().orElse(null);
            mUnsplicedPositions = new long[] { fragment.splicePositions()[SE_START], fragment.splicePositions()[SE_END] };
        }
    }

    public void cacheTranscriptData()
    {
        final FusionFragment fragment = mFragments.get(0);

        for (int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            final Map<RegionMatchType, List<TransExonRef>> transExonRefs = fragment.getTransExonRefs().get(fs);
            RegionMatchType highestMatchType = getHighestMatchType(transExonRefs.keySet());
            List<TransExonRef> teData = transExonRefs.get(highestMatchType);

            if(teData != null)
                mTransExonRefs.get(fs).addAll(teData);
        }
    }

    public boolean canAddDiscordantFragment(final FusionFragment fragment)
    {
        // the 2 reads' bounds need to fall within a correct intron relative to the SJs
        // and additionally be bound by any unspliced SJ positions
        if(mUnsplicedPositions != null)
        {
            for(ReadRecord read : fragment.getReads())
            {
                long[] readBounds = { read.getCoordsBoundary(true), read.getCoordsBoundary(false) };

                if(read.Chromosome.equals(mChromosomes[SE_START]) && read.getGeneCollecton() == mGeneCollections[SE_START])
                {
                    long minBound = min(mSjPositions[SE_START], mUnsplicedPositions[SE_START]);
                    long maxBound = max(mSjPositions[SE_START], mUnsplicedPositions[SE_START]);
                    if(!positionsWithin(readBounds[SE_START], readBounds[SE_END], minBound, maxBound))
                        return false;
                }
                else
                {
                    long minBound = min(mSjPositions[SE_END], mUnsplicedPositions[SE_END]);
                    long maxBound = max(mSjPositions[SE_END], mUnsplicedPositions[SE_END]);
                    if(!positionsWithin(readBounds[SE_START], readBounds[SE_END], minBound, maxBound))
                        return false;
                }
            }

            return true;
        }

        // otherwise check whether the correct exon is closest for each applicable transcript
        final List<TransExonRef> startTrans = fragment.getTransExonRefs().get(SE_START).get(INTRON);
        final List<TransExonRef> endTrans = fragment.getTransExonRefs().get(SE_END).get(INTRON);

        if(startTrans.isEmpty() || endTrans.isEmpty())
            return false;

        if(hasTranscriptExonMatch(startTrans, getTransExonRefsByStream(FS_UPSTREAM))
        && hasTranscriptNextExonMatch(endTrans, getTransExonRefsByStream(FS_DOWNSTREAM)))
        {
            return true;
        }
        else if(hasTranscriptExonMatch(endTrans, getTransExonRefsByStream(FS_UPSTREAM))
        && hasTranscriptNextExonMatch(startTrans, getTransExonRefsByStream(FS_DOWNSTREAM)))
        {
            return true;
        }

        return false;
    }

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
        return String.format("%d: chr(%s-%s) sj(%d-%d %d/%d %s) genes(%s-%s) frags(%d)",
                mId, mChromosomes[SE_START], mChromosomes[SE_END], mSjPositions[SE_START], mSjPositions[SE_END],
                mSjOrientations[SE_START], mSjOrientations[SE_END], getImpliedSvType(),
                getGeneName(FS_UPSTREAM), getGeneName(FS_DOWNSTREAM), mFragments.size());
    }

    public static String csvHeader()
    {
        return "FusionId,Valid,GeneIdUp,GeneNameUp,ChromosomeUp,PositionUp,UnsplicedPosUp,OrientUp,StrandUp"
                + ",GeneIdDown,GeneNameDown,ChromosomeDown,PositionDown,UnsplicedPosDown,OrientDown,StrandDown"
                + ",SvType,TotalFragments,SpliceFragments,UnspliceFragments,DiscordantFragments"
                + ",TransDataUp,TransDataDown,OtherGenesUp,OtherGenesDown,RelatedFusions";
    }

    public static String fusionId(int id) { return String.format("Id_%d", id); }

    public String toCsv()
    {
        StringJoiner csvData = new StringJoiner(DELIMITER);

        csvData.add(fusionId(mId));
        csvData.add(String.valueOf(hasViableGenes() && hasValidStreamData() && !hasIncompleteData()));

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

                final int[] streamIndices = hasValidStreamData() ? mFusionIndices : new int[] { SE_START, SE_END };
                csvData.add(mChromosomes[streamIndices[fs]]);
                csvData.add(String.valueOf(mSjPositions[streamIndices[fs]]));
                csvData.add(mUnsplicedPositions != null ? String.valueOf(mUnsplicedPositions[streamIndices[fs]]) : "-1");
                csvData.add(String.valueOf(mSjOrientations[streamIndices[fs]]));
                csvData.add(String.valueOf(geneData.Strand));
            }
            else
            {
                csvData.add("");
                csvData.add(mChromosomes[fs]);
                csvData.add(String.valueOf(mSjPositions[fs]));
                csvData.add("-1");
                csvData.add(String.valueOf(mSjOrientations[fs]));
                csvData.add("0");
            }
        }

        csvData.add(getImpliedSvType().toString());

        int totalFragments = 0;
        int splicedFragments = 0;
        int unsplicedFragments = 0;
        int discordantFragments = 0;
        for(final FusionFragment fragment : mFragments)
        {
            ++totalFragments;

            if(fragment.isSpliced())
                ++splicedFragments;
            else if(fragment.isUnspliced())
                ++unsplicedFragments;
            else if(fragment.type() == DISCORDANT)
                ++discordantFragments;
        }

        csvData.add(String.valueOf(totalFragments));
        csvData.add(String.valueOf(splicedFragments));
        csvData.add(String.valueOf(unsplicedFragments));
        csvData.add(String.valueOf(discordantFragments));

        for (int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            if(mTransExonRefs.get(fs).isEmpty())
            {
                csvData.add("NONE");
                continue;
            }

            String transData = "";
            for(final TransExonRef transExonRef : mTransExonRefs.get(fs))
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
        if(chromosome.equals("X"))
            return 23;
        else if(chromosome.equals("Y"))
            return 24;
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
