package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RnaUtils.impliedSvType;
import static com.hartwig.hmftools.isofox.fusion.FusionFragment.validPositions;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.BOTH_JUNCTIONS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
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
    private final byte[] mSjOrientations;

    private final List<List<EnsemblGeneData>> mCandidateGenes; // up and downstream genes
    private final String[] mFusionGeneIds;
    private final int[] mFusionIndices; // mapping of up & down stream to position data which is in SV terms

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

        mFragments = Lists.newArrayList(fragment);

        mLocationId = formLocationPair(mChromosomes, mGeneCollections);

        mRelatedFusions = Lists.newArrayList();
        mFusionGeneIds = new String[] {"", ""};
        mFusionIndices = new int[] {-1, -1};

        mCandidateGenes = Lists.newArrayListWithCapacity(FS_PAIR);
        mCandidateGenes.add(Lists.newArrayList());
        mCandidateGenes.add(Lists.newArrayList());
        mIncompleteData = false;
    }

    public int id() { return mId; }
    public String locationId() { return mLocationId; }
    public final String[] chromosomes() { return mChromosomes; }
    public final int[] geneCollections() { return mGeneCollections; }
    public final long[] splicePositions() { return mSjPositions; }
    public final byte[] spliceOrientations() { return mSjOrientations; }

    public boolean hasIncompleteData() { return mIncompleteData; }
    public void setIncompleteData() { mIncompleteData = true; }

    public boolean hasSplicedFragments() { return mFragments.stream().anyMatch(x -> x.isSpliced()); }
    public boolean hasUnsplicedFragments() { return mFragments.stream().anyMatch(x -> x.isUnspliced()); }

    public void setStreamData(final List<EnsemblGeneData> upstreamGenes, final List<EnsemblGeneData> downstreamGenes, boolean startIsUpstream)
    {
        mFusionIndices[FS_UPSTREAM] = startIsUpstream ? SE_START : SE_END;
        mFusionIndices[FS_DOWNSTREAM] = startIsUpstream ? SE_END : SE_START;
        mCandidateGenes.set(FS_UPSTREAM, upstreamGenes);
        mCandidateGenes.set(FS_DOWNSTREAM, downstreamGenes);

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
    public void addFusionFragment(final FusionFragment fragment) { mFragments.add(fragment); }

    public boolean spliceJunctionMatch(final FusionFragment fragment)
    {
        return validPositions(fragment.splicePositions()) && validPositions(mSjPositions)
                && mSjPositions[SE_START] == fragment.splicePositions()[SE_START] && mSjPositions[SE_END] == fragment.splicePositions()[SE_END]
                && mSjOrientations[SE_START] == fragment.spliceOrientations()[SE_START] && mSjOrientations[SE_END] == fragment.spliceOrientations()[SE_END];
    }

    public List<List<EnsemblGeneData>> getCandidateGenes() { return mCandidateGenes; }

    public boolean hasViableGenes() { return !mCandidateGenes.get(FS_UPSTREAM).isEmpty() && !mCandidateGenes.get(FS_DOWNSTREAM).isEmpty(); }
    public boolean hasValidStreamData() { return mFusionIndices[FS_UPSTREAM] >= 0 && mFusionIndices[FS_DOWNSTREAM] >= 0; }

    public boolean matchesTranscriptExons(final FusionReadData other)
    {
        final FusionFragment thisFragment = mFragments.get(0);
        final FusionFragment otherFragment = other.getFragments().get(0);

        for(int se = SE_START; se <= SE_END; ++se)
        {
            boolean hasMatch = false;

            for(List<TransExonRef> thisList : thisFragment.getTransExonRefs().get(se).values())
            {
                for(List<TransExonRef> otherList : otherFragment.getTransExonRefs().get(se).values())
                {
                    if(thisList.stream().anyMatch(x -> otherList.stream().anyMatch(y -> x.matches(y))))
                    {
                        hasMatch = true;
                        break;
                    }
                }

                if(hasMatch)
                    break;
            }

            if(!hasMatch)
                return false;
        }

        return true;
    }

    public void mergeFusionData(final FusionReadData other)
    {
        other.getFragments().forEach(x -> mFragments.add(x));
        other.getRelatedFusions().forEach(x -> mRelatedFusions.add(x));
    }

    public String getGeneName(int stream)
    {
        if(mCandidateGenes.get(stream).isEmpty())
            return "";

        if(mFusionGeneIds[stream].isEmpty())
            return mCandidateGenes.get(stream).get(0).GeneName;

        return mCandidateGenes.get(stream).stream()
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
        return "FusionId,Valid,GeneIdUp,GeneNameUp,ChromosomeUp,PositionUp,OrientUp,StrandUp"
                + ",GeneIdDown,GeneNameDown,ChromosomeDown,PositionDown,OrientDown,StrandDown"
                + ",SvType,TotalFragments,SpliceFragments,UnspliceFragments,DiscordantFragments"
                + ",TransDataUp,TransDataDown,OtherGenesUp,OtherGenesDown,RelatedFusions";
    }

    public String toCsv()
    {
        StringJoiner csvData = new StringJoiner(DELIMITER);

        csvData.add(String.format("Id_%d", mId));
        csvData.add(String.valueOf(hasViableGenes() && hasValidStreamData() && !hasIncompleteData()));

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            final String geneId = mFusionGeneIds[fs];
            final List<EnsemblGeneData> genes = mCandidateGenes.get(fs);

            csvData.add(geneId);

            final EnsemblGeneData geneData = genes.stream()
                    .filter(x -> x.GeneId.equals(geneId)).findFirst().map(x -> x).orElse(null);

            if(geneData != null)
            {
                csvData.add(geneData.GeneName);

                final int[] streamIndices = hasValidStreamData() ? mFusionIndices : new int[] { SE_START, SE_END };
                csvData.add(mChromosomes[streamIndices[fs]]);
                csvData.add(String.valueOf(mSjPositions[streamIndices[fs]]));
                csvData.add(String.valueOf(mSjOrientations[streamIndices[fs]]));
                csvData.add(String.valueOf(geneData.Strand));
            }
            else
            {
                csvData.add("");
                csvData.add(mChromosomes[fs]);
                csvData.add(String.valueOf(mSjPositions[fs]));
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

        csvData.add("TransUp");
        csvData.add("TransDowm");

        if(hasViableGenes())
        {
            String[] otherGenes = new String[] {"", ""};

            for (int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
            {
                for (final EnsemblGeneData geneData : mCandidateGenes.get(fs))
                {
                    if (!geneData.GeneId.equals(mFusionGeneIds[fs]))
                    {
                        otherGenes[fs] = appendStr(otherGenes[fs], geneData.GeneName, ';');
                    }
                }

                csvData.add(otherGenes[fs]);
            }
        }
        else
        {
            csvData.add("NONE");
            csvData.add("NONE");
        }

        if(!mRelatedFusions.isEmpty())
        {
            List<String> relatedFusions = mRelatedFusions.stream().map(x -> String.format("Id_%d", x)).collect(Collectors.toList());
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
                chromosomes[SE_START], geneCollectionIds[SE_START], chromosomes[SE_START], geneCollectionIds[SE_END]);
    }

    public static String formChromosomePair(final String chr1, final String chr2) { return chr1 + "_" + chr2; }
    public static String[] getChromosomePair(final String chrPair) { return chrPair.split("_"); }


}
