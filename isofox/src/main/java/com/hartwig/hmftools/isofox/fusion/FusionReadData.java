package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_MATCH;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.fusion.FusionFragment.validPositions;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.SPLICED_BOTH;
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
    private final List<FusionFragment> mFragments;

    private final List<Integer> mRelatedFusions;

    private final String[] mChromosomes;
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
        mSjPositions = new long[] { fragment.splicePositions()[SE_START], fragment.splicePositions()[SE_END] };
        mSjOrientations = new byte[]{ fragment.spliceOrientations()[SE_START], fragment.spliceOrientations()[SE_END] };

        mFragments = Lists.newArrayList(fragment);

        mRelatedFusions = Lists.newArrayList();
        mFusionGeneIds = new String[] {"", ""};
        mFusionIndices = new int[] {-1, -1};

        mCandidateGenes = Lists.newArrayListWithCapacity(FS_PAIR);
        mCandidateGenes.add(Lists.newArrayList());
        mCandidateGenes.add(Lists.newArrayList());
    }

    public int id() { return mId; }
    public final String[] chromosomes() { return mChromosomes; }
    public final long[] splicePositions() { return mSjPositions; }
    public final byte[] spliceOrientations() { return mSjOrientations; }

    public void setStreamData(final List<EnsemblGeneData> upstreamGenes, final List<EnsemblGeneData> downstreamGenes, boolean startIsUpstream)
    {
        mFusionIndices[FS_UPSTREAM] = startIsUpstream ? SE_START : SE_END;
        mCandidateGenes.set(FS_UPSTREAM, upstreamGenes);
        mCandidateGenes.set(FS_DOWNSTREAM, downstreamGenes);
    }

    public String chrPair() { return formChromosomePair(mChromosomes[SE_START], mChromosomes[SE_END]); }

    public final List<Integer> getRelatedFusions() { return mRelatedFusions; }

    public void addRelatedFusion(int id)
    {
        if(!mRelatedFusions.contains(id))
            mRelatedFusions.add(id);
    }

    public StructuralVariantType impliedSvType()
    {
        if(mChromosomes[SE_START].equals(mChromosomes[SE_END]))
        {
            if(mSjOrientations[SE_START] == mSjOrientations[SE_END])
            {
                return INV;
            }
            else if(mSjOrientations[SE_START] == 1)
            {
                return DEL;
            }
            else
            {
                return DUP;
            }
        }
        else
        {
            return BND;
        }
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
        return String.format("%d: chr(%s-%s) sj(%d-%d %d-%d %s) genes(%s-%s) frags(%d)",
                mId, mChromosomes[SE_START], mChromosomes[SE_END], mSjPositions[SE_START], mSjPositions[SE_END],
                mSjOrientations[SE_START], mSjOrientations[SE_END], impliedSvType(),
                getGeneName(FS_UPSTREAM), getGeneName(FS_DOWNSTREAM), mFragments.size());
    }

    public static String csvHeader()
    {
        return "FusionId,Valid,GeneIdUp,GeneNameUp,ChromosomeUp,PositionUp,OrientUp,StrandUp"
                + ",GeneIdDown,GeneNameDown,ChromosomeDown,PositionDown,OrientDown,StrandDown"
                + ",SvType,TotalFragments,SpliceFragments,DiscordantFragments"
                + ",TransDataUp,TransDataDown,OtherGenesUp,OtherGenesDown,RelatedFusions";
    }

    public String toCsv()
    {
        StringJoiner csvData = new StringJoiner(DELIMITER);

        csvData.add(String.format("Id_%d", mId));
        csvData.add(String.valueOf(hasViableGenes() && hasValidStreamData()));

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            final String geneId = mFusionGeneIds[fs];
            final List<EnsemblGeneData> genes = mCandidateGenes.get(fs);

            csvData.add(geneId);

            final String geneName = genes.stream()
                    .filter(x -> x.GeneId.equals(geneId)).findFirst().map(x -> x.GeneName).orElse("");

            csvData.add(geneName);

            final int[] streamIndices = hasValidStreamData() ? mFusionIndices : new int[] { SE_START, SE_END };
            csvData.add(mChromosomes[streamIndices[fs]]);
            csvData.add(String.valueOf(mSjPositions[streamIndices[fs]]));
            csvData.add(String.valueOf(mSjOrientations[streamIndices[fs]]));
        }

        csvData.add(impliedSvType().toString());

        int totalFragments = 0;
        int spliceFragments = 0;
        int discordantFragments = 0;
        for(final FusionFragment fragment : mFragments)
        {
            ++totalFragments;

            if(fragment.type() == SPLICED_BOTH)
                ++spliceFragments;
            else if(fragment.type() == DISCORDANT)
                ++discordantFragments;
        }

        csvData.add(String.valueOf(totalFragments));
        csvData.add(String.valueOf(spliceFragments));
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
                    if (!geneData.GeneId.equals(mFusionGeneIds[FS_UPSTREAM]))
                    {
                        otherGenes[fs] = appendStr(otherGenes[fs], geneData.GeneName, ';');
                    }
                }

                csvData.add(otherGenes[fs]);
            }
        }

        csvData.add("NONE");
        csvData.add("NONE");

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

    public static String formChromosomePair(final String chr1, final String chr2) { return chr1 + "_" + chr2; }
    public static String[] getChromosomePair(final String chrPair) { return chrPair.split("_"); }


}
