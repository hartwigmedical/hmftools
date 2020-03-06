package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.FragmentMatchType.MAX_FRAG_TYPE;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.typeAsInt;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.deriveCommonRegions;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionsOverlap;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GeneReadData
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions
    private final List<long[]> mCommonExonicRegions; // merge any overlapping exons, to form a set of exonic regions for the gene

    private final List<TranscriptData> mTranscripts;
    private final long[] mTranscriptsRange;

    private final List<long[]> mOtherGeneExonicRegions; // a set of exons from other genes which overlap this gene

    // summary results
    private final Map<Integer,int[][]> mTranscriptReadCounts; // count of fragments support types for each transcript, and whether unique
    private final Map<String,Double> mTranscriptAllocations; // results from the expected rate vs counts fit routine
    private double mFitResiduals;

    private final int[] mFragmentCounts;

    private static final Logger LOGGER = LogManager.getLogger(GeneReadData.class);

    public GeneReadData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;

        mExonRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
        mCommonExonicRegions = Lists.newArrayList();
        mOtherGeneExonicRegions = Lists.newArrayList();

        mFragmentCounts = new int[typeAsInt(GeneMatchType.MAX)];
        mTranscriptReadCounts = Maps.newHashMap();
        mTranscriptAllocations = Maps.newHashMap();
        mFitResiduals = 0;

        mTranscriptsRange = new long[SE_PAIR];
    }

    public String name() { return GeneData.GeneName;}

    public final List<TranscriptData> getTranscripts() { return mTranscripts; }
    public void setTranscripts(final List<TranscriptData> transDataList) { mTranscripts.addAll(transDataList); }

    public final List<RegionReadData> getExonRegions() { return mExonRegions; }

    public void addExonRegion(final RegionReadData exonRegionData)
    {
        if(hasExonRegion(exonRegionData.Region.start(), exonRegionData.Region.end()))
            return;

        mExonRegions.add(exonRegionData);
    }

    public boolean hasExonRegion(long posStart, long posEnd)
    {
        return mExonRegions.stream().anyMatch(x -> x.Region.start() == posStart && x.Region.end() == posEnd);
    }

    public RegionReadData findExonRegion(long posStart, long posEnd)
    {
        return mExonRegions.stream()
                .filter(x -> x.Region.start() == posStart && x.Region.end() == posEnd)
                .findFirst().orElse(null);
    }

    public final long[] getTranscriptsRange() { return mTranscriptsRange; }

    public void generateExonicRegions()
    {
        // form a genomic region for each unique exon amongst the transcripts
        for(final TranscriptData transData : mTranscripts)
        {
            RegionReadData prevRegionReadData = null;

            for(int i = 0; i < transData.exons().size(); ++ i)
            {
                ExonData exon = transData.exons().get(i);

                RegionReadData exonReadData = findExonRegion(exon.ExonStart, exon.ExonEnd);

                if (exonReadData == null)
                {
                    GenomeRegion region = GenomeRegions.create(GeneData.Chromosome, exon.ExonStart, exon.ExonEnd);
                    exonReadData = new RegionReadData(region);
                    addExonRegion(exonReadData);
                }

                exonReadData.addExonRef(transData.TransId, transData.TransName, exon.ExonRank);

                if(prevRegionReadData != null)
                {
                    prevRegionReadData.addPostRegion(exonReadData);
                    exonReadData.addPreRegion(prevRegionReadData);
                }

                prevRegionReadData = exonReadData;
            }

            mTranscriptsRange[SE_END] = max(transData.TransEnd, mTranscriptsRange[SE_END]);

            if(mTranscriptsRange[SE_START] == 0 || transData.TransStart < mTranscriptsRange[SE_START])
                mTranscriptsRange[SE_START] = transData.TransStart;
        }

        generateCommonExonicRegions();
    }

    public final int[] getCounts() { return mFragmentCounts; }
    public void addCount(GeneMatchType type, int count) { mFragmentCounts[typeAsInt(type)] += count; }

    public static final int TRANS_COUNT = 0;
    public static final int UNIQUE_TRANS_COUNT = 1;

    public int[][] getTranscriptReadCount(final String trans)
    {
        int[][] counts = mTranscriptReadCounts.get(trans);
        return counts != null ? counts : new int[MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
    }

    public void addTranscriptReadMatch(int transId, boolean isUnique, FragmentMatchType type)
    {
        int[][] counts = mTranscriptReadCounts.get(transId);
        if(counts == null)
        {
            counts = new int[MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
            mTranscriptReadCounts.put(transId,  counts);
        }

        if(isUnique)
        {
            ++counts[FragmentMatchType.typeAsInt(type)][UNIQUE_TRANS_COUNT];
        }

        ++counts[FragmentMatchType.typeAsInt(type)][TRANS_COUNT];
    }

    public Map<String,Double> getTranscriptAllocations() { return mTranscriptAllocations; }

    public void setFitResiduals(double residuals) { mFitResiduals = residuals; }
    public double getFitResiduals() { return mFitResiduals; }

    public double getTranscriptAllocation(final String transName)
    {
        Double allocation = mTranscriptAllocations.get(transName);
        return allocation != null ? allocation : 0;
    }

    private void generateCommonExonicRegions()
    {
        if(mExonRegions.isEmpty())
            return;

        List<long[]> commonExonicRegions = Lists.newArrayList(new long[] {mExonRegions.get(0).start(), mExonRegions.get(0).end()});

        for(int i = 1; i < mExonRegions.size(); ++i)
        {
            List<long[]> nextRegion = Lists.newArrayList(new long[] {mExonRegions.get(i).start(), mExonRegions.get(i).end()});
            commonExonicRegions = deriveCommonRegions(commonExonicRegions, nextRegion);
        }

        mCommonExonicRegions.addAll(commonExonicRegions);
    }

    public List<long[]> getCommonExonicRegions() { return mCommonExonicRegions; }
    public List<long[]> getOtherGeneExonicRegions() { return mOtherGeneExonicRegions; }

    public boolean overlapsOtherGeneExon(long posStart, long posEnd)
    {
        return mOtherGeneExonicRegions.stream().anyMatch(x -> positionsOverlap(posStart, posEnd, x[SE_START], x[SE_END]));
    }

    public long calcExonicRegionLength()
    {
        return mCommonExonicRegions.stream().mapToLong(x -> x[SE_END] - x[SE_START]).sum();
    }

    public static void markOverlappingGeneRegions(final List<GeneReadData> geneReadDataList, boolean logExonOverlaps)
    {
        // record against each gene any exon from another gene which overlaps it
        for(GeneReadData geneReadData : geneReadDataList)
        {
            long geneStart = geneReadData.GeneData.GeneStart;
            long geneEnd = geneReadData.GeneData.GeneEnd;

            int geneOverlaps = 0;
            for(int i = 0; i < geneReadDataList.size(); ++i)
            {
                final GeneReadData otherGeneReadData = geneReadDataList.get(i);

                if(otherGeneReadData == geneReadData)
                    continue;

                if(otherGeneReadData.GeneData.GeneStart > geneEnd)
                    break;

                if(otherGeneReadData.GeneData.GeneEnd < geneStart)
                    continue;

                geneReadData.getOtherGeneExonicRegions().addAll(otherGeneReadData.getCommonExonicRegions());
                ++geneOverlaps;

                if(logExonOverlaps)
                    logGeneExonOverlaps(geneReadData, otherGeneReadData);
            }

            /*
            if(geneOverlaps > 3)
            {
                LOGGER.info("gene({}) has {} overlaps and {} other-exon regions",
                        geneReadData.name(), geneOverlaps, geneReadData.getOtherGeneExonicRegions().size());
            }
           */
        }
    }

    public List<RegionReadData> findOverlappingRegions(final ReadRecord read)
    {
        return mExonRegions.stream()
                .filter(x -> read.overlapsMappedReads(x.Region.start(), x.Region.end()))
                .collect(Collectors.toList());
    }

    public static void logGeneExonOverlaps(final GeneReadData gene1, final GeneReadData gene2)
    {
        for(final long[] region1 : gene1.getCommonExonicRegions())
        {
            for(final long[] region2 : gene2.getCommonExonicRegions())
            {
                if (region1[SE_START] > region2[SE_END] || region1[SE_END] < region2[SE_START])
                    continue;

                // Time,Chromosome,GeneId1,GeneName1,PosStart1,PosEnd1,GeneId2,GeneName2,PosStart2,PosEnd2
                LOGGER.info("GENE_EXON_OVERLAP: {},{},{},{},{},{},{},{},{}",
                        gene1.GeneData.Chromosome, gene1.GeneData.GeneId, gene1.GeneData.GeneName, region1[SE_START], region1[SE_END],
                        gene2.GeneData.GeneId, gene2.GeneData.GeneName, region2[SE_START], region2[SE_END]);
            }
        }
    }

    public String toString()
    {
        return String.format("%s:%s location(%s:%d -> %d) trans(%d)",
                GeneData.GeneId, GeneData.GeneName, GeneData.Chromosome, GeneData.GeneStart, GeneData.GeneEnd, mTranscripts.size());
    }

    @VisibleForTesting
    public void clearCounts()
    {
        for(int i = 0; i < mFragmentCounts.length; ++i)
            mFragmentCounts[i] = 0;
    }

}
