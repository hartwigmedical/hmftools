package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.GeneReadData.generateCommonExonicRegions;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.isofox.common.RegionReadData.generateExonicRegions;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionsWithin;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class GeneCollection
{
    private final int mId;
    private final List<GeneReadData> mGenes;
    private final List<String> mGeneIds;
    private final String mChromosome;
    private final long[] mRegionBounds;

    private final Map<Integer,GeneReadData> mTransIdsGeneMap;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions
    private final List<long[]> mCommonExonicRegions; // merge any overlapping exons, to form a set of exonic regions for the gene
    private final List<TranscriptData> mTranscripts;

    private List<TranscriptData> mEnrichedTranscripts;
    private long[] mEnrichedRegion; // special regions of high read density

    // summary results
    private final Map<Integer,int[][]> mTranscriptReadCounts; // count of fragments support types for each transcript, and whether unique
    private final int[] mFragmentCounts; // counts by various classifications

    public GeneCollection(int id, final List<GeneReadData> genes)
    {
        mId = id;
        mGenes = genes;

        mGeneIds = Lists.newArrayList();
        mGenes.forEach(x -> mGeneIds.add(x.GeneData.GeneId));

        mRegionBounds = new long[SE_PAIR];

        mChromosome = genes.get(0).GeneData.Chromosome;

        mTransIdsGeneMap = Maps.newHashMap();

        mExonRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
        mCommonExonicRegions = Lists.newArrayList();

        buildCache();

        mTranscriptReadCounts = Maps.newHashMap();
        mFragmentCounts = new int[typeAsInt(FragmentType.MAX)];

        mEnrichedTranscripts = null;
        mEnrichedRegion = null;
    }

    public int id() { return mId; }
    public String chrId() { return String.format("%s_%d", mChromosome, mId); }
    public final String chromosome() { return mChromosome; }
    public final List<GeneReadData> genes() { return mGenes; }
    public final List<String> geneIds() { return mGeneIds; }
    public final long[] regionBounds() { return mRegionBounds; }

    public final List<TranscriptData> getTranscripts() { return mTranscripts; }

    public final List<RegionReadData> getExonRegions() { return mExonRegions; }
    public List<long[]> getCommonExonicRegions() { return mCommonExonicRegions; }

    public int getStrand(int transId)
    {
        final GeneReadData gene = mTransIdsGeneMap.get(transId);
        return gene != null ? gene.GeneData.Strand : 0;
    }

    public String geneNames() { return geneNames(10); }

    public String geneNames(int maxCount)
    {
        if(mGenes.size() == 1)
            return mGenes.get(0).name();

        StringBuilder geneNames = new StringBuilder(mGenes.get(0).name());

        for(int i = 1; i < min(mGenes.size(), maxCount); ++i)
        {
            geneNames.append(";" + mGenes.get(i).name());
        }

        return geneNames.toString();
    }

    public List<TranscriptData> getEnrichedTranscripts() { return mEnrichedTranscripts; }
    public long[] getEnrichedRegion() { return mEnrichedRegion; }

    public void setEnrichedTranscripts(final List<TranscriptData> transDataList)
    {
        mEnrichedTranscripts = Lists.newArrayList(transDataList);
        mEnrichedRegion = new long[SE_PAIR];

        for(TranscriptData transData : mEnrichedTranscripts)
        {
            for(ExonData exonData : transData.exons())
            {
                mEnrichedRegion[SE_START] = mEnrichedRegion[SE_START] > 0
                        ? min(mEnrichedRegion[SE_START], exonData.ExonStart) : exonData.ExonStart;

                mEnrichedRegion[SE_END] = max(mEnrichedRegion[SE_END], exonData.ExonEnd);
            }
        }
    }

    private void buildCache()
    {
        for(final GeneReadData gene : mGenes)
        {
            for(final TranscriptData transData : gene.getTranscripts())
            {
                mTransIdsGeneMap.put(transData.TransId, gene);
                mTranscripts.add(transData);

                mRegionBounds[SE_START] = mRegionBounds[SE_START] == 0 ? transData.TransStart : min(mRegionBounds[SE_START], transData.TransStart);
                mRegionBounds[SE_END] = max(mRegionBounds[SE_END], transData.TransEnd);
            }

            generateExonicRegions(gene.GeneData.GeneId, mChromosome, mExonRegions, gene.getTranscripts());

            // cache the relevant set of exon regions back into the gene for convenience
            for(final TranscriptData transData : gene.getTranscripts())
            {
                for (final ExonData exon : transData.exons())
                {
                    final RegionReadData exonReadData = findExonRegion(mExonRegions, exon.ExonStart, exon.ExonEnd);
                    if (exonReadData == null)
                    {
                        ISF_LOGGER.error("genes({}) failed to create exonic regions", geneNames());
                        return;
                    }

                    gene.addExonRegion(exonReadData);
                }
            }
        }

        generateCommonExonicRegions(mExonRegions, mCommonExonicRegions);
    }

    public boolean hasGeneData(int transId) { return mTransIdsGeneMap.containsKey(transId); }

    public GeneReadData findGeneData(int transId)
    {
        return mTransIdsGeneMap.get(transId);
    }

    public List<GeneReadData> findGenesCoveringRange(long posStart, long posEnd)
    {
        return mGenes.stream()
                .filter(x -> positionsOverlap(x.GeneData.GeneStart, x.GeneData.GeneEnd, posStart, posEnd))
                .collect(Collectors.toList());
    }

    public static final int TRANS_COUNT = 0;
    public static final int UNIQUE_TRANS_COUNT = 1;

    public int[][] getTranscriptReadCount(final int transId)
    {
        int[][] counts = mTranscriptReadCounts.get(transId);
        return counts != null ? counts : new int[FragmentMatchType.MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
    }

    public void addTranscriptReadMatch(int transId, boolean isUnique, FragmentMatchType type)
    {
        int[][] counts = mTranscriptReadCounts.get(transId);
        if(counts == null)
        {
            counts = new int[FragmentMatchType.MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
            mTranscriptReadCounts.put(transId,  counts);
        }

        if(isUnique)
        {
            ++counts[FragmentMatchType.typeAsInt(type)][UNIQUE_TRANS_COUNT];
        }

        ++counts[FragmentMatchType.typeAsInt(type)][TRANS_COUNT];
    }

    public final int[] getCounts() { return mFragmentCounts; }
    public void addCount(FragmentType type, int count) { mFragmentCounts[typeAsInt(type)] += count; }

    @VisibleForTesting
    public void clearCounts()
    {
        for(int i = 0; i < mFragmentCounts.length; ++i)
            mFragmentCounts[i] = 0;
    }

    public void logOverlappingGenes(final List<GeneReadData> overlappingGenes)
    {
        String geneNamesStr = "";
        int transcriptCount = 0;
        long minRange = -1;
        long maxRange = 0;

        for(GeneReadData geneReadData : overlappingGenes)
        {
            geneNamesStr = appendStr(geneNamesStr, geneReadData.GeneData.GeneId, ';');
            transcriptCount += geneReadData.getTranscripts().size();
            maxRange =  max(maxRange, geneReadData.GeneData.GeneEnd);
            minRange =  minRange < 0 ? geneReadData.GeneData.GeneStart : min(minRange, geneReadData.GeneData.GeneStart);
        }

        // Time,Chromosome,GeneCount,TranscriptCount,RangeStart,RangeEnd,GeneNames
        ISF_LOGGER.info("GENE_OVERLAP: {},{},{},{},{},{}", // chr({}) genes({}) transcripts({}) range({} -> {}),
                mChromosome, overlappingGenes.size(), transcriptCount, minRange, maxRange, geneNamesStr);
    }
}
