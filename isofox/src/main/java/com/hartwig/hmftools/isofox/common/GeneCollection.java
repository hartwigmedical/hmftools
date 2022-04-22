package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.FragmentType.typeAsInt;
import static com.hartwig.hmftools.isofox.common.GeneReadData.generateCommonExonicRegions;
import static com.hartwig.hmftools.isofox.common.RegionReadData.findExonRegion;
import static com.hartwig.hmftools.isofox.common.RegionReadData.generateExonicRegions;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class GeneCollection
{
    private final int mId;
    private final List<GeneReadData> mGenes;
    private final List<String> mGeneIds;
    private final String mChromosome;
    private final int[] mRegionBounds;
    private final int[] mNonGenicPositions; // position prior to this region beginning and since the last gene collection
    private boolean mEndOfChromosome;

    private final Map<Integer,GeneReadData> mTransIdsGeneMap;

    private final List<RegionReadData> mExonRegions; // set of unique exons ie with differing start and end positions
    private final List<int[]> mCommonExonicRegions; // merge any overlapping exons, to form a set of exonic regions for the gene
    private final List<TranscriptData> mTranscripts;

    private List<TranscriptData> mEnrichedTranscripts;
    private int[] mEnrichedRegion; // special regions of high read density

    // summary results
    private final Map<Integer,int[][]> mTranscriptReadCounts; // count of fragments support types for each transcript, and whether unique
    private final int[] mFragmentCounts; // counts by various classifications

    public GeneCollection(int id, final List<GeneReadData> genes)
    {
        mId = id;
        mGenes = genes;

        mGeneIds = Lists.newArrayList();
        mGenes.forEach(x -> mGeneIds.add(x.GeneData.GeneId));

        mRegionBounds = new int[SE_PAIR];
        mNonGenicPositions = new int[SE_PAIR];
        mEndOfChromosome = false;

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
    public final int[] regionBounds() { return mRegionBounds; }
    public void setNonGenicPosition(int se, int position) { mNonGenicPositions[se] = position; }
    public int[] getNonGenicPositions() { return mNonGenicPositions; }
    public boolean isEndOfChromosome() { return mEndOfChromosome; }
    public void setEndOfChromosome() { mEndOfChromosome = true; }

    public final List<TranscriptData> getTranscripts() { return mTranscripts; }

    public final List<RegionReadData> getExonRegions() { return mExonRegions; }
    public List<int[]> getCommonExonicRegions() { return mCommonExonicRegions; }

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
            geneNames.append(ITEM_DELIM + mGenes.get(i).name());
        }

        return geneNames.toString();
    }

    public boolean containsEnrichedRegion() { return mEnrichedRegion != null; }
    public List<TranscriptData> getEnrichedTranscripts() { return mEnrichedTranscripts; }
    public int[] getEnrichedRegion() { return mEnrichedRegion; }

    public boolean inEnrichedRegion(int posStart, int posEnd)
    {
        if(mEnrichedRegion == null)
            return false;

        return positionsOverlap(posStart, posEnd, mEnrichedRegion[SE_START], mEnrichedRegion[SE_END]);
    }

    public void markEnrichedAndExcludedGenes(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        if(config.Filters.EnrichedGeneIds.isEmpty())
            return;

        for(GeneReadData geneReadData : mGenes)
        {
            if(config.Filters.EnrichedGeneIds.contains(geneReadData.GeneData.GeneId))
            {
                mEnrichedTranscripts = Lists.newArrayList(geneTransCache.getTranscripts(geneReadData.GeneData.GeneId));
                mEnrichedRegion = new int[SE_PAIR];

                for(TranscriptData transData : mEnrichedTranscripts)
                {
                    for(ExonData exonData : transData.exons())
                    {
                        mEnrichedRegion[SE_START] = mEnrichedRegion[SE_START] > 0
                                ? min(mEnrichedRegion[SE_START], exonData.Start) : exonData.Start;

                        mEnrichedRegion[SE_END] = max(mEnrichedRegion[SE_END], exonData.End);
                    }
                }
            }
        }
    }

    public void setEnrichedTranscripts(final List<TranscriptData> transDataList)
    {
        mEnrichedTranscripts = Lists.newArrayList(transDataList);
        mEnrichedRegion = new int[SE_PAIR];

        for(TranscriptData transData : mEnrichedTranscripts)
        {
            for(ExonData exonData : transData.exons())
            {
                mEnrichedRegion[SE_START] = mEnrichedRegion[SE_START] > 0
                        ? min(mEnrichedRegion[SE_START], exonData.Start) : exonData.Start;

                mEnrichedRegion[SE_END] = max(mEnrichedRegion[SE_END], exonData.End);
            }
        }
    }

    public void setReadGeneCollections(final ReadRecord read, final int[] nonGenicBounds)
    {
        if(positionsWithin(read.PosStart, read.PosEnd, mRegionBounds[SE_START], mRegionBounds[SE_END]))
        {
            read.setGeneCollection(SE_START, mId, true);
            read.setGeneCollection(SE_END, mId, true);
        }
        else
        {
            // mark any read extending beyond this gene collection's bounds in part or full
            for (int se = SE_START; se <= SE_END; ++se)
            {
                if(positionWithin(
                        read.getCoordsBoundary(se), mRegionBounds[SE_START], mRegionBounds[SE_END]))
                {
                    read.setGeneCollection(se, mId, true);
                }
                else if(positionWithin(read.getCoordsBoundary(se), nonGenicBounds[SE_START], nonGenicBounds[SE_END]))
                {
                    read.setGeneCollection(se, mId, false);
                }

                // otherwise the gene collection will remain unset for now
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
                    final RegionReadData exonReadData = findExonRegion(mExonRegions, exon.Start, exon.End);
                    if (exonReadData == null)
                    {
                        ISF_LOGGER.error("genes({}) failed to create exonic regions", geneNames());
                        return;
                    }

                    gene.addExonRegion(exonReadData);
                }
            }

            gene.setHasUnsplicedRegions();
        }

        generateCommonExonicRegions(mExonRegions, mCommonExonicRegions);
    }

    public List<GeneReadData> findGenesCoveringRange(int posStart, int posEnd, boolean checkUnspliced)
    {
        if(checkUnspliced)
        {
            // check for the region being a) contained within the gene and b) the gene has unspliced sections
            return mGenes.stream()
                    .filter(x -> x.hasUnsplicedRegions())
                    .filter(x -> positionsWithin(posStart, posEnd, x.GeneData.GeneStart, x.GeneData.GeneEnd))
                    .collect(Collectors.toList());
        }
        else
        {
            // any overlap
            return mGenes.stream()
                    .filter(x -> positionsOverlap(x.GeneData.GeneStart, x.GeneData.GeneEnd, posStart, posEnd))
                    .collect(Collectors.toList());
        }
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

    public void addTranscriptReadMatch(int transId, FragmentMatchType type)
    {
        int[][] counts = mTranscriptReadCounts.get(transId);
        if(counts == null)
        {
            counts = new int[FragmentMatchType.MAX_FRAG_TYPE][UNIQUE_TRANS_COUNT+1];
            mTranscriptReadCounts.put(transId,  counts);
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

}
