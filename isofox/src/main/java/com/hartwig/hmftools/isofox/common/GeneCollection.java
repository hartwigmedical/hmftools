package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

public class GeneCollection
{
    private final List<GeneReadData> mGenes;
    private final String mChromosome;
    private final long[] mRegionBounds;

    private final Map<Integer,GeneReadData> mTransIdsGeneMap;

    public GeneCollection(final List<GeneReadData> genes)
    {
        mGenes = genes;

        mRegionBounds = new long[SE_PAIR];

        mRegionBounds[SE_START] = genes.stream().mapToLong(x -> x.GeneData.GeneStart).min().orElse(0);
        mRegionBounds[SE_END] = genes.stream().mapToLong(x -> x.GeneData.GeneEnd).max().orElse(0);

        mChromosome = genes.get(0).GeneData.Chromosome;

        mTransIdsGeneMap = Maps.newHashMap();
        buildTransGeneMap();
    }

    public final List<GeneReadData> genes() { return mGenes; }
    public final long[] regionBounds() { return mRegionBounds; }

    private void buildTransGeneMap()
    {
        for(final GeneReadData gene : mGenes)
        {
            for(final TranscriptData transData : gene.getTranscripts())
            {
                mTransIdsGeneMap.put(transData.TransId, gene);
            }
        }
    }

    public boolean hasGeneData(int transId) { return mTransIdsGeneMap.containsKey(transId); }

    public GeneReadData findGeneData(int transId)
    {
        return mTransIdsGeneMap.get(transId);
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
