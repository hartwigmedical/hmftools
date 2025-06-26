package com.hartwig.hmftools.common.driver.panel;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.GeneRegion;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.BaseRegion;

public final class DriverGeneRegions
{
    public static List<String> loadActionableGenes(final List<DriverGene> driverGenes)
    {
        return driverGenes.stream()
                .filter(x -> x.reportSomatic() || x.reportGermline() || x.reportPGX())
                .map(x -> x.gene()).collect(Collectors.toList());
    }

    public static Map<String,List<GeneRegion>> buildDriverGeneRegions(
            final EnsemblDataCache ensemblDataCache, final List<String> driverGenes)
    {
        Map<String,List<GeneRegion>> chrGeneRegionsMap = Maps.newHashMap();

        Map<String,List<GeneData>> chrGeneDataMap = ensemblDataCache.getChrGeneDataMap();

        Set<String> geneSet = driverGenes.stream().collect(Collectors.toSet());

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = ensemblDataCache.refGenomeVersion().versionedChromosome(chromosome.toString());
            List<GeneData> geneDataList = chrGeneDataMap.get(chromosomeStr);

            List<GeneRegion> chrPanelRegions = Lists.newArrayList();

            chrGeneRegionsMap.put(chromosomeStr, chrPanelRegions);

            for(GeneData geneData : geneDataList)
            {
                if(!geneSet.contains(geneData.GeneName))
                    continue;

                TranscriptData transData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                List<GeneRegion> transcriptRegions = getTranscriptRegions(geneData, transData, false, true);

                chrPanelRegions.addAll(transcriptRegions);
            }

            // sort and merge any overlaps
            Collections.sort(chrPanelRegions);
        }

        return chrGeneRegionsMap;
    }

    public static Map<String,List<BaseRegion>> buildDriverGeneBaseRegions(
            final EnsemblDataCache ensemblDataCache, final List<String> driverGenes, boolean includeUTR, boolean includeSplice)
    {
        Map<String,List<BaseRegion>> chrRegionsMap = Maps.newHashMap();

        Map<String,List<GeneData>> chrGeneDataMap = ensemblDataCache.getChrGeneDataMap();

        Set<String> geneSet = driverGenes.stream().collect(Collectors.toSet());

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = ensemblDataCache.refGenomeVersion().versionedChromosome(chromosome.toString());
            List<GeneData> geneDataList = chrGeneDataMap.get(chromosomeStr);

            List<GeneRegion> chrPanelRegions = Lists.newArrayList();

            for(GeneData geneData : geneDataList)
            {
                if(!geneSet.contains(geneData.GeneName))
                    continue;

                TranscriptData transData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                List<GeneRegion> transcriptRegions = getTranscriptRegions(geneData, transData, includeUTR, includeSplice);

                chrPanelRegions.addAll(transcriptRegions);
            }

            // sort and merge any overlaps
            Collections.sort(chrPanelRegions);

            int regionsRemoved = 0;

            int index = 0;

            // check for overlaps with the previous region
            while(index < chrPanelRegions.size() - 1)
            {
                GeneRegion region = chrPanelRegions.get(index);

                int nextIndex = index + 1;
                while(nextIndex < chrPanelRegions.size())
                {
                    GeneRegion nextRegion = chrPanelRegions.get(nextIndex);

                    if(region.end() >= nextRegion.start())
                    {
                        // LOGGER.trace("gene({}) merged region({}) with next({})", region.GeneName, region, nextRegion);

                        if(nextRegion.end() > region.end())
                        {
                            region.setEnd(nextRegion.end());
                        }
                        ++regionsRemoved;
                        chrPanelRegions.remove(nextIndex);
                    }
                    else
                    {
                        break;
                    }
                }

                ++index;
            }

            if(regionsRemoved > 0)
            {
                // LOGGER.debug("chr({}) merged {} regions from overlaps", chromosomeStr, regionsRemoved);
            }

            List<BaseRegion> chrRegions = chrPanelRegions.stream().map(x -> new BaseRegion(x.start(), x.end())).collect(Collectors.toList());
            chrRegionsMap.put(chromosomeStr, chrRegions);

        }

        return chrRegionsMap;
    }

    private static final int SPLICE_SIZE = 10;

    public static List<GeneRegion> getTranscriptRegions(
            final GeneData geneData, final TranscriptData transData, boolean includeUTR, boolean includeSplice)
    {
        int startPosition = includeUTR || transData.nonCoding() ? transData.TransStart : transData.CodingStart;
        int endPosition = includeUTR || transData.nonCoding() ? transData.TransEnd : transData.CodingEnd;

        final List<GeneRegion> regions = Lists.newArrayList();

        for(int i = 0; i < transData.exons().size(); i++)
        {
            ExonData exon = transData.exons().get(i);

            int exonStart = i == 0 ? exon.Start : exon.Start - SPLICE_SIZE;
            int exonEnd = i == transData.exons().size() - 1 ? exon.End : exon.End + SPLICE_SIZE;

            if(includeSplice)
            {
                exonStart -= SPLICE_SIZE;
                exonEnd += SPLICE_SIZE;
            }

            if(positionsOverlap(startPosition, endPosition, exonStart, exonEnd))
            {
                regions.add(new GeneRegion(
                        geneData.Chromosome, max(startPosition, exonStart), min(endPosition, exonEnd), geneData.GeneName, exon.Rank));
            }
        }

        return regions;
    }
}
