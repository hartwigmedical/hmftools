package com.hartwig.hmftools.common.genome.genepanel;

import java.io.InputStream;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

public final class HmfGenePanelSupplier {

    private HmfGenePanelSupplier() {
    }

    @NotNull
    public static List<HmfTranscriptRegion> allGeneList37() {
        InputStream inputStream = HmfGenePanelSupplier.class.getResourceAsStream("/ensembl/all_genes.37.tsv");
        return HmfGenomeFileLoader.fromInputStream(inputStream);
    }

    @NotNull
    public static List<HmfTranscriptRegion> allGeneList38() {
        InputStream inputStream = HmfGenePanelSupplier.class.getResourceAsStream("/ensembl/all_genes.38.tsv");
        return HmfGenomeFileLoader.fromInputStream(inputStream);
    }

    @NotNull
    public static Map<String, HmfTranscriptRegion> allGenesMap37() {
        return allGenesMap(allGeneList37());
    }

    @NotNull
    public static Map<String, HmfTranscriptRegion> allGenesMap38() {
        return allGenesMap(allGeneList38());
    }

    @NotNull
    public static SortedSetMultimap<String, HmfTranscriptRegion> allGenesPerChromosomeMap37() {
        return toSortedMap(allGeneList37());
    }

    @NotNull
    private static Map<String, HmfTranscriptRegion> allGenesMap(@NotNull List<HmfTranscriptRegion> regions) {
        Map<String, HmfTranscriptRegion> geneMap = Maps.newHashMap();
        for (HmfTranscriptRegion region : regions) {
            assert !geneMap.containsKey(region.gene());
            geneMap.put(region.gene(), region);
        }
        return geneMap;
    }

    @NotNull
    private static SortedSetMultimap<String, HmfTranscriptRegion> toSortedMap(@NotNull List<HmfTranscriptRegion> regions) {
        SortedSetMultimap<String, HmfTranscriptRegion> regionMap = TreeMultimap.create();
        for (HmfTranscriptRegion region : regions) {
            regionMap.put(region.chromosome(), region);
        }

        return regionMap;
    }
}