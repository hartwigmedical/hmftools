package com.hartwig.hmftools.serve.vicc.range;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

public class GeneRangeExtractor {

    @NotNull
    public Map<Feature, String> extractGeneRanges(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> geneRangesPerFeature = Maps.newHashMap();
        if (viccEntry.source() == ViccSource.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (feature.name().toLowerCase().contains("exon")) {
                    geneRangesPerFeature.put(feature, feature.name());
                }
            }
        }
        return geneRangesPerFeature;
    }
}
