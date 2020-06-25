package com.hartwig.hmftools.serve.vicc.range;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

public class GeneRangeExtractor {

    private static final Set<String> GENE_EXON = Sets.newHashSet("exon");
    private static final Set<String> GENE_MULTIPLE_CODONS = Sets.newHashSet("V600X ",
            "H1047X ",
            "V600E/K ",
            "Q61X ",
            "G12X ",
            "H1047X ",
            "G13X ",
            "E709X ",
            "D835X ",
            "G719X ",
            "D842X ",
            "R132X ");

    @NotNull
    public Map<Feature, String> extractGeneRanges(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> geneRangesPerFeature = Maps.newHashMap();
        for (Feature feature : viccEntry.features()) {
            String event = feature.name();
            if (viccEntry.source() == ViccSource.JAX) {
                String[] extractEvent = event.split(" ", 2);
                if (extractEvent.length == 2) {
                    event = extractEvent[1];
                }
            }

            if (GENE_EXON.contains(event.toLowerCase())) {
                geneRangesPerFeature.put(feature, feature.name());
            } else if (GENE_MULTIPLE_CODONS.contains(event)) {
                //TODO: possible using transvar for the codons
                geneRangesPerFeature.put(feature, feature.name());
            }
        }

        return geneRangesPerFeature;
    }
}
