package com.hartwig.hmftools.common.serve.classification;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class FusionPairAndExonRangeMatcher implements EventMatcher {

    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXON_RANGES_PER_GENE = Maps.newHashMap();

    static {
        Set<String> kitSet = Sets.newHashSet("EXON 11 MUTATION", "Exon 11 mutations", "Exon 11 deletions");
        Set<String> metSet = Sets.newHashSet("EXON 14 SKIPPING MUTATION");

        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("KIT", kitSet);
        FUSION_PAIR_AND_EXON_RANGES_PER_GENE.put("MET", metSet);
    }

    public FusionPairAndExonRangeMatcher() {
    }

    @Override
    public boolean matches(@NotNull String gene, @NotNull String event) {
        Set<String> entries = FUSION_PAIR_AND_EXON_RANGES_PER_GENE.get(gene);
        if (entries != null) {
            return entries.contains(event);
        }

        return false;
    }
}
