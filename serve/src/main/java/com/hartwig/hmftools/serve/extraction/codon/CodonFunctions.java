package com.hartwig.hmftools.serve.extraction.codon;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class CodonFunctions {

    private CodonFunctions() {
    }

    @NotNull
    public static Set<KnownCodon> consolidate(@NotNull Iterable<KnownCodon> codons) {
        Map<CodonAnnotation, Set<Knowledgebase>> sourcesPerAnnotation = Maps.newHashMap();
        for (KnownCodon codon : codons) {
            Set<Knowledgebase> sources = sourcesPerAnnotation.get(codon.annotation());
            if (sources == null) {
                sources = Sets.newHashSet();
            }
            sources.addAll(codon.sources());
            sourcesPerAnnotation.put(codon.annotation(), sources);
        }

        Set<KnownCodon> consolidated = Sets.newHashSet();
        for (Map.Entry<CodonAnnotation, Set<Knowledgebase>> entry : sourcesPerAnnotation.entrySet()) {
            consolidated.add(ImmutableKnownCodon.builder()
                    .annotation(entry.getKey())
                    .sources(entry.getValue())
                    .build());
        }
        return consolidated;
    }
}
