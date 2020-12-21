package com.hartwig.hmftools.serve.extraction.exon;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public final class ExonFunctions {

    private ExonFunctions() {
    }

    @NotNull
    public static Set<KnownExon> consolidate(@NotNull Iterable<KnownExon> exons) {
        Map<ExonAnnotation, Set<Knowledgebase>> sourcesPerAnnotation = Maps.newHashMap();
        for (KnownExon exon : exons) {
            Set<Knowledgebase> sources = sourcesPerAnnotation.get(exon.annotation());
            if (sources == null) {
                sources = Sets.newHashSet();
            }
            sources.addAll(exon.sources());
            sourcesPerAnnotation.put(exon.annotation(), sources);
        }

        Set<KnownExon> consolidated = Sets.newHashSet();
        for (Map.Entry<ExonAnnotation, Set<Knowledgebase>> entry : sourcesPerAnnotation.entrySet()) {
            consolidated.add(ImmutableKnownExon.builder()
                    .annotation(entry.getKey())
                    .sources(entry.getValue())
                    .build());
        }
        return consolidated;
    }
}
