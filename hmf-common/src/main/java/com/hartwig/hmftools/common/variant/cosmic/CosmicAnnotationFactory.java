package com.hartwig.hmftools.common.variant.cosmic;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class CosmicAnnotationFactory {

    private static final String COSMIC_IDENTIFIER = "COSM2ENST";

    private static final String FIELD_SEPARATOR = "\\|";
    private static final String GENE_TRANSCRIPT_SEPARATOR = "_";

    private CosmicAnnotationFactory() {
    }

    @NotNull
    public static List<CosmicAnnotation> fromContext(@NotNull final VariantContext context) {
        if (context.hasAttribute(COSMIC_IDENTIFIER)) {
            return context.getAttributeAsStringList(COSMIC_IDENTIFIER, "")
                    .stream()
                    .map(x -> x.trim().split(FIELD_SEPARATOR))
                    .map(CosmicAnnotationFactory::fromParts)
                    .collect(Collectors.toList());

        }
        return Collections.emptyList();
    }

    @NotNull
    private static CosmicAnnotation fromParts(@NotNull final String[] parts) {
        final String[] geneTranscriptFields = parts[1].split(GENE_TRANSCRIPT_SEPARATOR);

        return ImmutableCosmicAnnotation.builder()
                .id(parts[0])
                .gene(geneTranscriptFields[0])
                .transcript(geneTranscriptFields[1])
                .hgvsCoding(parts[2])
                .hgvsProtein(parts[3])
                .count(Integer.parseInt(parts[4]))
                .build();
    }
}
