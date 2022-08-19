package com.hartwig.hmftools.serve.extraction.snpeff;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class SnpEffAnnotationTestFactory
{

    private SnpEffAnnotationTestFactory()
    {
    }

    @NotNull
    static ImmutableSnpEffAnnotation.Builder builder()
    {
        return ImmutableSnpEffAnnotation.builder()
                .gene(Strings.EMPTY)
                .geneID(Strings.EMPTY)
                .featureType(Strings.EMPTY)
                .featureID(Strings.EMPTY)
                .rank(Strings.EMPTY)
                .hgvsCoding(Strings.EMPTY)
                .hgvsProtein(Strings.EMPTY)
                .effects(Strings.EMPTY);
    }
}
