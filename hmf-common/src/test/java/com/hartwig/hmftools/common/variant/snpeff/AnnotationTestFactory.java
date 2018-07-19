package com.hartwig.hmftools.common.variant.snpeff;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class AnnotationTestFactory {

    private AnnotationTestFactory() {
    }

    @NotNull
    public static ImmutableSnpEffAnnotation.Builder createVariantAnnotationBuilder(@NotNull VariantConsequence... consequences) {
        return ImmutableSnpEffAnnotation.builder()
                .allele(Strings.EMPTY)
                .severity(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .geneID(Strings.EMPTY)
                .featureType(Strings.EMPTY)
                .featureID(Strings.EMPTY)
                .transcriptBioType(Strings.EMPTY)
                .rank(Strings.EMPTY)
                .hgvsCoding(Strings.EMPTY)
                .hgvsProtein(Strings.EMPTY)
                .cDNAPosAndLength(Strings.EMPTY)
                .cdsPosAndLength(Strings.EMPTY)
                .aaPosAndLength(Strings.EMPTY)
                .distance(Strings.EMPTY)
                .addition(Strings.EMPTY)
                .effects(Strings.EMPTY)
                .consequences(Lists.newArrayList(consequences));
    }
}
