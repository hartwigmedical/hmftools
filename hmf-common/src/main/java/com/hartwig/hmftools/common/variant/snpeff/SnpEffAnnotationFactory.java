package com.hartwig.hmftools.common.variant.snpeff;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class SnpEffAnnotationFactory {

    private static final Logger LOGGER = LogManager.getLogger(SnpEffAnnotationFactory.class);

    public static final String SNPEFF_IDENTIFIER = "ANN";

    private static final String FIELD_SEPARATOR = "\\|";
    private static final String CONSEQUENCE_SEPARATOR = "&";

    private static final int EXPECTED_FIELD_SIZE_PER_ANNOTATION = 16;

    private SnpEffAnnotationFactory() {
    }

    @NotNull
    public static List<SnpEffAnnotation> fromContext(@NotNull final VariantContext context) {
        if (context.hasAttribute(SNPEFF_IDENTIFIER)) {
            return context.getAttributeAsStringList(SNPEFF_IDENTIFIER, "")
                    .stream()
                    .map(x -> enforceMinLength(x.trim().split(FIELD_SEPARATOR), EXPECTED_FIELD_SIZE_PER_ANNOTATION))
                    .filter(SnpEffAnnotationFactory::isCorrectNumberOfParts)
                    .map(SnpEffAnnotationFactory::fromParts)
                    .collect(Collectors.toList());

        }
        return Collections.emptyList();
    }

    @NotNull
    public static List<String> rawAnnotations(@NotNull final VariantContext context) {
        if (context.hasAttribute(SNPEFF_IDENTIFIER)) {
            return context.getAttributeAsStringList(SNPEFF_IDENTIFIER, "");
        }

        return Collections.emptyList();
    }

    private static boolean isCorrectNumberOfParts(@NotNull String[] parts) {
        if (parts.length == EXPECTED_FIELD_SIZE_PER_ANNOTATION) {
            return true;
        }

        final StringJoiner joiner = new StringJoiner("|");
        Stream.of(parts).forEach(joiner::add);

        LOGGER.warn("Annotation found with invalid field count: " + joiner.toString());
        return false;
    }

    @VisibleForTesting
    @NotNull
    static SnpEffAnnotation fromParts(@NotNull final String[] parts) {
        return ImmutableSnpEffAnnotation.builder()
                .allele(parts[0])
                .effects(parts[1])
                .consequences(toConsequences(toEffects(parts[1])))
                .severity(parts[2])
                .gene(parts[3])
                .geneID(parts[4])
                .featureType(parts[5])
                .featureID(parts[6])
                .transcriptBioType(parts[7])
                .rank(parts[8])
                .hgvsCoding(parts[9])
                .hgvsProtein(parts[10])
                .cDNAPosAndLength(parts[11])
                .cdsPosAndLength(parts[12])
                .aaPosAndLength(parts[13])
                .distance(parts[14])
                .addition(parts[15])
                .build();
    }

    @NotNull
    private static String[] enforceMinLength(@NotNull String[] parts, int minSize) {
        if (parts.length > minSize) {
            return parts;
        } else {
            final String[] values = new String[minSize];
            for (int i = 0; i < minSize; i++) {
                values[i] = i < parts.length ? parts[i] : Strings.EMPTY;
            }
            System.arraycopy(parts, 0, values, 0, parts.length);

            return values;
        }
    }

    @NotNull
    private static List<String> toEffects(@NotNull final String effectString) {
        return Lists.newArrayList(effectString.split(CONSEQUENCE_SEPARATOR));
    }

    @NotNull
    private static List<VariantConsequence> toConsequences(@NotNull final List<String> effects) {
        final List<VariantConsequence> consequences = Lists.newArrayList();
        for (final String part : effects) {
            boolean found = false;
            for (final VariantConsequence consequence : VariantConsequence.values()) {
                if (consequence.isParentTypeOf(part)) {
                    found = true;
                    consequences.add(consequence);
                }
            }
            if (!found) {
                consequences.add(VariantConsequence.OTHER);
            }
        }
        return consequences;
    }
}
