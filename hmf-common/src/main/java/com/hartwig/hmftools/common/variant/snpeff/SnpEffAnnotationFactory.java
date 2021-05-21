package com.hartwig.hmftools.common.variant.snpeff;

import static com.hartwig.hmftools.common.variant.VariantConsequence.sufferConsequences;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_SEQUENCE_FLAG;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public final class SnpEffAnnotationFactory {

    private static final Logger LOGGER = LogManager.getLogger(SnpEffAnnotationFactory.class);

    public static final String SNPEFF_IDENTIFIER = "ANN";

    static final String SPLICE_DONOR_VARIANT = "splice_donor_variant";

    private static final String FIELD_SEPARATOR = "\\|";
    static final String CONSEQUENCE_SEPARATOR = "&";

    private static final int EXPECTED_FIELD_SIZE_PER_ANNOTATION = 16;

    private SnpEffAnnotationFactory() {
    }

    @NotNull
    public static List<SnpEffAnnotation> fromContext(@NotNull final VariantContext context) {
        if (context.hasAttribute(SNPEFF_IDENTIFIER)) {
            return fromAnnotationList(context, context.getAttributeAsStringList(SNPEFF_IDENTIFIER, ""));
        }
        return Collections.emptyList();
    }

    @NotNull
    private static List<SnpEffAnnotation> fromAnnotationList(@NotNull final VariantContext context,
            @NotNull final List<String> annotation) {
        return annotation.stream()
                .map(x -> enforceMinLength(x.trim().split(FIELD_SEPARATOR), EXPECTED_FIELD_SIZE_PER_ANNOTATION))
                .filter(SnpEffAnnotationFactory::isCorrectNumberOfParts)
                .map(x -> fromParts(context, x))
                .collect(Collectors.toList());
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

    @NotNull
    private static SnpEffAnnotation fromParts(@NotNull VariantContext context, @NotNull String[] parts) {
        String effects = effect(context, parts);

        return ImmutableSnpEffAnnotation.builder()
                .allele(parts[0])
                .effects(effects)
                .consequences(sufferConsequences(toEffects(effects)))
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
    static String effect(@NotNull VariantContext variant, @NotNull String[] parts) {
        String hgvsCoding = parts[9];
        String effects = parts[1];
        if (!effects.contains("splice")) {
            return effects;
        }

        // Below is to support additional donor annotation for variants affecting +5 base. See also DEV-1650
        // Note, no longer strictly necessary to do here as have also got this in CodingEffectFactory.java
        boolean indel = variant.isIndel();
        if (!indel) {
            return hgvsCoding.contains("+5") ? SPLICE_DONOR_VARIANT + CONSEQUENCE_SEPARATOR + effects : effects;
        }

        final String hgvsCodingType;
        if (hgvsCoding.contains("ins")) {
            hgvsCodingType = "ins";
        } else if (hgvsCoding.contains("del")) {
            hgvsCodingType = "del";
        } else if (hgvsCoding.contains("dup")) {
            hgvsCodingType = "dup";
        } else {
            return effects;
        }

        int initialSpliceBase = initialIndelSpliceBase(hgvsCodingType.equals("ins"), hgvsCoding);
        if (initialSpliceBase == -1) {
            return effects;
        }

        final int adjustedSpliceBase;
        final String ref = variant.getReference().getBaseString();
        final String alt = variant.getAlternateAllele(0).getBaseString();
        if (isPositiveStrand(ref, alt, hgvsCoding)) {
            final String variantBases = ref.length() > alt.length() ? ref.substring(1) : alt.substring(1);
            final int microhomologyAdditionalBases = variant.getAttributeAsString(MICROHOMOLOGY_FLAG, Strings.EMPTY).length();
            final String repeatSequence = variant.getAttributeAsString(REPEAT_SEQUENCE_FLAG, Strings.EMPTY);
            final int repeatCount = variant.getAttributeAsInt(REPEAT_COUNT_FLAG, 0);
            final int repeatCountAdditionalBases = variantBases.equals(repeatSequence) ? repeatCount * repeatSequence.length() : 0;
            adjustedSpliceBase = initialSpliceBase + Math.max(microhomologyAdditionalBases, repeatCountAdditionalBases);
        } else {
            adjustedSpliceBase = initialSpliceBase;
        }

        return adjustedSpliceBase <= 5 ? SPLICE_DONOR_VARIANT + CONSEQUENCE_SEPARATOR + effects : effects;
    }

    @VisibleForTesting
    static int initialIndelSpliceBase(final boolean isInsert, final String hgvsCoding) {
        int firstIndexOfPlus = hgvsCoding.indexOf("+");
        if (firstIndexOfPlus < 0) {
            return -1;
        }

        try {
            int spliceLocation = Integer.parseInt(hgvsCoding.substring(firstIndexOfPlus + 1, firstIndexOfPlus + 2));
            int result = isInsert ? spliceLocation + 1 : spliceLocation;
            return result <= 5 ? result : -1;

        } catch (NumberFormatException e) {
            return -1;
        }
    }

    private static boolean isPositiveStrand(@NotNull final String ref, final String alt, final String hgvsCoding) {
        char lastBaseOfCoding = hgvsCoding.charAt(hgvsCoding.length() - 1);
        return ref.length() > alt.length()
                ? ref.charAt(ref.length() - 1) == lastBaseOfCoding
                : alt.charAt(alt.length() - 1) == lastBaseOfCoding;
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
}
