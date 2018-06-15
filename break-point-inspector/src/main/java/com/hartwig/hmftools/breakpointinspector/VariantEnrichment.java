package com.hartwig.hmftools.breakpointinspector;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.breakpointinspector.datamodel.EnrichedVariantContext;
import com.hartwig.hmftools.breakpointinspector.datamodel.HMFStructuralVariantType;
import com.hartwig.hmftools.breakpointinspector.datamodel.ImmutableEnrichedVariantContext;
import com.hartwig.hmftools.breakpointinspector.datamodel.ImmutableRange;
import com.hartwig.hmftools.breakpointinspector.datamodel.Range;

import org.apache.commons.lang3.ObjectUtils;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;

final class VariantEnrichment {

    private VariantEnrichment() {
    }

    @NotNull
    static EnrichedVariantContext enrich(@NotNull VariantContext variant, @NotNull VariantContext mateVariant,
            @NotNull SAMSequenceDictionary sequenceDictionary) {
        final HMFStructuralVariantType structuralVariantType = toHMFSVType(variant);
        final boolean imprecise = variant.hasAttribute("IMPRECISE");
        final boolean translocation = variant.getStructuralVariantType() == StructuralVariantType.BND;

        final Location location1 = Location.parseFromVariant(variant, sequenceDictionary);
        final Range uncertainty1 = extractCIPOS(variant);

        final Location location2 = extractLocation2(variant, location1, sequenceDictionary);
        final Range uncertainty2 = extractUncertainty2(variant,
                mateVariant,
                uncertainty1,
                imprecise,
                translocation,
                structuralVariantType.isInversion());

        final Set<String> filters = variant.getFilters().stream().filter(s -> !s.startsWith("BPI")).collect(Collectors.toSet());

        final String homologySequence = variant.getAttributeAsString("HOMSEQ", "");
        final String insertSequence;
        if (variant.hasAttribute("LEFT_SVINSSEQ") && variant.hasAttribute("RIGHT_SVINSSEQ")) {
            insertSequence = variant.getAttributeAsString("LEFT_SVINSSEQ", "") + "..." + variant.getAttributeAsString("RIGHT_SVINSSEQ", "");
        } else {
            insertSequence = variant.getAttributeAsString("SVINSSEQ", "");
        }

        return ImmutableEnrichedVariantContext.builder()
                .variant(variant)
                .type(structuralVariantType)
                .location1(location1)
                .location2(location2)
                .imprecise(imprecise)
                .uncertainty1(uncertainty1)
                .uncertainty2(uncertainty2)
                .filters(filters)
                .insertSequence(insertSequence)
                .homologySequence(homologySequence)
                .isTranslocation(translocation)
                .orientationBP1(structuralVariantType.orientationBP1())
                .orientationBP2(structuralVariantType.orientationBP2())
                .build();
    }

    @NotNull
    private static HMFStructuralVariantType toHMFSVType(@NotNull VariantContext variant) {
        switch (variant.getStructuralVariantType()) {
            case INS:
                return HMFStructuralVariantType.INS;
            case INV:
                if (variant.hasAttribute("INV3")) {
                    return HMFStructuralVariantType.INV3;
                } else if (variant.hasAttribute("INV5")) {
                    return HMFStructuralVariantType.INV5;
                } else {
                    throw new IllegalStateException(variant.getID() + " : expected either INV3 or INV5 flag");
                }
            case DEL:
                return HMFStructuralVariantType.DEL;
            case DUP:
                return HMFStructuralVariantType.DUP;
            case BND:
                final String[] leftSplit = leftSplit(variant);
                final String[] rightSplit = rightSplit(variant);

                if (leftSplit.length >= 2) {
                    if (leftSplit[0].length() > 0) {
                        return HMFStructuralVariantType.INV3;
                    } else {
                        return HMFStructuralVariantType.DUP;
                    }
                } else if (rightSplit.length >= 2) {
                    if (rightSplit[0].length() > 0) {
                        return HMFStructuralVariantType.DEL;
                    } else {
                        return HMFStructuralVariantType.INV5;
                    }
                } else {
                    throw new IllegalStateException(variant.getID() + " : Could not interpret breakpoint");
                }
            default:
                throw new IllegalStateException(variant.getID() + " : Unexpected SV Type:" + variant.getStructuralVariantType());
        }
    }

    @NotNull
    private static Location extractLocation2(@NotNull VariantContext variant, @NotNull Location location1,
            @NotNull SAMSequenceDictionary sequenceDictionary) {
        switch (variant.getStructuralVariantType()) {
            case INS:
                return location1.set(variant.getAttributeAsInt("END", 0));
            case INV:
            case DEL:
            case DUP:
                return location1.add(Math.abs(variant.getAttributeAsInt("SVLEN", 0)));
            case BND:
                final String[] leftSplit = leftSplit(variant);
                final String[] rightSplit = rightSplit(variant);

                if (leftSplit.length >= 2) {
                    return Location.parseLocationString(leftSplit[1], sequenceDictionary);
                } else if (rightSplit.length >= 2) {
                    return Location.parseLocationString(rightSplit[1], sequenceDictionary);
                } else {
                    throw new IllegalStateException(variant.getID() + " : Could not parse breakpoint");
                }
            default:
                throw new IllegalStateException(variant.getID() + " : Unexpected SV Type:" + variant.getStructuralVariantType());
        }
    }

    @NotNull
    private static Range extractUncertainty2(@NotNull VariantContext variant, @NotNull VariantContext mateVariant,
            @NotNull Range uncertainty1, boolean imprecise, boolean translocation, boolean inversion) {
        final List<Integer> ciEndList = variant.getAttributeAsIntList("CIEND", 0);
        Range uncertainty2 = ciEndList.size() == 2 ? ImmutableRange.of(ciEndList.get(0), ciEndList.get(1)) : null;

        if (translocation) {
            final String[] leftSplit = leftSplit(variant);
            final String[] rightSplit = rightSplit(variant);

            if (leftSplit.length >= 2) {
                if (leftSplit[0].length() > 0) {
                    uncertainty2 = Range.invert(uncertainty1);
                } else {
                    uncertainty2 = uncertainty1;
                }
            } else if (rightSplit.length >= 2) {
                if (rightSplit[0].length() > 0) {
                    uncertainty2 = uncertainty1;
                } else {
                    uncertainty2 = Range.invert(uncertainty1);
                }
            } else {
                throw new IllegalStateException(variant.getID() + " : Could not parse breakpoint");
            }

            if (imprecise) {
                uncertainty2 = extractCIPOS(mateVariant);
            }
        }

        return ObjectUtils.firstNonNull(uncertainty2, fixup(uncertainty1, imprecise, inversion));
    }

    @NotNull
    private static String[] leftSplit(@NotNull VariantContext variant) {
        return split(variant, "\\]");
    }

    @NotNull
    private static String[] rightSplit(@NotNull VariantContext variant) {
        return split(variant, "\\[");
    }

    @NotNull
    private static String[] split(@NotNull VariantContext variant, @NotNull String splitRegExp) {
        return variant.getAlternateAllele(0).getDisplayString().split(splitRegExp);
    }

    @NotNull
    private static Range extractCIPOS(@NotNull VariantContext variant) {
        final List<Integer> ciPosList = variant.getAttributeAsIntList("CIPOS", 0);
        return ciPosList.size() == 2 ? ImmutableRange.of(ciPosList.get(0), ciPosList.get(1)) : ImmutableRange.of(0, 0);
    }

    @NotNull
    private static Range fixup(@NotNull Range uncertainty1, boolean imprecise, boolean inversion) {
        if (imprecise) {
            return uncertainty1;
        } else {
            return inversion ? Range.invert(uncertainty1) : uncertainty1;
        }
    }
}
