package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ValidGermlineVariantPredicateTest {

    private static final GermlineSampleData VARIANT = new GermlineSampleData("1/1", 10, 5);
    private static final GermlineSampleData NOT_VARIANT = new GermlineSampleData("./.", 10, 5);

    @Test
    public void checkSNPVariant() {
        final GermlineVariant refSNPVariant = createVariant(VariantType.SNP, VARIANT, NOT_VARIANT);
        final GermlineVariant noVariant = createVariant(VariantType.SNP, NOT_VARIANT, NOT_VARIANT);

        final Predicate<GermlineVariant> snpPredicate = new ValidGermlineVariantPredicate(VariantType.SNP, true);
        assertTrue(snpPredicate.test(refSNPVariant));
        assertFalse(snpPredicate.test(noVariant));

        final Predicate<GermlineVariant> indelPredicate = new ValidGermlineVariantPredicate(VariantType.INDEL, true);
        assertFalse(indelPredicate.test(refSNPVariant));
        assertFalse(indelPredicate.test(noVariant));
    }

    @Test
    public void checkIndelVariant() {
        final GermlineVariant refIndelVariant = createVariant(VariantType.INDEL, VARIANT, NOT_VARIANT);
        final GermlineVariant noVariant = createVariant(VariantType.INDEL, NOT_VARIANT, NOT_VARIANT);

        final Predicate<GermlineVariant> indelPredicate = new ValidGermlineVariantPredicate(VariantType.INDEL, true);
        assertTrue(indelPredicate.test(refIndelVariant));
        assertFalse(indelPredicate.test(noVariant));

        final Predicate<GermlineVariant> snpPredicate = new ValidGermlineVariantPredicate(VariantType.SNP, true);
        assertFalse(snpPredicate.test(refIndelVariant));
        assertFalse(snpPredicate.test(noVariant));
    }

    @NotNull
    private static GermlineVariant createVariant(@NotNull final VariantType type,
            @NotNull final GermlineSampleData refData, @Nullable final GermlineSampleData tumorData) {
        return new GermlineVariant(type, "AnyFilter", refData, tumorData);
    }
}
