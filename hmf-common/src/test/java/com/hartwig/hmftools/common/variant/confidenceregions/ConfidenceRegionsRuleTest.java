package com.hartwig.hmftools.common.variant.confidenceregions;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConfidenceRegionsRuleTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void confidenceRegionsRuleWorks() {
        final Slicer highConfidence = SlicerFactory.fromSingleGenomeRegion(region(100, 1000));
        final Slicer extremeConfidence = SlicerFactory.fromSingleGenomeRegion(region(500, 600));

        final ConfidenceRegionsRule rule = ConfidenceRegionsRule.fromSlicers(highConfidence, extremeConfidence);

        final List<SomaticVariant> variants = Lists.newArrayList(
                // @formatter:off
                cosmicSNPVariantOnPositionWithCallers(300, 2),// MIVO: high-conf, cosmic => include
                dbsnpSNPVariantOnPositionWithCallers(400, 2), // MIVO: high-conf, dbsnp => exclude
                dbsnpSNPVariantOnPositionWithCallers(550, 2), // MIVO: extreme-conf => include
                dbsnpSNPVariantOnPositionWithCallers(2000, 4), // MIVO: low-conf => exclude
                indelVariantOnPositionWithCallers(550, 1), // MIVO: extreme-conf => include
                indelVariantOnPositionWithCallers(200, 1), // MIVO: high-conf => exclude
                indelVariantOnPositionWithCallers(250, 2) // MIVO: high-conf => exclude
                // @formatter:on
        );

        final List<SomaticVariant> filtered = rule.removeUnreliableVariants(variants);
        assertEquals(3, filtered.size());
    }

    @Test
    public void updateFilterFlagWorks() {
        final String filtered = "FILTERED";
        final String pass = "PASS";

        final SomaticVariant variant1 = new SomaticVariant.Builder().filter(pass).chromosome("1").build();
        final SomaticVariant variant2 = new SomaticVariant.Builder().filter(filtered).chromosome("1").build();
        final SomaticVariant variant3 = new SomaticVariant.Builder().filter(pass).chromosome("2").build();

        final Predicate<SomaticVariant> filter = variant -> variant.chromosome().equals("2");

        final ConfidenceRegionsRule confidenceRegionsRule = new ConfidenceRegionsRule(filter);

        final List<SomaticVariant> adjustedVariants =
                confidenceRegionsRule.updateFilterFlagForUnreliableVariants(Lists.newArrayList(variant1, variant2, variant3));

        assertEquals(ConfidenceRegionsRule.FILTER, adjustedVariants.get(0).filter());
        assertEquals(filtered + ";" + ConfidenceRegionsRule.FILTER, adjustedVariants.get(1).filter());
        assertEquals(pass, adjustedVariants.get(2).filter());
    }

    @NotNull
    private static GenomeRegion region(final long start, final long end) {
        return ImmutableBEDGenomeRegion.of(CHROMOSOME, start, end);
    }

    @NotNull
    private static SomaticVariant cosmicSNPVariantOnPositionWithCallers(long position, int numCallers) {
        final List<String> callers = Lists.newArrayList();
        for (int i = 0; i < numCallers; i++) {
            callers.add("any");
        }
        return new SomaticVariant.Builder().type(VariantType.SNP)
                .chromosome(CHROMOSOME)
                .position(position)
                .callers(callers)
                .cosmicID("any_id")
                .build();
    }

    @NotNull
    private static SomaticVariant dbsnpSNPVariantOnPositionWithCallers(long position, int numCallers) {
        final List<String> callers = Lists.newArrayList();
        for (int i = 0; i < numCallers; i++) {
            callers.add("any");
        }
        return new SomaticVariant.Builder().type(VariantType.SNP)
                .chromosome(CHROMOSOME)
                .position(position)
                .callers(callers)
                .dnsnpID("any_id")
                .build();
    }

    @NotNull
    private static SomaticVariant indelVariantOnPositionWithCallers(long position, int numCallers) {
        final List<String> callers = Lists.newArrayList();
        for (int i = 0; i < numCallers; i++) {
            callers.add("any");
        }
        return new SomaticVariant.Builder().type(VariantType.INDEL).chromosome(CHROMOSOME).position(position).callers(callers).build();
    }
}