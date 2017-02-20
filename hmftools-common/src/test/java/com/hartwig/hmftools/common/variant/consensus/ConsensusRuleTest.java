package com.hartwig.hmftools.common.variant.consensus;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.slicing.GenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConsensusRuleTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void consensusRuleWorks() {
        final Slicer highConfidence = SlicerFactory.forGenomeRegion(region(100, 1000));
        final Slicer cpctSlicing = SlicerFactory.forGenomeRegion(region(500, 600));

        final ConsensusRule rule = new ConsensusRule(highConfidence, cpctSlicing);

        final List<SomaticVariant> variants = Lists.newArrayList(cosmicSNPVariantOnPositionWithCallers(300, 2),
                // KODU: Include
                dbsnpSNPVariantOnPositionWithCallers(400, 2), // KODU: Exclude
                dbsnpSNPVariantOnPositionWithCallers(550, 2), // KODU: Include
                dbsnpSNPVariantOnPositionWithCallers(2000, 4), // KODU: Include
                indelVariantOnPositionWithCallers(550, 1), // KODU: Include
                indelVariantOnPositionWithCallers(200, 1), // KODU: Exclude
                indelVariantOnPositionWithCallers(250, 2) // KODU: Include
        );

        final List<SomaticVariant> filtered = rule.apply(variants);
        assertEquals(5, filtered.size());
    }

    @NotNull
    private static GenomeRegion region(final long start, final long end) {
        return new GenomeRegion(CHROMOSOME, start, end);
    }

    @NotNull
    private static SomaticVariant cosmicSNPVariantOnPositionWithCallers(long position, int numCallers) {
        final List<String> callers = Lists.newArrayList();
        for (int i = 0; i < numCallers; i++) {
            callers.add("any");
        }
        return new SomaticVariant.Builder().type(VariantType.SNP).chromosome(CHROMOSOME).position(position).callers(
                callers).cosmicID("any_id").build();
    }

    @NotNull
    private static SomaticVariant dbsnpSNPVariantOnPositionWithCallers(long position, int numCallers) {
        final List<String> callers = Lists.newArrayList();
        for (int i = 0; i < numCallers; i++) {
            callers.add("any");
        }
        return new SomaticVariant.Builder().type(VariantType.SNP).chromosome(CHROMOSOME).position(position).callers(
                callers).dnsnpID("any_id").build();
    }

    @NotNull
    private static SomaticVariant indelVariantOnPositionWithCallers(long position, int numCallers) {
        final List<String> callers = Lists.newArrayList();
        for (int i = 0; i < numCallers; i++) {
            callers.add("any");
        }
        return new SomaticVariant.Builder().type(VariantType.INDEL).chromosome(CHROMOSOME).position(position).callers(
                callers).build();
    }
}