package com.hartwig.hmftools.compar.mutation;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.ArrayList;
import java.util.Set;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticVariantComparerTest
{
    @Test
    public void emptyComparison()
    {
        var config = new ComparConfig();
        var victim = new SomaticVariantComparer(config);

        var diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        var sampleId = "TEST";
        var mismatches = new ArrayList<Mismatch>();
        var refVariants = new ArrayList<SomaticVariantData>();
        var newVariants = new ArrayList<SomaticVariantData>();
        var matchLevel = MatchLevel.DETAILED;

        var result = victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel);
        assertTrue(result);
        assertTrue(mismatches.isEmpty());
    }

    @Test
    public void simpleDetailedComparison()
    {
        var config = new ComparConfig();
        var victim = new SomaticVariantComparer(config);

        var diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        var sampleId = "TEST";
        var mismatches = new ArrayList<Mismatch>();
        var refVariants = new ArrayList<SomaticVariantData>();
        var newVariants = new ArrayList<SomaticVariantData>();
        var matchLevel = MatchLevel.DETAILED;

        final SomaticVariantData baseVariant = createSomaticVariantData();

        // match
        refVariants.add(baseVariant);
        newVariants.add(baseVariant);

        // value difference
        refVariants.add(baseVariant.withAlt("C").withReported(false));
        newVariants.add(baseVariant.withAlt("C").withCanonicalHgvsProteinImpact("p.Val600Gly").withReported(false));

        // ref only
        refVariants.add(baseVariant.withChromosome("8").withReported(false));
        newVariants.add(baseVariant.withChromosome("8").withReported(false).withIsFromUnfilteredVcf(true));

        //new only
        newVariants.add(baseVariant.withType(VariantType.INDEL));

        var result = victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel);
        assertTrue(result);
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }

    @Test
    public void simpleReportableComparison()
    {
        var config = new ComparConfig();
        var victim = new SomaticVariantComparer(config);

        var diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        var sampleId = "TEST";
        var mismatches = new ArrayList<Mismatch>();
        var refVariants = new ArrayList<SomaticVariantData>();
        var newVariants = new ArrayList<SomaticVariantData>();
        var matchLevel = MatchLevel.REPORTABLE;

        final SomaticVariantData baseVariant = createSomaticVariantData();

        // match
        refVariants.add(baseVariant);
        newVariants.add(baseVariant);

        // value difference
        refVariants.add(baseVariant.withAlt("C"));
        newVariants.add(baseVariant.withAlt("C").withCanonicalHgvsProteinImpact("p.Val600Gly"));

        // ref only
        refVariants.add(baseVariant.withChromosome("8"));
        newVariants.add(baseVariant.withChromosome("8").withReported(false));

        //new only
        refVariants.add(baseVariant.withType(VariantType.INDEL).withIsFromUnfilteredVcf(true).withReported(false));
        newVariants.add(baseVariant.withType(VariantType.INDEL));

        var result = victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel);
        assertTrue(result);
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }

    @NotNull
    private static SomaticVariantData createSomaticVariantData()
    {
        return new SomaticVariantData("7", 140453136, "A", "T", VariantType.SNP, "BRAF", true,
                Hotspot.HOTSPOT, VariantTier.HOTSPOT, false, "missense_variant", "MISSENSE",
                "c.1799T>A", "p.Val600Glu", null, false,
                275, 0., Set.of("PASS"), 1.1, 0.45,
                new AllelicDepth(116, 21), false, true,
                "7", 140453136);
    }
}
