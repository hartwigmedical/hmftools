package com.hartwig.hmftools.compar.mutation;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.ArrayList;

import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

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

        var baseVariant = SomaticVariantDataTestFactory.createDefault();

        // match
        refVariants.add(baseVariant.build());
        newVariants.add(baseVariant.build());

        // value difference
        refVariants.add(baseVariant.withAlt("C").withReported(false).build());
        newVariants.add(baseVariant.withAlt("C").withCanonicalHgvsProteinImpact("p.Val600Gly").withReported(false).build());

        // ref only
        refVariants.add(baseVariant.withChromosome("8").withReported(false).build());
        newVariants.add(baseVariant.withChromosome("8").withReported(false).withIsFromUnfilteredVcf(true).build());

        //new only
        newVariants.add(baseVariant.withType(VariantType.INDEL).build());

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

        var baseVariant = SomaticVariantDataTestFactory.createDefault();

        // match
        refVariants.add(baseVariant.build());
        newVariants.add(baseVariant.build());

        // value difference
        refVariants.add(baseVariant.withAlt("C").build());
        newVariants.add(baseVariant.withAlt("C").withCanonicalHgvsProteinImpact("p.Val600Gly").build());

        // ref only
        refVariants.add(baseVariant.withChromosome("8").build());
        newVariants.add(baseVariant.withChromosome("8").withReported(false).build());

        //new only
        refVariants.add(baseVariant.withType(VariantType.INDEL).withIsFromUnfilteredVcf(true).withReported(false).build());
        newVariants.add(baseVariant.withType(VariantType.INDEL).build());

        var result = victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel);
        assertTrue(result);
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }
}
