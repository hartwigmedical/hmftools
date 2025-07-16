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

        // match
        refVariants.add(TestSomaticVariantDataBuilder.create());
        newVariants.add(TestSomaticVariantDataBuilder.create());

        // value difference
        refVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.alt = "C";
            b.reported = false;
        }));
        newVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.alt = "C";
            b.reported = false;
            b.canonicalHgvsProteinImpact = "p.Val600Gly";
        }));

        // ref only
        refVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.chromosome = "8";
            b.reported = false;
        }));
        newVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.chromosome = "8";
            b.reported = false;
            b.isFromUnfilteredVcf = true;
        }));

        //new only
        newVariants.add(TestSomaticVariantDataBuilder.create(b -> b.type = VariantType.INDEL));

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

        // match
        refVariants.add(TestSomaticVariantDataBuilder.create());
        newVariants.add(TestSomaticVariantDataBuilder.create());

        // value difference
        refVariants.add(TestSomaticVariantDataBuilder.create(b -> b.alt = "C"));
        newVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.alt = "C";
            b.canonicalHgvsProteinImpact = "p.Val600Gly";
        }));

        // ref only
        refVariants.add(TestSomaticVariantDataBuilder.create(b -> b.chromosome = "8"));
        newVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.chromosome = "8";
            b.reported = false;
        }));

        //new only
        refVariants.add(TestSomaticVariantDataBuilder.create(b ->
        {
            b.type = VariantType.INDEL;
            b.isFromUnfilteredVcf = true;
            b.reported = false;
        }));
        newVariants.add(TestSomaticVariantDataBuilder.create(b -> b.type = VariantType.INDEL));

        var result = victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel);
        assertTrue(result);
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }
}
