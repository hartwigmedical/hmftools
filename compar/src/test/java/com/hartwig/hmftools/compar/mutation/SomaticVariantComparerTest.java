package com.hartwig.hmftools.compar.mutation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.ArrayList;
import java.util.List;

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
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;

        assertTrue(victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel));
        assertTrue(mismatches.isEmpty());
    }

    @Test
    public void simpleDetailedComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;

        // match
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());

        // value difference
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.alt = "C";
            b.reported = false;
        }));
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.alt = "C";
            b.reported = false;
            b.canonicalHgvsProteinImpact = "p.Val600Gly";
        }));

        // ref only
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "8";
            b.reported = false;
        }));
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "8";
            b.reported = false;
            b.isFromUnfilteredVcf = true;
        }));

        //new only
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.type = VariantType.INDEL));

        assertTrue(victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }

    @Test
    public void simpleReportableComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.REPORTABLE;

        // match
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());

        // value difference
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.alt = "C"));
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.alt = "C";
            b.canonicalHgvsProteinImpact = "p.Val600Gly";
        }));

        // ref only
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.chromosome = "8"));
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "8";
            b.reported = false;
        }));

        //new only
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.type = VariantType.INDEL;
            b.isFromUnfilteredVcf = true;
            b.reported = false;
        }));
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.type = VariantType.INDEL));

        assertTrue(victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }

    @Test
    public void uniqueOnChromosomeComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;

        // match
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());

        // ref only
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.chromosome = "12"));

        //new only
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.chromosome = "13"));

        assertTrue(victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.REF_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.NEW_ONLY).count());
        assertEquals(2, mismatches.size());
    }

    @Test
    public void bothInvalidComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;

        assertFalse(victim.identifyMismatches(sampleId, mismatches, null, null, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.INVALID_BOTH).count());
        assertEquals(1, mismatches.size());
    }

    @Test
    public void refInvalidComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;

        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.createWithAlternateDefaults());

        assertFalse(victim.identifyMismatches(sampleId, mismatches, null, newVariants, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.INVALID_REF).count());
        assertEquals(1, mismatches.size());
    }

    @Test
    public void newInvalidComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config);

        DiffThresholds diffThresholds = new DiffThresholds();
        victim.registerThresholds(diffThresholds);

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;

        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.createWithAlternateDefaults());

        assertFalse(victim.identifyMismatches(sampleId, mismatches, refVariants, null, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.MismatchType() == MismatchType.INVALID_NEW).count());
        assertEquals(1, mismatches.size());
    }
}
