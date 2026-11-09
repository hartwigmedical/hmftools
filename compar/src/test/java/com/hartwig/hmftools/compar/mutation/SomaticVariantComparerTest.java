
package com.hartwig.hmftools.compar.mutation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.ExpectedCache;
import com.hartwig.hmftools.compar.common.ExpectedValue;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;

import org.junit.Test;

public class SomaticVariantComparerTest
{
    @Test
    public void emptyComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.DETAILED;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();

        assertTrue(victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel));
        assertTrue(mismatches.isEmpty());
    }

    @Test
    public void simpleDetailedComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.DETAILED;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();

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
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.OLD_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }

    @Test
    public void simpleReportableComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.REPORTABLE;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();

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
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.OLD_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.NEW_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.VALUE).count());
        assertEquals(3, mismatches.size());
    }

    @Test
    public void uniqueOnChromosomeComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.DETAILED;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();

        // match
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());

        // ref only
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.chromosome = "12"));

        //new only
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create(b -> b.chromosome = "13"));

        assertTrue(victim.identifyMismatches(sampleId, mismatches, refVariants, newVariants, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.OLD_ONLY).count());
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.NEW_ONLY).count());
        assertEquals(2, mismatches.size());
    }

    @Test
    public void bothInvalidComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.DETAILED;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();

        assertFalse(victim.identifyMismatches(sampleId, mismatches, null, null, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.INVALID_BOTH).count());
        assertEquals(1, mismatches.size());
    }

    @Test
    public void refInvalidComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.DETAILED;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> newVariants = new ArrayList<>();

        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        newVariants.add(TestSomaticVariantDataBuilder.BUILDER.createWithAlternateDefaults());

        assertFalse(victim.identifyMismatches(sampleId, mismatches, null, newVariants, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.INVALID_OLD).count());
        assertEquals(1, mismatches.size());
    }

    @Test
    public void newInvalidComparison()
    {
        ComparConfig config = new ComparConfig();
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        MatchLevel matchLevel = MatchLevel.DETAILED;

        String sampleId = "TEST";
        List<Mismatch> mismatches = new ArrayList<>();
        List<SomaticVariantData> refVariants = new ArrayList<>();

        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.create());
        refVariants.add(TestSomaticVariantDataBuilder.BUILDER.createWithAlternateDefaults());

        assertFalse(victim.identifyMismatches(sampleId, mismatches, refVariants, null, matchLevel));
        assertEquals(1, mismatches.stream().filter(x -> x.Type == MismatchType.INVALID_NEW).count());
        assertEquals(1, mismatches.size());
    }

    @Test
    public void expectedVariantWithoutCallIsOldOnlyAtReportable()
    {
        ComparConfig config = new ComparConfig();
        config.Sources.set(0, SourceData.fromExpectedSource(SourceType.OLD, new ExpectedCache(), null));
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        SomaticVariantData expectedRow = SomaticVariantData.fromExpected(
                List.of(new ExpectedValue("7:140453136:A:T", "AF", "0.45", false)), victim.fieldsList());

        List<Mismatch> mismatches = new ArrayList<>();
        assertTrue(victim.identifyMismatches("TEST", mismatches, List.of(expectedRow), List.of(), MatchLevel.REPORTABLE));

        assertEquals(List.of(MismatchType.OLD_ONLY), mismatches.stream().map(x -> x.Type).toList());
    }

    @Test
    public void filteredCallForExpectedVariantIsOldOnly()
    {
        ComparConfig config = new ComparConfig();
        config.Sources.set(0, SourceData.fromExpectedSource(SourceType.OLD, new ExpectedCache(), null));
        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());

        SomaticVariantData expectedRow = SomaticVariantData.fromExpected(
                List.of(new ExpectedValue("7:140453136:A:T", "AF", "0.45", false)), victim.fieldsList());

        SomaticVariantData ponCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.reported = false;
            b.filters = Set.of("PON");
        });

        for(MatchLevel matchLevel : List.of(MatchLevel.REPORTABLE, MatchLevel.DETAILED))
        {
            List<Mismatch> mismatches = new ArrayList<>();
            assertTrue(victim.identifyMismatches("TEST", mismatches, List.of(expectedRow), List.of(ponCall), matchLevel));

            assertEquals(List.of(MismatchType.OLD_ONLY), mismatches.stream().map(x -> x.Type).toList());
            assertEquals(List.of("Filter(PASS/PON)"), mismatches.get(0).DiffValues);
        }
    }

    @Test
    public void unfilteredVcfRecordForExpectedVariantIsOldOnly()
    {
        ComparConfig config = new ComparConfig();
        config.Sources.set(0, SourceData.fromExpectedSource(SourceType.OLD, new ExpectedCache(), null));

        SomaticVariantData sageRecord = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.reported = false;
            b.isFromUnfilteredVcf = true;
        });

        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap())
        {
            @Override
            protected SomaticVariantData findUnfilteredVariant(
                    final SomaticVariantData testVariant, final SourceType sourceType, final String sourceSampleId)
            {
                return sageRecord;
            }
        };

        SomaticVariantData expectedRow = SomaticVariantData.fromExpected(
                List.of(new ExpectedValue("7:140453136:A:T", "AF", "0.45", false)), victim.fieldsList());

        List<Mismatch> mismatches = new ArrayList<>();
        assertTrue(victim.identifyMismatches("TEST", mismatches, List.of(expectedRow), List.of(), MatchLevel.DETAILED));

        assertEquals(List.of(MismatchType.OLD_ONLY), mismatches.stream().map(x -> x.Type).toList());
        assertSame(sageRecord, mismatches.get(0).NewItem);
    }
}
