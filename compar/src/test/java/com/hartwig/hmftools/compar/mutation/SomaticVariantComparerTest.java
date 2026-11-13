
package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.CommonUtils.compareAbsentEntries;
import static com.hartwig.hmftools.compar.mutation.VariantData.FLD_GENE;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.ExpectedCache;
import com.hartwig.hmftools.compar.common.ExpectedValue;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class SomaticVariantComparerTest
{
    @Rule
    public TemporaryFolder TempDir = new TemporaryFolder();

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
    public void absentGeneCountsPassCalls() throws IOException
    {
        SomaticVariantData brafCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.reported = false;
            b.canonicalCodingEffect = "NONE";
        });

        SomaticVariantData filteredBrafCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.reported = false;
            b.filters = Set.of("LOW_VAF");
        });

        SomaticVariantData krasCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "12";
            b.gene = "KRAS";
        });

        List<Mismatch> mismatches = compareWithAbsent(
                "BRAF", List.of(), List.of(brafCall, filteredBrafCall, krasCall), MatchLevel.DETAILED);

        assertEquals(List.of(MismatchType.ABSENT_FOUND, MismatchType.NEW_ONLY), mismatches.stream().map(x -> x.Type).toList());
        assertSame(brafCall, mismatches.get(0).NewItem);
        assertEquals("BRAF", mismatches.get(0).OldItem.fieldValues().get(FLD_GENE).displayValue());
        assertSame(krasCall, mismatches.get(1).NewItem);
    }

    @Test
    public void absentGeneAtReportableOnlyCountsReportedCalls() throws IOException
    {
        SomaticVariantData unreportedCall = TestSomaticVariantDataBuilder.BUILDER.create(b -> b.reported = false);
        SomaticVariantData reportedCall = TestSomaticVariantDataBuilder.BUILDER.create(b -> b.position = 140453140);

        List<Mismatch> mismatches = compareWithAbsent("BRAF", List.of(), List.of(unreportedCall, reportedCall), MatchLevel.REPORTABLE);

        assertEquals(List.of(MismatchType.ABSENT_FOUND), mismatches.stream().map(x -> x.Type).toList());
        assertSame(reportedCall, mismatches.get(0).NewItem);
    }

    @Test
    public void absentGeneWithoutCountedCallIsConfirmed() throws IOException
    {
        SomaticVariantData unreportedBrafCall = TestSomaticVariantDataBuilder.BUILDER.create(b -> b.reported = false);

        SomaticVariantData krasCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.chromosome = "12";
            b.gene = "KRAS";
        });

        List<Mismatch> mismatches = compareWithAbsent(
                "BRAF", List.of(), List.of(unreportedBrafCall, krasCall), MatchLevel.REPORTABLE);

        assertEquals(List.of(MismatchType.NEW_ONLY, MismatchType.ABSENT_CONFIRMED), mismatches.stream().map(x -> x.Type).toList());
        assertSame(krasCall, mismatches.get(0).NewItem);
        assertEquals("BRAF", mismatches.get(1).OldItem.key());
        assertNull(mismatches.get(1).NewItem);
    }

    @Test
    public void listedVariantTakesPriorityOverAbsentGene() throws IOException
    {
        SomaticVariantData listedCall = TestSomaticVariantDataBuilder.BUILDER.create();
        SomaticVariantData unlistedCall = TestSomaticVariantDataBuilder.BUILDER.create(b -> b.position = 140453140);

        List<SomaticVariantData> oldVariants = List.of(TestSomaticVariantDataBuilder.BUILDER.create());
        List<Mismatch> mismatches = compareWithAbsent("BRAF", oldVariants, List.of(listedCall, unlistedCall), MatchLevel.REPORTABLE);

        assertEquals(List.of(MismatchType.ABSENT_FOUND), mismatches.stream().map(x -> x.Type).toList());
        assertSame(unlistedCall, mismatches.get(0).NewItem);
    }

    @Test
    public void proteinKeyMatchesGeneAndProteinChange()
    {
        SomaticVariantComparer victim = new SomaticVariantComparer(new ComparConfig(), Collections.emptyMap());

        SomaticVariantData brafRow = SomaticVariantData.fromExpected(
                List.of(new ExpectedValue("BRAF:p.V600E", "AF", "0.10", false)), victim.fieldsList());

        SomaticVariantData egfrRow = SomaticVariantData.fromExpected(
                List.of(new ExpectedValue("EGFR:p.L858R", "AF", "0.10", false)), victim.fieldsList());

        SomaticVariantData v600kCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.alt = "C";
            b.canonicalHgvsProteinImpact = "p.Val600Lys";
        });

        SomaticVariantData v600eCall = TestSomaticVariantDataBuilder.BUILDER.create();

        List<Mismatch> mismatches = new ArrayList<>();
        assertTrue(victim.identifyMismatches(
                "TEST", mismatches, List.of(brafRow, egfrRow), List.of(v600kCall, v600eCall), MatchLevel.REPORTABLE));

        assertEquals(
                List.of(MismatchType.VALUE, MismatchType.OLD_ONLY, MismatchType.NEW_ONLY), mismatches.stream().map(x -> x.Type).toList());
        assertEquals("BRAF p.Val600Glu", mismatches.get(0).OldItem.key());
        assertSame(v600eCall, mismatches.get(0).NewItem);
        assertEquals(List.of("AF(0.10/0.45)"), mismatches.get(0).DiffValues);
        assertEquals("EGFR p.Leu858Arg", mismatches.get(1).OldItem.key());
        assertSame(v600kCall, mismatches.get(2).NewItem);
    }

    @Test
    public void expectedItemsCarryTheirRowRemarks()
    {
        SomaticVariantComparer victim = new SomaticVariantComparer(new ComparConfig(), Collections.emptyMap());

        Map<String,List<ExpectedValue>> valuesByKey = new TreeMap<>(Map.of(
                "7:140453136:A:T", List.of(
                        new ExpectedValue("7:140453136:A:T", "AF", "0.45", false, "confirmed by ddPCR"),
                        new ExpectedValue("7:140453136:A:T", "Reported", "true", false, "confirmed by ddPCR"),
                        new ExpectedValue("7:140453136:A:T", "Gene", "BRAF", false, "second assay")),
                "EGFR:p.L858R", List.of(new ExpectedValue("EGFR:p.L858R", "AF", "0.10", false)),
                "KRAS", List.of(new ExpectedValue("KRAS", "Status", "", true, "wild type by panel"))));

        List<String> remarks = victim.loadFromExpected(valuesByKey).stream().map(ComparableItem::remark).toList();

        assertEquals(List.of("confirmed by ddPCR;second assay", "", "wild type by panel"), remarks);
    }

    @Test
    public void codingKeyIsInvalid()
    {
        SomaticVariantComparer victim = new SomaticVariantComparer(new ComparConfig(), Collections.emptyMap());
        assertNull(SomaticVariantData.fromExpected(List.of(new ExpectedValue("BRAF:c.1799T>A", "AF", "0.10", false)), victim.fieldsList()));
    }

    @Test
    public void absentProteinKeyCountsOnlyThatChange() throws IOException
    {
        SomaticVariantData v600eCall = TestSomaticVariantDataBuilder.BUILDER.create();

        SomaticVariantData v600kCall = TestSomaticVariantDataBuilder.BUILDER.create(b ->
        {
            b.alt = "C";
            b.canonicalHgvsProteinImpact = "p.Val600Lys";
        });

        List<Mismatch> mismatches = compareWithAbsent("BRAF:p.V600E", List.of(), List.of(v600eCall, v600kCall), MatchLevel.DETAILED);

        assertEquals(List.of(MismatchType.ABSENT_FOUND, MismatchType.NEW_ONLY), mismatches.stream().map(x -> x.Type).toList());
        assertSame(v600eCall, mismatches.get(0).NewItem);
        assertTrue(mismatches.get(0).OldItem.fieldValues().isEmpty());
        assertSame(v600kCall, mismatches.get(1).NewItem);
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

    private List<Mismatch> compareWithAbsent(
            final String absentKey, final List<SomaticVariantData> oldVariants, final List<SomaticVariantData> newVariants,
            final MatchLevel matchLevel) throws IOException
    {
        File expectedFile = TempDir.newFile();
        Files.write(
                expectedFile.toPath(), List.of(
                        String.join(TSV_DELIM, "Category", "Key", "Status", "FieldValues"),
                        String.join(TSV_DELIM, "SOMATIC_VARIANT", absentKey, "Absent", "")));

        ExpectedCache expectedCache = new ExpectedCache();
        expectedCache.loadFiles(expectedFile.getAbsolutePath());

        ComparConfig config = new ComparConfig();
        config.Sources.set(0, SourceData.fromExpectedSource(SourceType.OLD, expectedCache, null));

        SomaticVariantComparer victim = new SomaticVariantComparer(config, Collections.emptyMap());
        List<Mismatch> mismatches = new ArrayList<>();
        assertTrue(victim.identifyMismatches("TEST", mismatches, oldVariants, newVariants, matchLevel));

        compareAbsentEntries(victim, "TEST", mismatches);
        return mismatches;
    }
}
