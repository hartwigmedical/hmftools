package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverType;
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

public class GeneCopyNumberComparerTest
{
    @Rule
    public TemporaryFolder TempDir = new TemporaryFolder();

    @Test
    public void expectedCopyNumbersAreComparedOnGivenFields() throws IOException
    {
        File expectedFile = TempDir.newFile();
        Files.write(
                expectedFile.toPath(), List.of(
                        String.join(TSV_DELIM, "Category", "Key", "MinCopyNumber", "MaxCopyNumber", "DriverType", "Remark"),
                        String.join(TSV_DELIM, "GENE_COPY_NUMBER", "ERBB2", "N/A", "8", "AMP", "amplified by FISH"),
                        String.join(TSV_DELIM, "GENE_COPY_NUMBER", "CDKN2A", "0", "N/A", "N/A", ""),
                        String.join(TSV_DELIM, "GENE_COPY_NUMBER", "PTEN", "N/A", "N/A", "DEL", ""),
                        String.join(TSV_DELIM, "GENE_COPY_NUMBER", "MYC", "N/A", "6", "N/A", "")));

        ExpectedCache expectedCache = new ExpectedCache();
        expectedCache.loadFiles(expectedFile.getAbsolutePath());

        Map<String, List<ExpectedValue>> expectedValuesByKey = expectedCache.sampleExpectedEntries("TEST", GENE_COPY_NUMBER).stream()
                .collect(Collectors.groupingBy(x -> x.Key));

        ComparConfig config = new ComparConfig();
        config.Sources.set(0, SourceData.fromExpectedSource(SourceType.OLD, expectedCache, null));
        GeneCopyNumberComparer comparer = new GeneCopyNumberComparer(config, Collections.emptyMap());
        List<ComparableItem> expectedItems = comparer.loadFromExpected(expectedValuesByKey);

        List<ComparableItem> newItems = Lists.newArrayList(
                TestGeneCopyNumberDataBuilder.BUILDER.create(b ->
                {
                    b.gene = "ERBB2";
                    b.maxCopyNumber = 8.2;
                }),
                TestGeneCopyNumberDataBuilder.BUILDER.create(b ->
                {
                    b.gene = "CDKN2A";
                    b.minCopyNumber = 1.5;
                }),
                TestGeneCopyNumberDataBuilder.BUILDER.create(b ->
                {
                    b.gene = "PTEN";
                    b.driverType = DriverType.HET_DEL;
                }),
                TestGeneCopyNumberDataBuilder.BUILDER.create(b -> b.gene = "KRAS"));

        List<Mismatch> mismatches = Lists.newArrayList();
        comparer.compareItems(mismatches, MatchLevel.REPORTABLE, true, true, expectedItems, newItems);

        Map<String, Mismatch> mismatchesByGene = mismatches.stream()
                .collect(Collectors.toMap(x -> x.OldItem != null ? x.OldItem.key() : x.NewItem.key(), x -> x));

        assertEquals(4, mismatchesByGene.size());
        assertEquals(MismatchType.FULL_MATCH, mismatchesByGene.get("ERBB2").Type);
        assertEquals(MismatchType.VALUE, mismatchesByGene.get("CDKN2A").Type);
        assertEquals(1, mismatchesByGene.get("CDKN2A").DiffValues.size());
        assertTrue(mismatchesByGene.get("CDKN2A").DiffValues.get(0).startsWith("MinCopyNumber("));
        assertEquals(MismatchType.VALUE, mismatchesByGene.get("PTEN").Type);
        assertTrue(mismatchesByGene.get("PTEN").DiffValues.get(0).startsWith("DriverType("));
        assertEquals(MismatchType.OLD_ONLY, mismatchesByGene.get("MYC").Type);
    }
}
