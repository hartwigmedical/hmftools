package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;
import static com.hartwig.hmftools.compar.common.CategoryType.RNA_FUSION;
import static com.hartwig.hmftools.compar.common.CategoryType.SOMATIC_VARIANT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ExpectedCacheTest
{
    @Rule
    public TemporaryFolder TempDir = new TemporaryFolder();

    @Test
    public void invalidStatusOrAbsentKeyIsRejected()
    {
        ExpectedCache.checkStatus(SOMATIC_VARIANT, "", "7:140453136:A:T");
        ExpectedCache.checkStatus(SOMATIC_VARIANT, "Present", "7:140453136:A:T");
        ExpectedCache.checkStatus(SOMATIC_VARIANT, "absent", "BRAF");
        ExpectedCache.checkStatus(RNA_FUSION, "Absent", "TMPRSS2:ERG");

        assertThrows(IllegalArgumentException.class, () -> ExpectedCache.checkStatus(SOMATIC_VARIANT, "Absnt", "BRAF"));
        assertThrows(IllegalArgumentException.class, () -> ExpectedCache.checkStatus(SOMATIC_VARIANT, "Absent ", "BRAF"));
        assertThrows(IllegalArgumentException.class, () -> ExpectedCache.checkStatus(PURITY, "Absent", "PURITY"));
        assertThrows(IllegalArgumentException.class, () -> ExpectedCache.checkStatus(SOMATIC_VARIANT, "Absent", ""));
        assertThrows(IllegalArgumentException.class, () -> ExpectedCache.checkStatus(SOMATIC_VARIANT, "Absent", "7:140453136:A:T"));
        assertThrows(IllegalArgumentException.class, () -> ExpectedCache.checkStatus(RNA_FUSION, "Absent", "EML4:ALK:2:3000:2:4000"));
    }

    @Test
    public void presentAndAbsentRowsAreReturnedSeparately() throws IOException
    {
        File expectedFile = TempDir.newFile();
        Files.write(
                expectedFile.toPath(), List.of(
                        String.join(TSV_DELIM, "Category", "Key", "Status", "FieldValues"),
                        String.join(TSV_DELIM, "SOMATIC_VARIANT", "7:140453136:A:T", "Present", "AF=0.3"),
                        String.join(TSV_DELIM, "SOMATIC_VARIANT", "KRAS", "Absent", "")));

        ExpectedCache expectedCache = new ExpectedCache();
        expectedCache.loadFiles(expectedFile.getAbsolutePath());

        List<ExpectedValue> presentValues = expectedCache.sampleExpectedEntries("TEST", SOMATIC_VARIANT);
        List<ExpectedValue> absentValues = expectedCache.sampleAbsentEntries("TEST", SOMATIC_VARIANT);

        assertEquals(1, presentValues.size());
        assertEquals("7:140453136:A:T", presentValues.get(0).Key);
        assertEquals(1, absentValues.size());
        assertEquals("KRAS", absentValues.get(0).Key);
        assertTrue(expectedCache.sampleAbsentEntries("TEST", RNA_FUSION).isEmpty());
    }

    @Test
    public void remarkIsKeptAsNotesNotAsField() throws IOException
    {
        File columnsFile = TempDir.newFile();
        Files.write(
                columnsFile.toPath(), List.of(
                        String.join(TSV_DELIM, "Category", "Key", "Status", "AF", "Remark"),
                        String.join(TSV_DELIM, "SOMATIC_VARIANT", "7:140453136:A:T", "Present", "0.3", "confirmed by ddPCR")));

        File fieldValuesFile = TempDir.newFile();
        Files.write(
                fieldValuesFile.toPath(), List.of(
                        String.join(TSV_DELIM, "Category", "Key", "Status", "FieldValues", "Remark"),
                        String.join(TSV_DELIM, "SOMATIC_VARIANT", "7:140453136:A:T", "Present", "AF=0.3", "confirmed by ddPCR")));

        for(File expectedFile : List.of(columnsFile, fieldValuesFile))
        {
            ExpectedCache expectedCache = new ExpectedCache();
            expectedCache.loadFiles(expectedFile.getAbsolutePath());

            List<ExpectedValue> values = expectedCache.sampleExpectedEntries("TEST", SOMATIC_VARIANT);

            assertEquals(1, values.size());
            assertEquals("AF", values.get(0).FieldName);
            assertEquals("0.3", values.get(0).Value);
            assertEquals("confirmed by ddPCR", values.get(0).Remark);
        }
    }
}
