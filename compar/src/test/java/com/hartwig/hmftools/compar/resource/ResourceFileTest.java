package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.compar.resource.ResourceTestUtils.tsv;
import static com.hartwig.hmftools.compar.resource.ResourceTestUtils.writeFile;
import static com.hartwig.hmftools.compar.resource.ResourceType.CLINVAR;
import static com.hartwig.hmftools.compar.resource.ResourceType.DRIVER_GENE_PANEL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ResourceFileTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    @Test
    public void rejectsFilesThatCannotBeCompared() throws IOException
    {
        assertEquals("missing header", readError(DRIVER_GENE_PANEL, "DriverGenePanel.38.tsv"));

        assertEquals("missing header", readError(CLINVAR, "clinvar.38.vcf.gz", "##fileformat=VCFv4.1"));

        assertEquals("missing columns(gene)", readError(DRIVER_GENE_PANEL, "DriverGenePanel.38.tsv", tsv("name", "reportMissense")));

        assertEquals(
                "line(2) has 1 of 2 columns",
                readError(DRIVER_GENE_PANEL, "DriverGenePanel.38.tsv", tsv("gene", "reportMissense"), "KRAS"));

        assertEquals(
                "line(3) has 2 of 8 columns",
                readError(CLINVAR, "clinvar.38.vcf.gz", "##fileformat=VCFv4.1",
                        tsv("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"), tsv("1", "100")));

        File notGzip = temp.newFolder().toPath().resolve("clinvar.38.vcf.gz").toFile();
        Files.writeString(notGzip.toPath(), "not gzip");

        assertThrows(IOException.class, () -> new ResourceFile(ResourceDefinition.create(CLINVAR), notGzip.toPath()));
    }

    private String readError(final ResourceType type, final String fileName, final String... lines) throws IOException
    {
        File file = writeFile(temp.newFolder(), fileName, lines);
        IOException exception = assertThrows(IOException.class, () -> new ResourceFile(ResourceDefinition.create(type), file.toPath()));
        return exception.getMessage();
    }
}
