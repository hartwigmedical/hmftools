package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class ContigMappingsTest
{
    @Test
    public void testRoundTrip() throws IOException
    {
        File tempFile = Files.createTempFile("contig_mappings_", ".tsv").toFile();
        tempFile.deleteOnExit();

        List<ContigEntry> entries = List.of(
                new ContigEntry(
                        "ensENSG0000001_GENEA_ENST0000001", "ENSG0000001", "GENEA", "ENST0000001", CHR_1,
                        List.of(new BaseRegion(100, 199), new BaseRegion(300, 399))),
                new ContigEntry(
                        "ensENSG0000002_GENEB_ENST0000002", "ENSG0000002", "GENEB", "ENST0000002", CHR_2,
                        List.of(new BaseRegion(1000, 1050))));

        ContigMappings.write(tempFile.getAbsolutePath(), entries);

        List<ContigEntry> readBack = ContigMappings.read(tempFile.getAbsolutePath());

        assertEquals(2, readBack.size());

        ContigEntry first = readBack.get(0);
        assertEquals("ensENSG0000001_GENEA_ENST0000001", first.contigName());
        assertEquals("ENSG0000001", first.geneId());
        assertEquals("GENEA", first.geneName());
        assertEquals("ENST0000001", first.transName());
        assertEquals(CHR_1, first.chromosome());
        assertEquals(2, first.exonSpans().size());
        assertEquals(new BaseRegion(100, 199), first.exonSpans().get(0));
        assertEquals(new BaseRegion(300, 399), first.exonSpans().get(1));

        ContigEntry second = readBack.get(1);
        assertEquals("ENST0000002", second.transName());
        assertEquals(1, second.exonSpans().size());
        assertEquals(new BaseRegion(1000, 1050), second.exonSpans().get(0));
    }

    @Test
    public void testVaryingSpanCounts() throws IOException
    {
        File tempFile = Files.createTempFile("contig_mappings_", ".tsv").toFile();
        tempFile.deleteOnExit();

        List<ContigEntry> entries = List.of(
                new ContigEntry("ens1", "GENE1", "NAME1", "TRANS1", CHR_1,
                        List.of(new BaseRegion(100, 200))),
                new ContigEntry("ens2", "GENE2", "NAME2", "TRANS2", CHR_1,
                        List.of(new BaseRegion(100, 200), new BaseRegion(500, 600))),
                new ContigEntry("ens3", "GENE3", "NAME3", "TRANS3", CHR_2,
                        List.of(new BaseRegion(10, 50), new BaseRegion(100, 150), new BaseRegion(200, 250), new BaseRegion(300, 400))));

        ContigMappings.write(tempFile.getAbsolutePath(), entries);
        List<ContigEntry> readBack = ContigMappings.read(tempFile.getAbsolutePath());

        assertEquals(3, readBack.size());
        assertEquals(1, readBack.get(0).exonSpans().size());
        assertEquals(2, readBack.get(1).exonSpans().size());
        assertEquals(4, readBack.get(2).exonSpans().size());
    }
}
