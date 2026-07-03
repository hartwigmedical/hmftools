package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BaseRegion;

import org.junit.Test;

public class ContigSidecarTest
{
    @Test
    public void testRoundTrip() throws IOException
    {
        File tempFile = Files.createTempFile("contig_sidecar_", ".tsv").toFile();
        tempFile.deleteOnExit();

        List<ContigEntry> entries = List.of(
                new ContigEntry(
                        "chr1_tx", 1, 200, "ENSG0000001", "GENEA", "ENST0000001", CHR_1, 1,
                        List.of(new BaseRegion(100, 199), new BaseRegion(300, 399))),
                new ContigEntry(
                        "chr2_tx", 1, 51, "ENSG0000002", "GENEB", "ENST0000002", CHR_2, -1,
                        List.of(new BaseRegion(1000, 1050))));

        ContigSidecar.write(tempFile.getAbsolutePath(), entries);

        List<ContigEntry> readBack = ContigSidecar.read(tempFile.getAbsolutePath());

        assertEquals(2, readBack.size());

        ContigEntry first = readBack.get(0);
        assertEquals("chr1_tx", first.contigName());
        assertEquals(1, first.altStart());
        assertEquals(200, first.altEnd());
        assertEquals("ENSG0000001", first.geneId());
        assertEquals("GENEA", first.geneName());
        assertEquals("ENST0000001", first.transName());
        assertEquals(CHR_1, first.chromosome());
        assertEquals(1, first.strand());
        assertEquals(2, first.exonSpans().size());
        assertEquals(new BaseRegion(100, 199), first.exonSpans().get(0));
        assertEquals(new BaseRegion(300, 399), first.exonSpans().get(1));

        ContigEntry second = readBack.get(1);
        assertEquals("chr2_tx", second.contigName());
        assertEquals(1, second.altStart());
        assertEquals(51, second.altEnd());
        assertEquals("ENST0000002", second.transName());
        assertEquals(-1, second.strand());
        assertEquals(1, second.exonSpans().size());
        assertEquals(new BaseRegion(1000, 1050), second.exonSpans().get(0));
    }

    @Test
    public void testHeaderOnlySidecarReadsAsEmpty() throws IOException
    {
        // a sidecar with only its header row (no entries) loads as an empty list rather than failing
        File tempFile = Files.createTempFile("contig_sidecar_empty_", ".tsv").toFile();
        tempFile.deleteOnExit();

        String header = Arrays.stream(ContigSidecar.Column.values())
                .map(Enum::name)
                .collect(Collectors.joining(TSV_DELIM));
        Files.writeString(tempFile.toPath(), header + "\n");

        List<ContigEntry> readBack = ContigSidecar.read(tempFile.getAbsolutePath());
        assertTrue(readBack.isEmpty());
    }
}
