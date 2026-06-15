package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

public class AltContigPackerTest
{
    private static final int SPACER = 5;

    private static TranscriptContigBuilder.TranscriptContigResult txResult(
            final String transName, final String sequence, final List<BaseRegion> exonSpans)
    {
        return new TranscriptContigBuilder.TranscriptContigResult(
                "G_" + transName, "GENE_" + transName, transName, CHR_1, 1, exonSpans, sequence);
    }

    @Test
    public void testEmptyInputProducesEmptyResult()
    {
        final AltContigPacker packer = new AltContigPacker(SPACER);
        final AltContigPacker.PackResult result = packer.pack(CHR_1, List.of());

        assertEquals("1_tx", result.altContig());
        assertEquals("", result.sequence());
        assertTrue(result.entries().isEmpty());
        assertTrue(result.isEmpty());
    }

    @Test
    public void testSingleTranscriptHasNoSpacer()
    {
        final AltContigPacker packer = new AltContigPacker(SPACER);
        final AltContigPacker.PackResult result = packer.pack(
                CHR_1, List.of(txResult("T1", "ACGTACGT", List.of(new BaseRegion(100, 107)))));

        assertEquals("ACGTACGT", result.sequence());
        assertEquals(1, result.entries().size());

        final ContigEntry entry = result.entries().get(0);
        assertEquals("1_tx", entry.contigName());
        assertEquals(1, entry.altStart());
        assertEquals(8, entry.altEnd());
        assertEquals("T1", entry.transName());
    }

    @Test
    public void testTwoTranscriptsAreSeparatedBySpacer()
    {
        final AltContigPacker packer = new AltContigPacker(SPACER);
        final AltContigPacker.PackResult result = packer.pack(CHR_1, List.of(
                txResult("T1", "AAAA", List.of(new BaseRegion(100, 103))),
                txResult("T2", "CCC", List.of(new BaseRegion(200, 202)))));

        assertEquals("AAAA" + "NNNNN" + "CCC", result.sequence());
        assertEquals(2, result.entries().size());

        final ContigEntry first = result.entries().get(0);
        assertEquals(1, first.altStart());
        assertEquals(4, first.altEnd());

        final ContigEntry second = result.entries().get(1);
        assertEquals(10, second.altStart());
        assertEquals(12, second.altEnd());
    }

    @Test
    public void testThreeTranscriptsHaveContiguousNonOverlappingRanges()
    {
        final AltContigPacker packer = new AltContigPacker(SPACER);
        final AltContigPacker.PackResult result = packer.pack(CHR_1, List.of(
                txResult("T1", "AAAAAA", List.of(new BaseRegion(100, 105))),
                txResult("T2", "CCC", List.of(new BaseRegion(200, 202))),
                txResult("T3", "GGGGGGGG", List.of(new BaseRegion(300, 307)))));

        assertEquals(6 + 5 + 3 + 5 + 8, result.sequence().length());
        for(int i = 0; i < result.entries().size() - 1; ++i)
        {
            final ContigEntry a = result.entries().get(i);
            final ContigEntry b = result.entries().get(i + 1);
            assertEquals(a.altEnd() + SPACER + 1, b.altStart());
        }
    }

    @Test
    public void testEntryMetadataIsCarriedThrough()
    {
        final AltContigPacker packer = new AltContigPacker(SPACER);
        final AltContigPacker.PackResult result = packer.pack(
                CHR_1, List.of(txResult("T1", "ACGT", List.of(new BaseRegion(100, 103)))));

        final ContigEntry entry = result.entries().get(0);
        assertEquals("G_T1", entry.geneId());
        assertEquals("GENE_T1", entry.geneName());
        assertEquals(CHR_1, entry.chromosome());
        assertEquals(List.of(new BaseRegion(100, 103)), entry.exonSpans());
    }
}
