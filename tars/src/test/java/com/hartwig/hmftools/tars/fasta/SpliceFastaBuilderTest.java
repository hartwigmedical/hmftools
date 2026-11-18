package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.fasta.SpliceFastaBuilder.packChromosomeContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

import org.junit.Test;

public class SpliceFastaBuilderTest
{
    private static final int SPACER = 5;

    private static TranscriptContigBuilder.TranscriptContigResult txResult(
            final String transName, final String sequence, final List<BaseRegion> exonSpans)
    {
        return new TranscriptContigBuilder.TranscriptContigResult(
                "G_" + transName, "GENE_" + transName, transName, CHR_1, 1, exonSpans, sequence);
    }

    @Test
    public void testEmptyInputProducesEmptySequence()
    {
        List<ContigEntry> entries = new ArrayList<>();
        String sequence = packChromosomeContig("1_tx", List.of(), entries, SPACER);

        assertEquals("", sequence);
        assertTrue(entries.isEmpty());
    }

    @Test
    public void testSingleTranscriptHasNoSpacer()
    {
        List<ContigEntry> entries = new ArrayList<>();
        String sequence = packChromosomeContig(
                "1_tx", List.of(txResult("T1", "ACGTACGT", List.of(new BaseRegion(100, 107)))), entries, SPACER);

        assertEquals("ACGTACGT", sequence);
        assertEquals(1, entries.size());

        ContigEntry entry = entries.get(0);
        assertEquals("1_tx", entry.contigName());
        assertEquals(1, entry.contigStart());
        assertEquals(8, entry.contigEnd());
        assertEquals("T1", entry.transName());
    }

    @Test
    public void testTwoTranscriptsAreSeparatedBySpacer()
    {
        List<ContigEntry> entries = new ArrayList<>();
        String sequence = packChromosomeContig(
                "1_tx", List.of(
                        txResult("T1", "AAAA", List.of(new BaseRegion(100, 103))),
                        txResult("T2", "CCC", List.of(new BaseRegion(200, 202)))), entries, SPACER);

        assertEquals("AAAA" + "NNNNN" + "CCC", sequence);
        assertEquals(2, entries.size());

        assertEquals(1, entries.get(0).contigStart());
        assertEquals(4, entries.get(0).contigEnd());
        assertEquals(10, entries.get(1).contigStart());
        assertEquals(12, entries.get(1).contigEnd());
    }

    @Test
    public void testEntryMetadataIsCarriedThrough()
    {
        List<ContigEntry> entries = new ArrayList<>();
        packChromosomeContig("1_tx", List.of(txResult("T1", "ACGT", List.of(new BaseRegion(100, 103)))), entries, SPACER);

        ContigEntry entry = entries.get(0);
        assertEquals("G_T1", entry.geneId());
        assertEquals("GENE_T1", entry.geneName());
        assertEquals(CHR_1, entry.chromosome());
        assertEquals(List.of(new BaseRegion(100, 103)), entry.exonSpans());
    }
}
