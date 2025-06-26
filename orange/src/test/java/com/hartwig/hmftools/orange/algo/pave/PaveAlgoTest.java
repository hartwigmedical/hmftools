package com.hartwig.hmftools.orange.algo.pave;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class PaveAlgoTest
{
    private static final String MATCHING_GENE = "gene";
    private static final String MATCHING_TRANSCRIPT = "trans";

    @Test
    public void canCreatePaveEntries()
    {
        // "gene", coding range: 300-800, reverse strand. exon1: 600-900, exon2: 450-550, exon3: 250:350
        PaveAlgo pave = new PaveAlgo(TestEnsemblDataCacheFactory.loadTestCache(), false);

        // Missing gene
        assertNull(pave.run("other gene", MATCHING_TRANSCRIPT, 0));

        // Missing transcript
        assertNull(pave.run(MATCHING_GENE, "other trans", 0));

        PaveEntry entryBeyond = pave.run(MATCHING_GENE, MATCHING_TRANSCRIPT, 1000);
        assertNull(entryBeyond.affectedCodon());
        assertNull(entryBeyond.affectedExon());

        PaveEntry entryUtr5 = pave.run(MATCHING_GENE, MATCHING_TRANSCRIPT, 850);
        assertNull(entryUtr5.affectedCodon());
        assertEquals(1, (int) entryUtr5.affectedExon());

        PaveEntry entryCoding = pave.run(MATCHING_GENE, MATCHING_TRANSCRIPT, 500);
        assertEquals(84, (int) entryCoding.affectedCodon());
        assertEquals(2, (int) entryCoding.affectedExon());

        PaveEntry entrySplice = pave.run(MATCHING_GENE, MATCHING_TRANSCRIPT, 445);
        assertNull(entrySplice.affectedCodon());
        assertEquals(2, (int) entrySplice.affectedExon());

        PaveEntry entryIntronic = pave.run(MATCHING_GENE, MATCHING_TRANSCRIPT, 400);
        assertNull(entryIntronic.affectedCodon());
        assertNull(entryIntronic.affectedExon());

        PaveEntry entryUtr3 = pave.run(MATCHING_GENE, MATCHING_TRANSCRIPT, 250);
        assertNull(entryUtr3.affectedCodon());
        assertEquals(3, (int) entryUtr3.affectedExon());
    }

    @Test
    public void canFindAffectedExon()
    {
        ExonData exon1 = createExon(1, 30, 50);
        ExonData exon2 = createExon(2, 110, 150);
        ExonData exon3 = createExon(3, 180, 200);
        List<ExonData> exons = Lists.newArrayList(exon1, exon2, exon3);

        assertNull(PaveAlgo.findAffectedExon(exons, 10));
        assertEquals(exon1, PaveAlgo.findAffectedExon(exons, 25));
        assertEquals(exon1, PaveAlgo.findAffectedExon(exons, 30));
        assertEquals(exon1, PaveAlgo.findAffectedExon(exons, 58));
        assertNull(PaveAlgo.findAffectedExon(exons, 80));
        assertEquals(exon2, PaveAlgo.findAffectedExon(exons, 105));
        assertEquals(exon2, PaveAlgo.findAffectedExon(exons, 150));
        assertEquals(exon2, PaveAlgo.findAffectedExon(exons, 155));
        assertNull(PaveAlgo.findAffectedExon(exons, 170));
        assertEquals(exon3, PaveAlgo.findAffectedExon(exons, 175));
        assertEquals(exon3, PaveAlgo.findAffectedExon(exons, 190));
        assertNull(PaveAlgo.findAffectedExon(exons, 250));
    }

    @Test
    public void canFindAffectedCodonNegativeStrand()
    {
        ExonData exon1 = createExon(3, 3, 5);
        ExonData exon2 = createExon(2, 11, 15);
        ExonData exon3 = createExon(1, 18, 20);
        TranscriptData transcript = createTranscript(Lists.newArrayList(exon1, exon2, exon3), 4, 19, Strand.REVERSE);

        assertNull(PaveAlgo.findAffectedCodon(transcript, 2));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 3));
        assertEquals(3, (int) PaveAlgo.findAffectedCodon(transcript, 4));
        assertEquals(3, (int) PaveAlgo.findAffectedCodon(transcript, 5));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 6));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 10));
        assertEquals(3, (int) PaveAlgo.findAffectedCodon(transcript, 11));
        assertEquals(2, (int) PaveAlgo.findAffectedCodon(transcript, 12));
        assertEquals(2, (int) PaveAlgo.findAffectedCodon(transcript, 13));
        assertEquals(2, (int) PaveAlgo.findAffectedCodon(transcript, 14));
        assertEquals(1, (int) PaveAlgo.findAffectedCodon(transcript, 15));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 16));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 17));
        assertEquals(1, (int) PaveAlgo.findAffectedCodon(transcript, 18));
        assertEquals(1, (int) PaveAlgo.findAffectedCodon(transcript, 19));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 20));
    }

    @Test
    public void canHandleNonCodingTranscripts()
    {
        ExonData exon1 = createExon(1, 3, 5);
        ExonData exon2 = createExon(2, 11, 15);
        ExonData exon3 = createExon(3, 18, 20);
        TranscriptData transcript = createTranscript(Lists.newArrayList(exon1, exon2, exon3), null, null, Strand.FORWARD);

        assertNull(PaveAlgo.findAffectedCodon(transcript, 4));
        assertNull(PaveAlgo.findAffectedCodon(transcript, 5));
    }

    @NotNull
    private static ExonData createExon(int rank, int start, int end)
    {
        return new ExonData(-1, start, end, rank, -1, -1);
    }

    @NotNull
    private static TranscriptData createTranscript(@NotNull List<ExonData> exons, @Nullable Integer codingStart,
            @Nullable Integer codingEnd, @NotNull Strand strand)
    {
        byte strandEnum = strand == Strand.FORWARD ? Strand.POS_STRAND : Strand.NEG_STRAND;
        TranscriptData transcript =
                new TranscriptData(-1, Strings.EMPTY, Strings.EMPTY, false, strandEnum, -1, -1, codingStart, codingEnd, Strings.EMPTY, null);
        transcript.setExons(exons);
        return transcript;
    }
}