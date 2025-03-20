package com.hartwig.hmftools.pavereverse.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.junit.Test;

public class GeneTranscriptTest extends ReversePaveTestBase
{
    int[] exonStarts = { 10, 30, 50, 70, 90, 110 };
    int codingStart = 55; // inclusive
    int codingEnd = 94; // inclusive
    // 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    // 1        10        20        30        40        50        60        70        80        90        100       110       120  125
    // ____-----++++++++++----------++++++++++----------+++++*****----------**********----------*****+++++----------++++++++++-----_____

    GeneData geneData = createEnsemblGeneData("id132", "BLAH", "1",  POS_STRAND, 5, 125);
    TranscriptData transcript = createTransExons(geneData.GeneId, 123, POS_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
    TranscriptData rsTranscript = createTransExons(geneData.GeneId, 124, NEG_STRAND, exonStarts, 9, codingStart, codingEnd, false, "whatever");
    GeneTranscript geneTranscript = new GeneTranscript(geneData, transcript);
    GeneTranscript geneTranscriptRS = new GeneTranscript(geneData, rsTranscript);

    @Test
    public void totalTranslatedLength()
    {
        assertEquals(20, geneTranscript.totalTranslatedLength());
        assertEquals(20, geneTranscriptRS.totalTranslatedLength());
    }

    @Test
    public void codingRegionLengths()
    {
        assertEquals(List.of(5, 10, 5), geneTranscript.codingRegionLengths());
        assertEquals(List.of(5, 10, 5), geneTranscriptRS.codingRegionLengths());
    }

    @Test
    public void absolutePositionOfTranslatedBase()
    {
        assertEquals(55, geneTranscript.absolutePositionOfTranslatedBase(1));
        assertEquals(56, geneTranscript.absolutePositionOfTranslatedBase(2));
        assertEquals(59, geneTranscript.absolutePositionOfTranslatedBase(5));
        assertEquals(70, geneTranscript.absolutePositionOfTranslatedBase(6));
        assertEquals(71, geneTranscript.absolutePositionOfTranslatedBase(7));
        assertEquals(79, geneTranscript.absolutePositionOfTranslatedBase(15));
        assertEquals(90, geneTranscript.absolutePositionOfTranslatedBase(16));
        assertEquals(94, geneTranscript.absolutePositionOfTranslatedBase(20));
    }

    @Test
    public void absolutePositionOfTranslatedBaseReverseStrand()
    {
        assertEquals(94, geneTranscriptRS.absolutePositionOfTranslatedBase(1));
        assertEquals(93, geneTranscriptRS.absolutePositionOfTranslatedBase(2));
        assertEquals(90, geneTranscriptRS.absolutePositionOfTranslatedBase(5));
        assertEquals(79, geneTranscriptRS.absolutePositionOfTranslatedBase(6));
        assertEquals(78, geneTranscriptRS.absolutePositionOfTranslatedBase(7));
        assertEquals(70, geneTranscriptRS.absolutePositionOfTranslatedBase(15));
        assertEquals(59, geneTranscriptRS.absolutePositionOfTranslatedBase(16));
        assertEquals(55, geneTranscriptRS.absolutePositionOfTranslatedBase(20));
    }

    @Test
    public void absolutePositionOf5PrimeUtrExonicBase()
    {
        assertEquals(54, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-1));
        assertEquals(53, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-2));
        assertEquals(50, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-5));
        assertEquals(39, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-6));
        assertEquals(30, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-15));
        assertEquals(19, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-16));
        assertEquals(10, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-25));
    }

    @Test
    public void absolutePositionOf5PrimeUtrExonicBaseReverseStrand()
    {
        assertEquals(95, geneTranscriptRS.absolutePositionOf5PrimeUtrExonicBase(-1));
        assertEquals(96, geneTranscriptRS.absolutePositionOf5PrimeUtrExonicBase(-2));
        assertEquals(99, geneTranscriptRS.absolutePositionOf5PrimeUtrExonicBase(-5));
        assertEquals(110, geneTranscriptRS.absolutePositionOf5PrimeUtrExonicBase(-6));
        assertEquals(119, geneTranscriptRS.absolutePositionOf5PrimeUtrExonicBase(-15));
    }

    @Test
    public void absolutePositionOf3PrimeUtrExonicBase()
    {
        assertEquals(95, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(1));
        assertEquals(96, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(2));
        assertEquals(99, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(5));
        assertEquals(110, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(6));
        assertEquals(119, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(15));
    }

    @Test
    public void absolutePositionOf3PrimeUtrExonicBaseReverseStrand()
    {
        assertEquals(54, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(1));
        assertEquals(53, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(2));
        assertEquals(50, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(5));
        assertEquals(39, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(6));
        assertEquals(30, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(15));
        assertEquals(19, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(16));
        assertEquals(10, geneTranscriptRS.absolutePositionOf3PrimeUtrExonicBase(25));
    }
}
