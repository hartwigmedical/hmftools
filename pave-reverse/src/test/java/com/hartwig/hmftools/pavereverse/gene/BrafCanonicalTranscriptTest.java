package com.hartwig.hmftools.pavereverse.gene;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.junit.Before;
import org.junit.Test;

public class BrafCanonicalTranscriptTest extends ReversePaveTestBase
{
    private GeneTranscript geneTranscript;

    @Before
    public void setUp()
    {
        GeneData gene = reversePave.EnsemblCache.getGeneDataByName(braf);
        TranscriptData transcriptData = reversePave.EnsemblCache.getTranscriptData(gene.GeneId, brafCanonical);
        geneTranscript = new GeneTranscript(gene, transcriptData);
    }

    @Test
    public void totalTranslatedLength()
    {
        assertEquals(767*3, geneTranscript.totalTranslatedLength());
    }

    @Test
    public void codingRegionLengths()
    {
        List<Integer> lengths = geneTranscript.codingRegionLengths();
        assertEquals(18, lengths.size());
        assertEquals(138, (int) lengths.get(0));
        assertEquals(102, (int) lengths.get(1));
        assertEquals(264, (int) lengths.get(2));
    }

    @Test
    public void absolutePositionOfTranslatedBase()
    {
        assertEquals(140_924_703, geneTranscript.absolutePositionOfTranslatedBase(1));
        assertEquals(140_924_702, geneTranscript.absolutePositionOfTranslatedBase(2));
        assertEquals(140_924_566, geneTranscript.absolutePositionOfTranslatedBase(138));

        assertEquals(140_850_212, geneTranscript.absolutePositionOfTranslatedBase(139));
        assertEquals(140_850_211, geneTranscript.absolutePositionOfTranslatedBase(140));
        assertEquals(140_850_111, geneTranscript.absolutePositionOfTranslatedBase(240));

        assertEquals(140_834_872, geneTranscript.absolutePositionOfTranslatedBase(241));
        assertEquals(140_834_609, geneTranscript.absolutePositionOfTranslatedBase(504));

        assertEquals(140_734_597, geneTranscript.absolutePositionOfTranslatedBase(767*3));
    }

    @Test
    public void absolutePositionOf5PrimeUtrExonicBase()
    {
        assertEquals(140_924_704, geneTranscript.absolutePositionOf5PrimeUtrExonicBase(-1));
    }

    @Test
    public void absolutePositionOf3PrimeUtrExonicBase()
    {
        assertEquals(140_734_596, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(1));
        assertEquals(140_734_595, geneTranscript.absolutePositionOf3PrimeUtrExonicBase(2));
    }
}
