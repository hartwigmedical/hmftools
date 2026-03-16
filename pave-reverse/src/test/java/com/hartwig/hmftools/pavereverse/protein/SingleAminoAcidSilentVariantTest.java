package com.hartwig.hmftools.pavereverse.protein;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.pavereverse.base.BaseSequence;
import com.hartwig.hmftools.pavereverse.base.CodonWithinExons;

import org.junit.Test;

public class SingleAminoAcidSilentVariantTest extends VariantTest
{
    @Test
    public void variantSequenceTest()
    {
        SingleAminoAcidSilentVariant variant = new SingleAminoAcidSilentVariant(gene, transcript, taa, 5);
        assertEquals("MAAQVAPAAS", variant.variantSequence().sequence());
    }

    @Test
    public void possibleVariantsTest()
    {
        SingleAminoAcidSilentVariant variant = new SingleAminoAcidSilentVariant(gene, transcript, taa, 4);
        BaseSequence bases = new BaseSequence(0, "CAG", true);
        CodonWithinExons codon = CodonWithinExons.factory("", bases, "");
        Set<CodonChange> changes = variant.possibleVariants(codon);
        assertEquals(1, changes.size());
        assertEquals("CAA", changes.iterator().next().AlternateCodon);
        assertEquals("CAG", changes.iterator().next().ReferenceCodon);

        bases = new BaseSequence(0, "CAA", true);
        codon = CodonWithinExons.factory("", bases, "");
        changes = variant.possibleVariants(codon);
        assertEquals(1, changes.size());
        assertEquals("CAG", changes.iterator().next().AlternateCodon);
        assertEquals("CAA", changes.iterator().next().ReferenceCodon);

        variant = new SingleAminoAcidSilentVariant(gene, transcript, taa, 5);
        bases = new BaseSequence(0, "GTC", true);
        codon = CodonWithinExons.factory("", bases, "");
        changes = variant.possibleVariants(codon);
        assertEquals(3, changes.size());
        for(CodonChange change : changes)
        {
            assertEquals("GTC", change.ReferenceCodon);
            assertNotEquals("GTC", change.AlternateCodon);
            assertEquals("V", AminoAcids.findAminoAcidForCodon(change.AlternateCodon));
        }
    }

    @Test
    public void codingRegionLengthsTest()
    {
        SingleAminoAcidSilentVariant variant = this.saasv("ADCK2:p.V2V");
        assertEquals(627, variant.AminoAcidsTranscript.AminoAcids.length()); // sanity
        List<Integer> returned = variant.codingRegionLengths();
        assertEquals(8, returned.size());
        assertEquals(3 * 627, returned.stream().mapToInt(Integer::intValue).sum());
        assertEquals(933, returned.get(0).intValue());
        assertEquals(147, returned.get(1).intValue());
        assertEquals(54, returned.get(6).intValue());
        assertEquals(141, returned.get(7).intValue());
    }

    @Test
    public void codingRegionExons()
    {
        SingleAminoAcidSilentVariant variant = this.saasv("ZYX:p.A2A");
        List<Integer> returned = variant.codingRegionLengths();
        assertEquals(9, returned.size());
        assertEquals(573, variant.AminoAcidsTranscript.AminoAcids.length()); // sanity
        assertEquals(3 * 573, returned.stream().mapToInt(Integer::intValue).sum());
        assertEquals(208, returned.get(0).intValue()); // coding start > end of this exon
        assertEquals(200, returned.get(1).intValue());
        assertEquals(121, returned.get(7).intValue());
        assertEquals(105, returned.get(8).intValue());
    }

    @Test
    public void codingRegionLengthsReverseStrandTest()
    {
        SingleAminoAcidSilentVariant variant = this.saasv("BRAF:p.V600V");
        List<Integer> returned = variant.codingRegionLengths();
        assertEquals(18, returned.size());
        assertEquals(3 * 767, returned.stream().mapToInt(Integer::intValue).sum());
        assertEquals(138, returned.get(0).intValue());
        assertEquals(102, returned.get(1).intValue());
        assertEquals(135, returned.get(16).intValue());
        assertEquals(174, returned.get(17).intValue());
    }

    @Test
    public void altValueTest()
    {
        assertEquals("P", this.saasv("VHL:p.P2P").altValue());
    }

    @Test
    public void positionOfFirstAlteredCodonTest()
    {
        assertEquals(2, saasv("VHL:p.P2P").positionOfFirstAlteredCodon());
        assertEquals(2230, saasv("MTOR:p.L2230L").positionOfFirstAlteredCodon());
    }
}
