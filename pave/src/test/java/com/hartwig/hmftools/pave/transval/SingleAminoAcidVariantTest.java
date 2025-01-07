package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

public class SingleAminoAcidVariantTest extends TransvalTestBase
{

    @Test
    public void referenceAminoAcidTest()
    {
        // Forward strand examples:
        Assert.assertEquals("M", this.variant("ADCK2:p.1E").referenceAminoAcid());
        Assert.assertEquals("V", this.variant("ADCK2:p.2E").referenceAminoAcid());
        Assert.assertEquals("R", this.variant("ADCK2:p.71E").referenceAminoAcid());
        Assert.assertEquals("K", this.variant("ADCK2:p.421E").referenceAminoAcid());
        Assert.assertEquals("P", this.variant("ADCK2:p.625E").referenceAminoAcid());
        Assert.assertEquals("P", this.variant("ADCK2:p.626E").referenceAminoAcid());

        // Reverse strand examples:
        Assert.assertEquals("M", this.variant("BRAF:p.1E").referenceAminoAcid());
        Assert.assertEquals("A", this.variant("BRAF:p.2E").referenceAminoAcid());
        Assert.assertEquals("V", this.variant("BRAF:p.600E").referenceAminoAcid());
        Assert.assertEquals("V", this.variant("BRAF:p.765E").referenceAminoAcid());
        Assert.assertEquals("H", this.variant("BRAF:p.766E").referenceAminoAcid());
    }

    @Test
    public void codingRegionLengthsTest()
    {
        // TODO cases where there are 1 or 2 exons
        SingleAminoAcidVariant variant = this.variant("ADCK2:p.3K");
        Assert.assertEquals(627, variant.AminoAcidSequence.AminoAcids.length()); // sanity
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(8, returned.size());
        Assert.assertEquals(3 * 627, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(933, returned.get(0).intValue());
        Assert.assertEquals(147, returned.get(1).intValue());
        Assert.assertEquals(54, returned.get(6).intValue());
        Assert.assertEquals(141, returned.get(7).intValue());
    }

    @Test
    public void codingRegionExons()
    {
        SingleAminoAcidVariant variant = this.variant("ZYX:p.1E");
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(9, returned.size());
        Assert.assertEquals(573, variant.AminoAcidSequence.AminoAcids.length()); // sanity
        Assert.assertEquals(3 * 573, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(208, returned.get(0).intValue()); // coding start > end of this exon
        Assert.assertEquals(200, returned.get(1).intValue());
        Assert.assertEquals(121, returned.get(7).intValue());
        Assert.assertEquals(105, returned.get(8).intValue());
    }

    @Test
    public void codingRegionLengthsReverseStrandTest()
    {
        SingleAminoAcidVariant variant = this.variant("BRAF:p.600E");
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(18, returned.size());
        Assert.assertEquals(3 * 767, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(138, returned.get(0).intValue());
        Assert.assertEquals(102, returned.get(1).intValue());
        Assert.assertEquals(135, returned.get(16).intValue());
        Assert.assertEquals(174, returned.get(17).intValue());
    }

    @Test
    public void variantAminoAcidTest()
    {
        Assert.assertEquals("E", this.variant("BRAF:p.1E").variantAminoAcid());
    }

    @Test
    public void referenceCodonTest()
    {
        Assert.assertEquals("ATG", variant("ADCK2:p.1E").referenceCodon(genome));
        Assert.assertEquals("GTG", variant("ADCK2:p.2E").referenceCodon(genome));
        Assert.assertEquals("CGC", variant("ADCK2:p.71E").referenceCodon(genome));
        Assert.assertEquals("TCA", variant("ADCK2:p.141E").referenceCodon(genome));
        Assert.assertEquals("GCT", variant("ADCK2:p.281E").referenceCodon(genome));
        Assert.assertEquals("GAG", variant("ADCK2:p.301E").referenceCodon(genome));
        Assert.assertEquals("AAA", variant("ADCK2:p.311E").referenceCodon(genome));
        Assert.assertEquals("GTG", variant("ADCK2:p.312E").referenceCodon(genome));
        Assert.assertEquals("TTG", variant("ADCK2:p.313E").referenceCodon(genome));
        Assert.assertEquals("E", variant("ADCK2:p.351E").referenceAminoAcid()); // sanity
        Assert.assertEquals("GAG", variant("ADCK2:p.351E").referenceCodon(genome));

        Assert.assertEquals("Q", variant("ADCK2:p.359E").referenceAminoAcid()); // sanity
        Assert.assertEquals("CAA", variant("ADCK2:p.359E").referenceCodon(genome));
        Assert.assertEquals("Q", variant("ADCK2:p.360E").referenceAminoAcid()); // sanity
        Assert.assertEquals("CAG", variant("ADCK2:p.360E").referenceCodon(genome));
        Assert.assertEquals("ATT", variant("ADCK2:p.361E").referenceCodon(genome));
        Assert.assertEquals("AAG", variant("ADCK2:p.580E").referenceCodon(genome));
        Assert.assertEquals("AAG", variant("ADCK2:p.580E").referenceCodon(genome));
        Assert.assertEquals("GTA", variant("ADCK2:p.581E").referenceCodon(genome));
        Assert.assertEquals("AAG", variant("ADCK2:p.582E").referenceCodon(genome));
        Assert.assertEquals("CCC", variant("ADCK2:p.625E").referenceCodon(genome));
        Assert.assertEquals("CCG", variant("ADCK2:p.626E").referenceCodon(genome));
        Assert.assertEquals("TGA", variant("ADCK2:p.627E").referenceCodon(genome)); // stop
    }

    @Test
    public void referenceCodonAcrossExonBoundariesTest()
    {
        Assert.assertEquals("ATG", variant("ZYX:p.1E").referenceCodon(genome));
        Assert.assertEquals("GCG", variant("ZYX:p.2E").referenceCodon(genome));
        Assert.assertEquals("CCG", variant("ZYX:p.68E").referenceCodon(genome));
        Assert.assertEquals("GAA", variant("ZYX:p.69E").referenceCodon(genome));
        Assert.assertEquals("GAC", variant("ZYX:p.70E").referenceCodon(genome)); // first nuke in exon 2, remainder in exon 3
        Assert.assertEquals("GTC", variant("ZYX:p.380E").referenceCodon(genome));
        Assert.assertEquals("AAC", variant("ZYX:p.381E").referenceCodon(genome));
        Assert.assertEquals("GAA", variant("ZYX:p.382E").referenceCodon(genome)); // first nuke in exon 4, remainder in exon 5
        Assert.assertEquals("TAC", variant("ZYX:p.496E").referenceCodon(genome));
        Assert.assertEquals("CAC", variant("ZYX:p.497E").referenceCodon(genome));
        Assert.assertEquals("AAG", variant("ZYX:p.498E").referenceCodon(genome)); // first 2 nukes in exon 8, last in exon 9
    }

    @Test
    public void referenceCodonReverseStrandGeneTest()
    {
        Assert.assertEquals("ATG", variant("BRAF:p.1E").referenceCodon(genome));
        Assert.assertEquals("GCG", variant("BRAF:p.2E").referenceCodon(genome));
        Assert.assertEquals("GCG", variant("BRAF:p.3E").referenceCodon(genome));
        Assert.assertEquals("CCG", variant("BRAF:p.44E").referenceCodon(genome));
        Assert.assertEquals("GAG", variant("BRAF:p.45E").referenceCodon(genome));
        Assert.assertEquals("GAG", variant("BRAF:p.46E").referenceCodon(genome)); // end of first exon
        Assert.assertEquals("GTG", variant("BRAF:p.47E").referenceCodon(genome));
        Assert.assertEquals("TGG", variant("BRAF:p.48E").referenceCodon(genome));
        Assert.assertEquals("CTG", variant("BRAF:p.79E").referenceCodon(genome));
        Assert.assertEquals("GAG", variant("BRAF:p.80E").referenceCodon(genome)); // end of second exon

        Assert.assertEquals("CAG", variant("BRAF:p.201E").referenceCodon(genome));
        Assert.assertEquals("GAT", variant("BRAF:p.202E").referenceCodon(genome));
        Assert.assertEquals("GGA", variant("BRAF:p.203E").referenceCodon(genome)); // 2 nukes in 4th exon, 1 in 5th

        Assert.assertEquals("TCT", variant("BRAF:p.325E").referenceCodon(genome));
        Assert.assertEquals("ATT", variant("BRAF:p.326E").referenceCodon(genome));
        Assert.assertEquals("GGG", variant("BRAF:p.327E").referenceCodon(genome)); // 2 nukes in 7th exon, 1 in 8th
        Assert.assertEquals("CCC", variant("BRAF:p.328E").referenceCodon(genome));

        Assert.assertEquals("TGG", variant("BRAF:p.476E").referenceCodon(genome));
        Assert.assertEquals("CAT", variant("BRAF:p.477E").referenceCodon(genome));
        Assert.assertEquals("GGT", variant("BRAF:p.478E").referenceCodon(genome)); // 1 nuke in an exon, 2 in the next
        Assert.assertEquals("GAT", variant("BRAF:p.479E").referenceCodon(genome));

        Assert.assertEquals("GTG", variant("BRAF:p.600E").referenceCodon(genome));
        Assert.assertEquals("AAA", variant("BRAF:p.601E").referenceCodon(genome));
        Assert.assertEquals("TCT", variant("BRAF:p.602E").referenceCodon(genome));
        Assert.assertEquals("GGA", variant("BRAF:p.670E").referenceCodon(genome));
    }

    @Test
    public void mtorTest()
    {
        // See TransvarConverterTest in the serve codebase
        Assert.assertEquals("TTA", variant("MTOR:p.L2230V").referenceCodon(genome));
    }

    @Test
    public void codonIsInSingleExon()
    {
        Assert.assertTrue(variant("ADCK2:p.1E").codonIsInSingleExon());
        Assert.assertTrue(variant("BRAF:p.46E").codonIsInSingleExon());
        Assert.assertFalse(variant("BRAF:p.327E").codonIsInSingleExon());
        Assert.assertFalse(variant("BRAF:p.478E").codonIsInSingleExon());
    }
}
