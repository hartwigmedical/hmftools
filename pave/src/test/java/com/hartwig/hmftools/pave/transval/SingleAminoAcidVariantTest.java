package com.hartwig.hmftools.pave.transval;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class SingleAminoAcidVariantTest extends TransvalTestBase
{

    @Test
    public void referenceValueTest()
    {
        // Forward strand examples:
        Assert.assertEquals("M", this.saav("ADCK2:p.1E").referenceAminoAcids());
        Assert.assertEquals("V", this.saav("ADCK2:p.2E").referenceAminoAcids());
        Assert.assertEquals("R", this.saav("ADCK2:p.71E").referenceAminoAcids());
        Assert.assertEquals("K", this.saav("ADCK2:p.421E").referenceAminoAcids());
        Assert.assertEquals("P", this.saav("ADCK2:p.625E").referenceAminoAcids());
        Assert.assertEquals("P", this.saav("ADCK2:p.626E").referenceAminoAcids());

        // Reverse strand examples:
        Assert.assertEquals("M", this.saav("BRAF:p.1E").referenceAminoAcids());
        Assert.assertEquals("A", this.saav("BRAF:p.2E").referenceAminoAcids());
        Assert.assertEquals("V", this.saav("BRAF:p.600E").referenceAminoAcids());
        Assert.assertEquals("V", this.saav("BRAF:p.765E").referenceAminoAcids());
        Assert.assertEquals("H", this.saav("BRAF:p.766E").referenceAminoAcids());
    }

    @Test
    public void codingRegionLengthsTest()
    {
        // TODO cases where there are 1 or 2 exons
        SingleAminoAcidVariant variant = this.saav("ADCK2:p.3K");
        Assert.assertEquals(627, variant.mAminoAcidSequence.AminoAcids.length()); // sanity
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
        SingleAminoAcidVariant variant = this.saav("ZYX:p.1E");
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(9, returned.size());
        Assert.assertEquals(573, variant.mAminoAcidSequence.AminoAcids.length()); // sanity
        Assert.assertEquals(3 * 573, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(208, returned.get(0).intValue()); // coding start > end of this exon
        Assert.assertEquals(200, returned.get(1).intValue());
        Assert.assertEquals(121, returned.get(7).intValue());
        Assert.assertEquals(105, returned.get(8).intValue());
    }

    @Test
    public void codingRegionLengthsReverseStrandTest()
    {
        SingleAminoAcidVariant variant = this.saav("BRAF:p.600E");
        List<Integer> returned = variant.codingRegionLengths();
        Assert.assertEquals(18, returned.size());
        Assert.assertEquals(3 * 767, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(138, returned.get(0).intValue());
        Assert.assertEquals(102, returned.get(1).intValue());
        Assert.assertEquals(135, returned.get(16).intValue());
        Assert.assertEquals(174, returned.get(17).intValue());
    }

    @Test
    public void altValueTest()
    {
        Assert.assertEquals("E", this.saav("BRAF:p.1E").altValue());
    }

    @Test
    public void referenceCodonTest()
    {
        Assert.assertEquals("ATG", saav("ADCK2:p.1E").referenceCodon(genome));
        Assert.assertEquals("GTG", saav("ADCK2:p.2E").referenceCodon(genome));
        Assert.assertEquals("CGC", saav("ADCK2:p.71E").referenceCodon(genome));
        Assert.assertEquals("TCA", saav("ADCK2:p.141E").referenceCodon(genome));
        Assert.assertEquals("GCT", saav("ADCK2:p.281E").referenceCodon(genome));
        Assert.assertEquals("GAG", saav("ADCK2:p.301H").referenceCodon(genome));
        Assert.assertEquals("AAA", saav("ADCK2:p.311E").referenceCodon(genome));
        Assert.assertEquals("GTG", saav("ADCK2:p.312E").referenceCodon(genome));
        Assert.assertEquals("TTG", saav("ADCK2:p.313E").referenceCodon(genome));
        Assert.assertEquals("E", saav("ADCK2:p.351H").referenceAminoAcids()); // sanity
        Assert.assertEquals("GAG", saav("ADCK2:p.351H").referenceCodon(genome));

        Assert.assertEquals("Q", saav("ADCK2:p.359E").referenceAminoAcids()); // sanity
        Assert.assertEquals("CAA", saav("ADCK2:p.359E").referenceCodon(genome));
        Assert.assertEquals("Q", saav("ADCK2:p.360E").referenceAminoAcids()); // sanity
        Assert.assertEquals("CAG", saav("ADCK2:p.360E").referenceCodon(genome));
        Assert.assertEquals("ATT", saav("ADCK2:p.361E").referenceCodon(genome));
        Assert.assertEquals("AAG", saav("ADCK2:p.580E").referenceCodon(genome));
        Assert.assertEquals("AAG", saav("ADCK2:p.580E").referenceCodon(genome));
        Assert.assertEquals("GTA", saav("ADCK2:p.581E").referenceCodon(genome));
        Assert.assertEquals("AAG", saav("ADCK2:p.582E").referenceCodon(genome));
        Assert.assertEquals("CCC", saav("ADCK2:p.625E").referenceCodon(genome));
        Assert.assertEquals("CCG", saav("ADCK2:p.626E").referenceCodon(genome));
        Assert.assertEquals("TGA", saav("ADCK2:p.627E").referenceCodon(genome)); // stop
    }

    @Test
    public void referenceCodonAcrossExonBoundariesTest()
    {
        Assert.assertEquals("ATG", saav("ZYX:p.1E").referenceCodon(genome));
        Assert.assertEquals("GCG", saav("ZYX:p.2E").referenceCodon(genome));
        Assert.assertEquals("CCG", saav("ZYX:p.68E").referenceCodon(genome));
        Assert.assertEquals("GAA", saav("ZYX:p.69W").referenceCodon(genome));
        Assert.assertEquals("GAC", saav("ZYX:p.70E").referenceCodon(genome)); // first nuke in exon 2, remainder in exon 3
        Assert.assertEquals("GTC", saav("ZYX:p.380E").referenceCodon(genome));
        Assert.assertEquals("AAC", saav("ZYX:p.381E").referenceCodon(genome));
        Assert.assertEquals("GAA", saav("ZYX:p.382W").referenceCodon(genome)); // first nuke in exon 4, remainder in exon 5
        Assert.assertEquals("TAC", saav("ZYX:p.496E").referenceCodon(genome));
        Assert.assertEquals("CAC", saav("ZYX:p.497E").referenceCodon(genome));
        Assert.assertEquals("AAG", saav("ZYX:p.498E").referenceCodon(genome)); // first 2 nukes in exon 8, last in exon 9
    }

    @Test
    public void referenceCodonReverseStrandGeneTest()
    {
        Assert.assertEquals("ATG", saav("BRAF:p.1E").referenceCodon(genome));
        Assert.assertEquals("GCG", saav("BRAF:p.2E").referenceCodon(genome));
        Assert.assertEquals("GCG", saav("BRAF:p.3E").referenceCodon(genome));
        Assert.assertEquals("CCG", saav("BRAF:p.44E").referenceCodon(genome));
        Assert.assertEquals("GAG", saav("BRAF:p.45W").referenceCodon(genome));
        Assert.assertEquals("GAG", saav("BRAF:p.46W").referenceCodon(genome)); // end of first exon
        Assert.assertEquals("GTG", saav("BRAF:p.47E").referenceCodon(genome));
        Assert.assertEquals("TGG", saav("BRAF:p.48E").referenceCodon(genome));
        Assert.assertEquals("CTG", saav("BRAF:p.79E").referenceCodon(genome));
        Assert.assertEquals("GAG", saav("BRAF:p.80Q").referenceCodon(genome)); // end of second exon

        Assert.assertEquals("CAG", saav("BRAF:p.201E").referenceCodon(genome));
        Assert.assertEquals("GAT", saav("BRAF:p.202E").referenceCodon(genome));
        Assert.assertEquals("GGA", saav("BRAF:p.203E").referenceCodon(genome)); // 2 nukes in 4th exon, 1 in 5th

        Assert.assertEquals("TCT", saav("BRAF:p.325E").referenceCodon(genome));
        Assert.assertEquals("ATT", saav("BRAF:p.326E").referenceCodon(genome));
        Assert.assertEquals("GGG", saav("BRAF:p.327E").referenceCodon(genome)); // 2 nukes in 7th exon, 1 in 8th
        Assert.assertEquals("CCC", saav("BRAF:p.328E").referenceCodon(genome));

        Assert.assertEquals("TGG", saav("BRAF:p.476E").referenceCodon(genome));
        Assert.assertEquals("CAT", saav("BRAF:p.477E").referenceCodon(genome));
        Assert.assertEquals("GGT", saav("BRAF:p.478E").referenceCodon(genome)); // 1 nuke in an exon, 2 in the next
        Assert.assertEquals("GAT", saav("BRAF:p.479E").referenceCodon(genome));

        Assert.assertEquals("GTG", saav("BRAF:p.600E").referenceCodon(genome));
        Assert.assertEquals("AAA", saav("BRAF:p.601E").referenceCodon(genome));
        Assert.assertEquals("TCT", saav("BRAF:p.602E").referenceCodon(genome));
        Assert.assertEquals("GGA", saav("BRAF:p.670E").referenceCodon(genome));
    }

    @Test
    public void mtorTest()
    {
        // See TransvarConverterTest in the serve codebase
        Assert.assertEquals("TTA", saav("MTOR:p.L2230V").referenceCodon(genome));
    }

    @Test
    public void codonIsInSingleExon()
    {
        Assert.assertTrue(saav("ADCK2:p.1E").codonIsInSingleExon());
        Assert.assertTrue(saav("BRAF:p.46A").codonIsInSingleExon());
        Assert.assertFalse(saav("BRAF:p.327E").codonIsInSingleExon());
        Assert.assertFalse(saav("BRAF:p.478E").codonIsInSingleExon());
    }

    @Test
    public void positionOfFirstAlteredCodonTest()
    {
        Assert.assertEquals(1, saav("ZYX:p.1E").positionOfFirstAlteredCodon());
        Assert.assertEquals(2230, saav("MTOR:p.L2230V").positionOfFirstAlteredCodon());
    }

    @Test
    public void positionOfLastAlteredCodonTest()
    {
        Assert.assertEquals(1, saav("ZYX:p.1E").positionOfLastAlteredCodon());
        Assert.assertEquals(2230, saav("MTOR:p.L2230V").positionOfLastAlteredCodon());
    }
}
