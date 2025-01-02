package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

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
    public void exonCodingLengthsTest()
    {
        // TODO cases where there are 1 or 2 exons
        SingleAminoAcidVariant variant = this.variant("ADCK2:p.3K");
        Assert.assertEquals(627, variant.aminoAcidSequence.AminoAcids.length()); // sanity
        List<Integer> returned = variant.exonCodingLengths();
        Assert.assertEquals(8, returned.size());
        Assert.assertEquals(3 * 627, returned.stream().mapToInt(Integer::intValue).sum());
        Assert.assertEquals(933, returned.get(0).intValue());
        Assert.assertEquals(147, returned.get(1).intValue());
        Assert.assertEquals(54, returned.get(6).intValue());
        Assert.assertEquals(141, returned.get(7).intValue());
    }

    @Test
    public void exonCodingLengthsReverseStrandTest()
    {
        SingleAminoAcidVariant variant = this.variant("BRAF:p.600E");
        List<Integer> returned = variant.exonCodingLengths();
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
    public void possibleVariantCodonsTest()
    {
        Assert.assertEquals(Set.of("TTC", "TTT"), this.variant("BRAF:p.1F").possibleVariantCodons());
        Assert.assertEquals(Set.of("GGA", "GGC", "GGG", "GGT"), this.variant("BRAF:p.7G").possibleVariantCodons());
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

    }

    @Test
    public void referenceCodonReverseStrandGeneTest() throws FileNotFoundException
    {

        RefGenomeSource refGenomeSource = new RefGenomeSource(new IndexedFastaSequenceFile(new File("/Users/timlavers/work/data/reference_genome_no_alts/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")));

//        Assert.assertEquals("GTG", variant("BRAF:p.600E").referenceCodon(genome));
//        Assert.assertEquals("AAA", variant("BRAF:p.601E").referenceCodon(genome));
        Assert.assertEquals("TCT", variant("BRAF:p.602E").referenceCodon(refGenomeSource));
//        Assert.assertEquals("GGA", variant("BRAF:p.670E").referenceCodon(genome));
    }
}
