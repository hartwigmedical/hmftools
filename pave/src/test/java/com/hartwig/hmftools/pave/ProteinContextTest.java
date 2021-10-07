package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.pave.CodingUtils.getExtraBases;

import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class ProteinContextTest
{
    @Test
    public void testExtraCodingBases()
    {
        final MockRefGenome refGenome = new MockRefGenome();

        String chr1Bases = generateRandomBases(100);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        int[] exonStarts = { 10, 30, 50, 70};

        // codons start on at 10, 13, 16 etc
        Integer codingStart = new Integer(15);
        Integer codingEnd = new Integer(75);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        ExonData exon = transDataPosStrand.exons().get(1);

        String bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 32, 1, false);
        String refBases = chr1Bases.substring(31, 32);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 32, 2, false);
        refBases = chr1Bases.substring(30, 32);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 32, 3, false);
        refBases = chr1Bases.substring(20, 21) + chr1Bases.substring(30, 32);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 38, 4, true);
        refBases = chr1Bases.substring(39, 41) + chr1Bases.substring(50, 52);
        assertTrue(refBases.equals(bases));

        TranscriptData transDataNegStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, NEG_STRAND, exonStarts, 10, codingStart, codingEnd, false, "");

        exon = transDataPosStrand.exons().get(2);

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 52, 4, false);
        refBases = chr1Bases.substring(39, 41) + chr1Bases.substring(50, 52);
        assertTrue(refBases.equals(bases));

        bases = getExtraBases(transDataPosStrand, refGenome, CHR_1, exon, 58, 6, true);
        refBases = chr1Bases.substring(59, 61) + chr1Bases.substring(70, 74);
        assertTrue(refBases.equals(bases));

    }

    @Test
    public void testBaseMutations()
    {
        // SNVs and MNVs
        final MockRefGenome refGenome = new MockRefGenome();

        int[] exonStarts = { 0, 100, 200 };

        // codons start on at 10, 13, 16 etc
        Integer codingStart = new Integer(10);
        Integer codingEnd = new Integer(250);

        TranscriptData transDataPosStrand = createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 80, codingStart, codingEnd, false, "");

        String chr1Bases = generateRandomBases(300);

        // set the specific AAs for a region of this mock ref genome
        // S: TCA -> TCG - so last base of codon can change
        // I: ATC -> ATT
        // L: TTA -> CTA - first base changes

        //                     40      43      46
        String aminoAcidSeq = "TCA" + "ATC" + "TTA";
        chr1Bases = chr1Bases.substring(0, 40) + aminoAcidSeq + chr1Bases.substring(49);
        refGenome.RefGenomeMap.put(CHR_1, chr1Bases);

        // test SNVs at the 1st and 3rd codon bases

        // last base of a codon changes
        int codonPos = 40;
        String codon = chr1Bases.substring(codonPos, codonPos + 3);
        String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);

        String alt = "G";
        String synCodon = chr1Bases.substring(codonPos, codonPos + 2) + alt;
        assertTrue(aminoAcid.equals(AminoAcids.findAminoAcidForCodon(synCodon)));

        int pos = 42;
        String ref = chr1Bases.substring(pos, pos + 1);
        VariantData var = new VariantData(CHR_1, pos, ref, alt);

        /*
        ProteinContext proteinContext = ProteinContext.determineContext(var, transDataPosStrand, refGenome);
        assertTrue(proteinContext.hasCodingBases());
        assertEquals("S", proteinContext.WildtypeAA);
        assertEquals("S", proteinContext.NovelAA);
         */
    }

}
