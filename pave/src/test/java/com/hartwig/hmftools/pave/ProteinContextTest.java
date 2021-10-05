package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class ProteinContextTest
{
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
