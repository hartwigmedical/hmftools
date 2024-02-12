package com.hartwig.hmftools.neo.missense;

import static com.hartwig.hmftools.common.codon.AminoAcids.AMINO_ACID_TO_CODON_MAP;
import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.codon.Codons.STOP_CODON_1;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class MissensePeptideTest
{
    @Test
    public void testMissensePepetides()
    {
        MissenseConfig config = new MissenseConfig(3, 0);
        MockRefGenome refGenome = new MockRefGenome();

        String aminoAcids = "MACDEFGHIK";
        String codingBases = getCodingBases(aminoAcids) + STOP_CODON_1;
        String revCodingBases = Nucleotides.reverseComplementBases(codingBases);
        String refBases = generateRandomBases(16) + revCodingBases.substring(0, 5) + generateRandomBases(9) +
                revCodingBases.substring(5, 16) + generateRandomBases(9) + revCodingBases.substring(16, 27)
                + generateRandomBases(9) + revCodingBases.substring(27, 33) + generateRandomBases(30);

        refGenome.RefGenomeMap.put(CHR_1, refBases);
        MissenseCalcs calcs = new MissenseCalcs(config, refGenome);

        GeneData testGene = new GeneData(GENE_ID_1, GENE_NAME_1, CHR_1, NEG_STRAND, 10, 80, "");

        int codingStart = 16;
        int codingEnd = 75;
        TranscriptData transDataNeg = GeneTestUtils.createTransExons(
                testGene.GeneId, TRANS_ID_1, NEG_STRAND,
                new int[] {10, 30, 50, 70}, 10, codingStart, codingEnd, true, "");

        calcs.processTranscript(testGene, transDataNeg);

        assertEquals(178, calcs.peptideData().size());
        MissensePeptide first = calcs.peptideData().get(0);
        assertEquals(1, first.CodonIndex);
        assertEquals(75, first.Position);
        assertEquals("CAT", first.Context);
        assertEquals('T', first.RefBase);
        assertEquals('G', first.AltBase);
        assertEquals("LAC", first.Peptide);

        MissensePeptide last = calcs.peptideData().get(calcs.peptideData().size() - 1);
        assertEquals(11, last.CodonIndex);
        assertEquals(16, last.Position);
        assertEquals("TTA", last.Context);
        assertEquals('T', last.RefBase);
        assertEquals('A', last.AltBase);
        assertEquals("IKY", last.Peptide);

        for(MissensePeptide missensePeptide : calcs.peptideData())
        {
            assertEquals(1, mismatchCount(aminoAcids + STOP_AMINO_ACID, missensePeptide));
        }

        // now with flanks
        config = new MissenseConfig(3, 3);

        calcs = new MissenseCalcs(config, refGenome);
        calcs.processTranscript(testGene, transDataNeg);

        assertEquals(178, calcs.peptideData().size());

        for(MissensePeptide missensePeptide : calcs.peptideData())
        {
            assertEquals(1, mismatchCount(aminoAcids + STOP_AMINO_ACID, missensePeptide));
        }
    }

    private static int mismatchCount(final String allAminoAcids, final MissensePeptide missensePeptide)
    {
        String fullSequence = missensePeptide.UpFlank + missensePeptide.Peptide + missensePeptide.DownFlank;

        int mismatchCount = 0;
        int aaOffset = missensePeptide.PeptideStartIndex - 1 - missensePeptide.UpFlank.length();

        for(int pepIndex = 0; pepIndex < fullSequence.length(); ++pepIndex)
        {
            int aaIndex = pepIndex + aaOffset;
            if(allAminoAcids.charAt(aaIndex) != fullSequence.charAt(pepIndex))
                ++mismatchCount;
        }

        return mismatchCount;
    }

    private static String getCodingBases(final String aminoAcids)
    {
        String bases = "";

        for(int i = 0; i < aminoAcids.length(); ++i)
        {
            bases += AMINO_ACID_TO_CODON_MAP.get(aminoAcids.substring(i, i + 1)).get(0);
        }

        return bases;
    }
}
