package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Test;

public class DeletionInsertionTest extends TransvalTestBase
{

    @Test
    public void referenceAminoAcidsTest()
    {
        assertEquals("M", di("ADCK2:p.M1_M1delinsKQ").referenceAminoAcids());
        assertEquals("MVAP", di("ADCK2:p.M1_P4delinsK").referenceAminoAcids());
        assertEquals("EAT", di("ADCK2:p.E301_T303delinsQQ").referenceAminoAcids());
        assertEquals("PX", di("ADCK2:p.E626_T627delinsQQ").referenceAminoAcids());
    }

    @Test
    public void altAminoAcidsTest()
    {
        assertEquals("KQ", di("ADCK2:p.M1_M1delinsKQ").altAminoAcidSequence());
        assertEquals("K", di("ADCK2:p.M1_P4delinsK").altAminoAcidSequence());
    }

    @Test
    public void positionTest()
    {
        assertEquals(1, di("ADCK2:p.M1_M1delinsKQ").positionOfFirstAlteredCodon());
        assertEquals(1, di("ADCK2:p.M1_M1delinsKQ").positionOfLastAlteredCodon());
        assertEquals(111, di("VHL:p.S111_L116delinsKQ").positionOfFirstAlteredCodon());
        assertEquals(116, di("VHL:p.S111_L116delinsKQ").positionOfLastAlteredCodon());
    }

    @Test
    public void referenceBasesTest()
    {
        assertEquals(seq("ATG", null), di("VHL:p.M1_M1delinsKQ").referenceBases(genome));
        assertEquals(seq("ATGCCC", null), di("VHL:p.M1_P2delinsKQ").referenceBases(genome));
        assertEquals(seq("ATGCCCCGG", null), di("VHL:p.M1_R3delinsKQ").referenceBases(genome));
        assertEquals(seq("CCCCGG", null), di("VHL:p.P2_R3delinsKQ").referenceBases(genome));
        assertEquals(seq("GAGAACTGGGAC", null), di("VHL:p.E6_D9delinsKQ").referenceBases(genome));

        assertEquals(seq("AGCTACCGA", null), di("VHL:p.S111_R113delinsKQ").referenceBases(genome));
        assertEquals(seq("AGCTACCGAG", "GT"), di("VHL:p.S111_G114delinsKQ").referenceBases(genome));
        assertEquals(seq("AGCTACCGAG", "GTCAC"), di("VHL:p.S111_H115delinsKQ").referenceBases(genome));
        assertEquals(seq("AGCTACCGAG", "GTCACCTT"), di("VHL:p.S111_L116delinsKQ").referenceBases(genome));
        assertEquals(seq("TACCGAG", "GTCACCTT"), di("VHL:p.Y112_L116delinsKQ").referenceBases(genome));
        assertEquals(seq("CGAG", "GTCACCTT"), di("VHL:p.R113_L116delinsKQ").referenceBases(genome));
        assertEquals(seq("G", "GTCACCTT"), di("VHL:p.G114_L116delinsKQ").referenceBases(genome));

        assertEquals(seq("CACTGTGTCCCCGACTAC", null), di("ZYX:p.H491_Y496delinsW").referenceBases(genome));
        assertEquals(seq("TGTGTCCCCGACTACCAC", null), di("ZYX:p.C492_H497delinsW").referenceBases(genome));
        assertEquals(seq("GTCCCCGACTACCACAA", "G"), di("ZYX:p.V493_K498delinsW").referenceBases(genome));
        assertEquals(seq("CCCGACTACCACAA", "GCAG"), di("ZYX:p.P494_Q499delinsW").referenceBases(genome));
        assertEquals(seq("GACTACCACAA", "GCAGTAC"), di("ZYX:p.D495_Y500delinsW").referenceBases(genome));
    }

    @Test
    public void candidateAlternativeCodonsTest()
    {
        checkAltCodons(di("VHL:p.M1_M1delinsW"), "TGG");
        checkAltCodons(di("VHL:p.M1_M1delinsWW"), "TGG:TGG");
        checkAltCodons(di("VHL:p.M1_M1delinsAW"), "GCT,GCC,GCA,GCG:TGG");
        checkAltCodons(di("VHL:p.M1_M1delinsWA"), "TGG:GCT,GCC,GCA,GCG");
        checkAltCodons(di("VHL:p.M1_M1delinsAA"), "GCT,GCC,GCA,GCG:GCT,GCC,GCA,GCG");
        checkAltCodons(di("VHL:p.R113_L116delinsKP"), "AAA,AAG:CCT,CCC,CCA,CCG");
        checkAltCodons(di("VHL:p.R113_L116delinsKPY"), "AAA,AAG:CCT,CCC,CCA,CCG:TAT,TAC");
    }

    @Test
    public void candidateAlternativeNucleotideSequencesTest()
    {
        checkAltSequences(di("VHL:p.M1_M1delinsW"), Collections.singleton("TGG"));

        checkAltSequences(di("VHL:p.M1_M1delinsWW"), Collections.singleton("TGGTGG"));

        Set<String> expectedAW = new HashSet<>();
        expectedAW.add("GCTTGG");
        expectedAW.add("GCCTGG");
        expectedAW.add("GCATGG");
        expectedAW.add("GCGTGG");
        checkAltSequences(di("VHL:p.M1_M1delinsAW"), expectedAW);

        Set<String> expectedKPY = new HashSet<>();
        expectedKPY.add("AAACCTTAT");
        expectedKPY.add("AAACCTTAC");
        expectedKPY.add("AAACCCTAT");
        expectedKPY.add("AAACCCTAC");
        expectedKPY.add("AAACCATAT");
        expectedKPY.add("AAACCATAC");
        expectedKPY.add("AAACCGTAT");
        expectedKPY.add("AAACCGTAC");
        expectedKPY.add("AAGCCTTAT");
        expectedKPY.add("AAGCCTTAC");
        expectedKPY.add("AAGCCCTAT");
        expectedKPY.add("AAGCCCTAC");
        expectedKPY.add("AAGCCATAT");
        expectedKPY.add("AAGCCATAC");
        expectedKPY.add("AAGCCGTAT");
        expectedKPY.add("AAGCCGTAC");
        checkAltSequences(di("VHL:p.R113_L116delinsKPY"), expectedKPY);
    }

    private void checkAltSequences(DeletionInsertion di, Set<String> expected)
    {
        Set<String> actual = di.candidateAlternativeNucleotideSequences();
        assertEquals(expected, actual);
    }

    private void checkAltCodons(DeletionInsertion di, String expectedAltCodonsSeparatedByColonsAndCommas)
    {
        List<Set<String>> actual = di.candidateAlternativeCodons();
        String[] list = expectedAltCodonsSeparatedByColonsAndCommas.split(":");
        assertEquals(actual.size(), list.length);
        for(int i = 0; i < actual.size(); ++i) {
            Set<String> expected = Arrays.stream(list[i].split(",")).collect(Collectors.toSet());
            assertEquals(expected, actual.get(i));
        }
    }

    private DeletionInsertion di(String definition)
    {
        return transval.variationParser().parseDeletionInsertion(definition);
    }

}
