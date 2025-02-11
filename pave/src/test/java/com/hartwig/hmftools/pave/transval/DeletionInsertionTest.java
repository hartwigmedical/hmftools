package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

public class DeletionInsertionTest extends TransvalTestBase
{

    @Test
    public void referenceAminoAcidsTest()
    {
        assertEquals("MV", di("ADCK2:p.M1_V2delinsKQ").referenceAminoAcids());
        assertEquals("MVAP", di("ADCK2:p.M1_P4delinsK").referenceAminoAcids());
        assertEquals("EAT", di("ADCK2:p.E301_T303delinsQQ").referenceAminoAcids());
    }

    @Test
    public void altAminoAcidsTest()
    {
        assertEquals("KQ", di("ADCK2:p.M1_V2delinsKQ").altAminoAcidSequence());
        assertEquals("K", di("ADCK2:p.M1_P4delinsK").altAminoAcidSequence());
    }

    @Test
    public void positionTest()
    {
        assertEquals(1, di("ADCK2:p.M1_V2delinsKQ").positionOfFirstAlteredCodon());
        assertEquals(2, di("ADCK2:p.M1_V2delinsKQ").positionOfLastAlteredCodon());
        assertEquals(111, di("VHL:p.S111_L116delinsKQ").positionOfFirstAlteredCodon());
        assertEquals(116, di("VHL:p.S111_L116delinsKQ").positionOfLastAlteredCodon());
    }

    @Test
    public void referenceBasesTest()
    {
        check(seq("ATGCCC", ""), di("VHL:p.M1_P2delinsKQ").referenceBases(genome));
        check(seq("ATGCCC", ""), di("VHL:p.M1_P2delinsKQ").referenceBases(genome));
        check(seq("ATGCCCCGG", ""), di("VHL:p.M1_R3delinsKQ").referenceBases(genome));
        check(seq("CCCCGG", ""), di("VHL:p.P2_R3delinsKQ").referenceBases(genome));
        check(seq("GAGAACTGGGAC", ""), di("VHL:p.E6_D9delinsKQ").referenceBases(genome));

        check(seq("AGCTACCGA", ""), di("VHL:p.S111_R113delinsKQ").referenceBases(genome));
        check(seq("AGCTACCGAG", "GT"), di("VHL:p.S111_G114delinsKQ").referenceBases(genome));

        check(seq("CACTGTGTCCCCGACTAC", ""), di("ZYX:p.H491_Y496delinsW").referenceBases(genome));
        check(seq("TGTGTCCCCGACTACCAC", ""), di("ZYX:p.C492_H497delinsW").referenceBases(genome));
        check(seq("GTCCCCGACTACCACAA", "G"), di("ZYX:p.V493_K498delinsW").referenceBases(genome));
    }

    private void check(SplitCodonSequence expected, SplitCodonSequence actual)
    {
        assertEquals(expected.completeSequence(), actual.completeSequence());
        assertEquals(expected.retainedPrefix(), actual.retainedPrefix());
        assertEquals(expected.retainedSuffix(), actual.retainedSuffix());
    }

    @Test
    public void candidateAlternativeNucleotideSequencesTest()
    {
        checkAltSequences(di("VHL:p.M1_P2delinsW"), Collections.singleton("TGG"));

        checkAltSequences(di("VHL:p.M1_P2delinsWW"), Collections.singleton("TGGTGG"));

        Set<String> expectedAW = new HashSet<>();
        expectedAW.add("GCTTGG");
        expectedAW.add("GCCTGG");
        expectedAW.add("GCATGG");
        expectedAW.add("GCGTGG");
        checkAltSequences(di("VHL:p.M1_P2delinsAW"), expectedAW);

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

    @Test
    public void candidateAlternativeNucleotideSequencesWithFixedPrefixTest()
    {
        checkAltSequences(di("VHL:p.M1_P2delinsW"), Collections.singleton("TGG"));

        checkAltSequences(di("VHL:p.M1_P2delinsWW"), Collections.singleton("TGGTGG"));

        Set<String> expectedAW = new HashSet<>();
        expectedAW.add("GCTTGG");
        expectedAW.add("GCCTGG");
        expectedAW.add("GCATGG");
        expectedAW.add("GCGTGG");
        checkAltSequences(di("VHL:p.M1_P2delinsAW"), expectedAW);

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

    @Test
    public void egfrTest()
    {
        // This test arose from a bug that lead to the incorrect handling of EGFR:p.L747_K754delinsSPQ.
        Set<String> candidates = di("EGFR:p.L747_K754delinsSPQ").candidateAlternativeNucleotideSequences("", "");
        AminoAcidSequence spq = aaSeq("SPQ");
        candidates.forEach(nucleotides -> assertEquals(spq, AminoAcidSequence.fromNucleotides(nucleotides)));
        assertEquals(6 * 4 * 2, candidates.size());
    }

    private void checkAltSequences(DeletionInsertion di, Set<String> expected)
    {
        Set<String> actual = di.candidateAlternativeNucleotideSequences("", "");
        assertEquals(expected, actual);
    }

    private DeletionInsertion di(String definition)
    {
        return transval.variationParser().parseDeletionInsertion(definition);
    }
}
