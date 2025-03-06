package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Assert;
import org.junit.Test;

public class NucleotidesCalculatorTest extends ReversePaveTestBase
{
    @Test
    public void someBaseSequenceTest()
    {
        Assert.assertEquals("TGG", calculator("W", "", "").anyBaseSequence());
        Assert.assertEquals("TGGTGG", calculator("WW", "", "").anyBaseSequence());
        String s = calculator("FEAST", "", "").anyBaseSequence();
        Assert.assertEquals("FEAST", AminoAcidSequence.fromNucleotides(s).sequence());
    }

    @Test
    public void oneAA()
    {
        checkPossibilities(calculator("W", "C", ""), "");
        checkPossibilities(calculator("W", "", "C"), "");
        checkPossibilities(calculator("W", "T", "C"), "");

        checkPossibilities(calculator("W", "", ""), "TGG");
        checkPossibilities(calculator("W", "", "G"), "TG");
        checkPossibilities(calculator("W", "T", ""), "GG");
        checkPossibilities(calculator("W", "TG", ""), "G");
        checkPossibilities(calculator("W", "T", "G"), "G");

        checkPossibilities(calculator("Y", "", ""), "TAT,TAC");
        checkPossibilities(calculator("Y", "T", ""), "AT,AC");
        checkPossibilities(calculator("Y", "TA", ""), "T,C");
        checkPossibilities(calculator("Y", "", "C"), "TA");
        checkPossibilities(calculator("Y", "", "AC"), "T");
        checkPossibilities(calculator("Y", "T", "C"), "A");
    }

    @Test
    public void twoAAs()
    {
        checkPossibilities(calculator("WW", "C", ""), "");
        checkPossibilities(calculator("WW", "", "C"), "");
        checkPossibilities(calculator("WW", "", ""), "TGGTGG");
        checkPossibilities(calculator("WW", "T", ""), "GGTGG");
        checkPossibilities(calculator("WW", "", "G"), "TGGTG");
        checkPossibilities(calculator("WW", "T", "G"), "GGTG");
        
        checkPossibilities(calculator("YW", "", ""), "TATTGG,TACTGG");
        checkPossibilities(calculator("YW", "T", ""), "ATTGG,ACTGG");
        checkPossibilities(calculator("YW", "", "G"), "TATTG,TACTG");
        checkPossibilities(calculator("YW", "T", "G"), "ATTG,ACTG");
        checkPossibilities(calculator("YW", "TA", "G"), "TTG,CTG");
        checkPossibilities(calculator("YW", "TA", "GG"), "TT,CT");

        checkPossibilities(calculator("EF", "", ""), "GAATTT,GAATTC,GAGTTT,GAGTTC");
        checkPossibilities(calculator("EF", "G", ""), "AATTT,AATTC,AGTTT,AGTTC");
        checkPossibilities(calculator("EF", "GA", ""), "ATTT,ATTC,GTTT,GTTC");
        checkPossibilities(calculator("EF", "GA", "C"), "ATT,GTT");
    }

    @Test
    public void manyAAs()
    {
        checkPossibilities(calculator("KND", "", "T"), "AAAAATGA,AAAAACGA,AAGAATGA,AAGAACGA");
        checkPossibilities(calculator("KND", "A", "T"), "AAAATGA,AAAACGA,AGAATGA,AGAACGA");

        Set<String> expected = new HashSet<>();
        expected.add("AAAAATGAT");
        expected.add("AAAAACGAT");
        expected.add("AAGAATGAT");
        expected.add("AAGAACGAT");
        expected.add("AAAAATGAC");
        expected.add("AAAAACGAC");
        expected.add("AAGAATGAC");
        expected.add("AAGAACGAC");
        Assert.assertEquals(expected, calculator("KND", "", "").allPossibleBaseSequences());
    }

    @Test
    public void candidateAlternativeCodonsTest()
    {
        checkAltCodons(calculator("W", "T", ""), "GG");
        checkAltCodons(calculator("W", "T", "G"), "G");
        checkAltCodons(calculator("W", "TG", ""), "G");
        checkAltCodons(calculator("W", "", "G"), "TG");
        checkAltCodons(calculator("W", "G", ""), "");
        checkAltCodons(calculator("W", "", "T"), "");

        checkAltCodons(calculator("WW", "", ""), "TGG:TGG");
        checkAltCodons(calculator("WW", "T", "G"), "GG:TG");
        checkAltCodons(calculator("WW", "T", "GG"), "GG:T");
        checkAltCodons(calculator("WW", "TG", "GG"), "G:T");
        checkAltCodons(calculator("WW", "TG", "C"), "G: ");
        checkAltCodons(calculator("WW", "C", ""), " :TGG");
        checkAltCodons(calculator("WW", "C", "C"), " : ");

        checkAltCodons(calculator("AW", "", ""), "GCT,GCC,GCA,GCG:TGG");
        checkAltCodons(calculator("AW", "G", ""), "CT,CC,CA,CG:TGG");
        checkAltCodons(calculator("AW", "GC", ""), "T,C,A,G:TGG");
        checkAltCodons(calculator("AW", "GC", "G"), "T,C,A,G:TG");
        checkAltCodons(calculator("AW", "GC", "GG"), "T,C,A,G:T");

        checkAltCodons(calculator("WA", "", ""), "TGG:GCT,GCC,GCA,GCG");
        checkAltCodons(calculator("WA", "T", ""), "GG:GCT,GCC,GCA,GCG");
        checkAltCodons(calculator("WA", "TG", ""), "G:GCT,GCC,GCA,GCG");
        checkAltCodons(calculator("WA", "TG", "G"), "G:GC");
        checkAltCodons(calculator("WA", "TG", "CG"), "G:G");

        checkAltCodons(calculator("AA", "", ""), "GCT,GCC,GCA,GCG:GCT,GCC,GCA,GCG");
        checkAltCodons(calculator("AA", "G", ""), "CT,CC,CA,CG:GCT,GCC,GCA,GCG");
        checkAltCodons(calculator("AA", "G", "G"), "CT,CC,CA,CG:GC");

        checkAltCodons(calculator("KP", "", ""), "AAA,AAG:CCT,CCC,CCA,CCG");
        checkAltCodons(calculator("KP", "A", ""), "AA,AG:CCT,CCC,CCA,CCG");
        checkAltCodons(calculator("KP", "AA", ""), "A,G:CCT,CCC,CCA,CCG");
        checkAltCodons(calculator("KP", "AA", "A"), "A,G:CC");

        checkAltCodons(calculator("KPY", "", ""), "AAA,AAG:CCT,CCC,CCA,CCG:TAT,TAC");
        checkAltCodons(calculator("KPY", "A", ""), "AA,AG:CCT,CCC,CCA,CCG:TAT,TAC");
        checkAltCodons(calculator("KPY", "AA", ""), "A,G:CCT,CCC,CCA,CCG:TAT,TAC");
        checkAltCodons(calculator("KPY", "", "T"), "AAA,AAG:CCT,CCC,CCA,CCG:TA");
        checkAltCodons(calculator("KPY", "AA", "T"), "A,G:CCT,CCC,CCA,CCG:TA");
    }

    private NucleotidesCalculator calculator(String aminoAcids, String prefix, String suffix)
    {
        return new NucleotidesCalculator(aaSeq(aminoAcids), prefix, suffix);
    }

    private void checkPossibilities(NucleotidesCalculator calculator, String expectedSeparatedByCommas)
    {
        if(expectedSeparatedByCommas.trim().isBlank())
        {
            assertTrue(calculator.allPossibleBaseSequences().isEmpty());
            return;
        }
        Set<String> actual = calculator.allPossibleBaseSequences();
        Set<String> expected = css(expectedSeparatedByCommas);
        assertEquals(expected, actual);
    }

    private void checkAltCodons(NucleotidesCalculator di, String expectedAltCodonsSeparatedByColonsAndCommas)
    {
        List<Set<String>> actual = di.candidateAlternativeTruncatedCodons();
        String[] list = expectedAltCodonsSeparatedByColonsAndCommas.split(":");
        assertEquals(actual.size(), list.length);
        for(int i = 0; i < actual.size(); ++i) {
            if (list[i].isBlank())
            {
                assertEquals(0, actual.get(i).size());
            }
            else
            {
                Set<String> expected = Arrays.stream(list[i].trim().split(",")).collect(Collectors.toSet());
                assertEquals(expected, actual.get(i));
            }
        }
    }

}
