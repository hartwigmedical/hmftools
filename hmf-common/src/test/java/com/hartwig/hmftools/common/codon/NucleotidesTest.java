package com.hartwig.hmftools.common.codon;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBasesInPlace;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class NucleotidesTest
{
    @Test
    public void testReverseComplementBasesString()
    {
        assertEquals("", reverseComplementBases(""));
        assertEquals("T", reverseComplementBases("A"));
        assertEquals("N", reverseComplementBases("N"));
        assertEquals("NNAAGGCCTT", reverseComplementBases("AAGGCCTTNN"));
        assertEquals("NAAGGCCTT", reverseComplementBases("AAGGCCTTN"));
    }

    @Test
    public void testReverseComplementBasesBytes()
    {
        assertEquals("", new String(reverseComplementBases("".getBytes())));
        assertEquals("T", new String(reverseComplementBases("A".getBytes())));
        assertEquals("N", new String(reverseComplementBases("N".getBytes())));
        assertEquals("NNAAGGCCTT", new String(reverseComplementBases("AAGGCCTTNN".getBytes())));
        assertEquals("NAAGGCCTT", new String(reverseComplementBases("AAGGCCTTN".getBytes())));
    }

    @Test
    public void testReverseComplementBasesInPlace()
    {
        {
            byte[] sequence = "".getBytes();
            reverseComplementBasesInPlace(sequence, 0, 0);
            assertEquals("", new String(sequence));
        }

        {
            byte[] sequence = "A".getBytes();
            reverseComplementBasesInPlace(sequence, 0, 1);
            assertEquals("T", new String(sequence));
        }

        {
            byte[] sequence = "N".getBytes();
            reverseComplementBasesInPlace(sequence, 0, 1);
            assertEquals("N", new String(sequence));
        }

        {
            byte[] sequence = "AAGGCCTTNN".getBytes();
            reverseComplementBasesInPlace(sequence, 0, 10);
            assertEquals("NNAAGGCCTT", new String(sequence));
        }

        {
            byte[] sequence = "AAGGCCTTN".getBytes();
            reverseComplementBasesInPlace(sequence, 0, 9);
            assertEquals("NAAGGCCTT", new String(sequence));
        }

        {
            byte[] sequence = "AAGGCCTTNN".getBytes();
            reverseComplementBasesInPlace(sequence, 1, 4);
            assertEquals("AGCCTCTTNN", new String(sequence));
        }

        {
            byte[] sequence = "AAGGCCTTNN".getBytes();
            reverseComplementBasesInPlace(sequence, 4, 4);
            assertEquals("AAGGAAGGNN", new String(sequence));
        }

        {
            byte[] sequence = "AAGGCCTTNN".getBytes();
            reverseComplementBasesInPlace(sequence, 1, 5);
            assertEquals("AGGCCTTTNN", new String(sequence));
        }

        {
            byte[] sequence = "AAGGCCTTNN".getBytes();
            reverseComplementBasesInPlace(sequence, 2, 5);
            assertEquals("AAAGGCCTNN", new String(sequence));
        }
    }
}
