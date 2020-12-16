package com.hartwig.hmftools.common.utils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.codon.Codons;

import org.junit.Test;

public class SuffixTreeTest {

    @Test
    public void testHavanaBanana() {
        SuffixTree tree = new SuffixTree("havanabanana");
        assertTrue(tree.contains("banana"));
        assertFalse(tree.contains("bananad"));
    }

    @Test
    public void testRealExample() {

        final String dna = "AATGTAGGTGCTGCTGTGAAGGGATTTAGCAGATATAATTAAGGGTCTCAATTAGTTGACTTTATGCTGCGTTTATCCTGCTTGGACTGTCCTAATCAGGTGAGCCCTTGAAAGGACTGGGTTCTTCATGAGCATAGAGACTTACAGTGTG";
        final String aminoAcids = Codons.asCodonString(dna);

        SuffixTree tree = new SuffixTree(aminoAcids);
        for (int i = 0; i < aminoAcids.length() - 10; i++) {
            assertTrue(tree.contains(aminoAcids.substring(i, i + 10)));
        }
    }


    @Test
    public void testEndsWith() {
        final String read = "NVTENCGAAXAVKGFSRYNXTENGAA";
        SuffixTree tree = new SuffixTree(read);
        assertEquals(0, tree.endsWith("GA"));
        assertEquals(0, tree.endsWith("GA"));

        assertEquals(1, tree.endsWith("A"));
        assertEquals(2, tree.endsWith("AA"));
        assertEquals(3, tree.endsWith("GAA"));

        assertEquals(1, tree.endsWith("AT"));
        assertEquals(2, tree.endsWith("AAG"));
        assertEquals(3, tree.endsWith("GAAH"));

        assertEquals(4, tree.endsWith("NGAATAF"));
        assertEquals(0, tree.endsWith("TEN"));
        assertEquals(0, tree.endsWith("AXA"));
    }

}
