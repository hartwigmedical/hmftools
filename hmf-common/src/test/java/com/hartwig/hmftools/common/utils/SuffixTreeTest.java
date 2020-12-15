package com.hartwig.hmftools.common.utils;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

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
        final String read = "AATGTAGGTGCTGCTGTGAAGGGATTTAGCAGATATAATTAAGGGTCTCAATTAGTTGACTTTATGCTGCGTTTATCCTGCTTGGACTGTCCTAATCAGGTGAGCCCTTGAAAGGACTGGGTTCTTCATGAGCATAGAGACTTACAGTGTG";
        SuffixTree tree = new SuffixTree(read);
        for (int i = 0; i < read.length() - 10; i++) {
            assertTrue(tree.contains(read.substring(i, i + 10)));
        }
    }

}
