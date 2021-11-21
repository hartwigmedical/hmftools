package com.hartwig.hmftools.teal;

import static com.hartwig.hmftools.teal.TeloConstants.CANONICAL_TELOMERE_SEQ;
import static com.hartwig.hmftools.teal.TeloUtils.hasTelomericContent;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class TestTeloUtils
{
    @Test
    public void testTelomericContent()
    {
        String readBases = "AGCT" + CANONICAL_TELOMERE_SEQ + "AGCT" + CANONICAL_TELOMERE_SEQ + CANONICAL_TELOMERE_SEQ + "GG";
        assertTrue(hasTelomericContent(readBases));

        readBases = "AGCT" + CANONICAL_TELOMERE_SEQ + "AGCT" + CANONICAL_TELOMERE_SEQ + "GG";
        assertFalse(hasTelomericContent(readBases));
    }

    @Test
    public void testLikelyTelomeric()
    {
        String readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG";
        String readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC";

        assertTrue(TeloUtils.isLikelyGTelomeric(readBasesG));
        assertFalse(TeloUtils.isLikelyCTelomeric(readBasesG));

        assertFalse(TeloUtils.isLikelyGTelomeric(readBasesC));
        assertTrue(TeloUtils.isLikelyCTelomeric(readBasesC));

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";

    }

    @Test
    public void testGTelomHexamer()
    {
        // canonical
        assertTrue(TeloUtils.isGTeloHexamer("TTAGGG", 0));
        assertTrue(TeloUtils.isGTeloHexamer("GGGGTTAGGGTT", 4));
        assertFalse(TeloUtils.isGTeloHexamer("GGGGTTAGGGTT", 3));

        // non canonical
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518223/
        // TCAGGG is most important type for both ALT maintained telomere
        assertTrue(TeloUtils.isGTeloHexamer("TCAGGG", 0));
        assertTrue(TeloUtils.isGTeloHexamer("GTAGGG", 0));
        assertTrue(TeloUtils.isGTeloHexamer("GGGTCAGGG", 3));
        assertTrue(TeloUtils.isGTeloHexamer("GGGGTAGGG", 3));

        // but to protect again poly G tail we want to check for number of G is not too many
        assertFalse(TeloUtils.isGTeloHexamer("GGAGGG", 0));
        assertFalse(TeloUtils.isGTeloHexamer("GTGGGG", 0));
    }

    @Test
    public void testCTelomHexamer()
    {
        // canonical
        assertTrue(TeloUtils.isCTeloHexamer("CCCTAA", 0));
        assertTrue(TeloUtils.isCTeloHexamer("GGGGCCCTAAG", 4));
        assertFalse(TeloUtils.isCTeloHexamer("GGGGCCCTAAG", 3));

        // non canonical
        // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3518223/
        // TCAGGG is most important type for both ALT maintained telomere
        assertTrue(TeloUtils.isCTeloHexamer("CCCTGA", 0));
        assertTrue(TeloUtils.isCTeloHexamer("CCCTAC", 0));
        assertTrue(TeloUtils.isCTeloHexamer("CCCCCCTGA", 3));
        assertTrue(TeloUtils.isCTeloHexamer("CCCCCCTACC", 3));

        // since for the sense strand we protect against poly G we also try to be consistent
        // and not let too many C on the other strand
        assertFalse(TeloUtils.isCTeloHexamer("CCCTCC", 0));
        assertFalse(TeloUtils.isCTeloHexamer("CCCCAC", 0));
    }

    @Test
    public void testStrictlyGTeloermic()
    {
        // canonical
        assertTrue(TeloUtils.isStrictlyGTelomeric("TTAGGGTTAGGGTTAGGG"));

        // non canonical
        assertTrue(TeloUtils.isStrictlyGTelomeric("TCAGGGTGAGGGTTAGGG"));

        // add 1 incomplete hexamer on both ends
        assertTrue(TeloUtils.isStrictlyGTelomeric("GTTAGGGTTAGGGTTAGGGT"));

        // 5 incomplete on both ends
        assertTrue(TeloUtils.isStrictlyGTelomeric("TAGGGTTAGGGTTAGGGTTAGGGTTAGG"));

        // test one hexamer is not telomeric
        assertFalse(TeloUtils.isStrictlyGTelomeric("TTAGGGTTAGGGTGGGGG"));

        // test just some random sequence
        assertFalse(TeloUtils.isStrictlyGTelomeric("TAGTGGATTTGGAGGTTATG"));
    }

    @Test
    public void testStrictlyCTeloermic()
    {
        // canonical
        assertTrue(TeloUtils.isStrictlyCTelomeric("CCCTAACCCTAACCCTAA"));

        // add 1 incomplete hexamer on both ends
        assertTrue(TeloUtils.isStrictlyCTelomeric("ACCCTAACCCTAACCCTAAC"));

        // 5 incomplete on both ends
        assertTrue(TeloUtils.isStrictlyCTelomeric("CCTAACCCTAACCCTAACCCTAACCCTA"));

        // test one hexamer is not telomeric
        assertFalse(TeloUtils.isStrictlyCTelomeric("CCCTAACCCTAACCCTCC"));

        // test just some random sequence
        assertFalse(TeloUtils.isStrictlyCTelomeric("TAGTGGATTTGGAGGTTATG"));
    }
}
