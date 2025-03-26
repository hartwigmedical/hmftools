package com.hartwig.hmftools.pavereverse;

import org.junit.Test;

public final class DnaVariantsTest extends ReversePaveTestBase
{
    @Test
    public void baseIsUpstreamOfFirstUtrExonReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant("TERT", "ENST00000310581", "c.-146C>T");
        check(bsc, "G", "A", "chr5", 1295135);
    }

    @Test
    public void baseIsUpstreamOfFirstUtrExon()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-80G>A");
        check(bsc, "G", "A", "chr7", 143381345); // Sanity check, this is actually the first 5' UTR exonic base.

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-81C>A");
        check(bsc, "C", "A", "chr7", 143381344);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-85G>A");
        check(bsc, "G", "A", "chr7", 143381340);
    }

    @Test
    public void calculateVariantsForSuppliedTranscriptId()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.5C>A");
        check(bsc, "C", "A", "chr7", 143381572 + 5 - 1);

        // Now with a different transcript.
        bsc = reversePave.calculateDnaVariant(zyx, "ENST00000354434", "c.2C>A");
        check(bsc, "C", "A", "chr7", 143381575 + 2 - 1);
    }

    @Test
    public void mnvSubstitution()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "210_212delCTTinsGAAC");
        check(bsc, "CTT", "GAAC", "chr7", 143382249);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "210delCinsGAAC");
        check(bsc, "C", "GAAC", "chr7", 143382249);
    }

    @Test
    public void mnvSubstitutionReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.135_138delGGAGinsTTTT");
        check(bsc, "CTCC", "AAAA", "chr7", 140_924_566);
    }

    @Test
    public void complexMNV()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "211_218delinsGAAC");
        check(bsc, "TTTCCCCT", "GAAC", "chr7", 143382250);
    }

    @Test
    public void complexMnvReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.135_138delinsTTTT");
        check(bsc, "CTCC", "AAAA", "chr7", 140_924_566);
    }

    @Test
    public void substitutionAtStartOfExon()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.209A>T");
        check(bsc, "A", "T", "chr7", 143382248);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.409C>G");
        check(bsc, "C", "G", "chr7", 143382593);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.1615G>T");
        check(bsc, "G", "T", "chr7", 143390578);
    }

    @Test
    public void substitutionAtStartOfExonReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.1A>C");
        check(bsc, "T", "G", "chr7", 140_924_703);

        bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.10C>T");
        check(bsc, "G", "A", "chr7", 140_924_694);
    }

    @Test
    public void substitutionAtEndOfExon()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.208G>T");
        check(bsc, "G", "T", "chr7", 143381779);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.408G>T");
        check(bsc, "G", "T", "chr7", 143382447);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.1719A>G");
        check(bsc, "A", "G", "chr7", 143390682);
    }

    @Test
    public void substitutionAtEndOfExonReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.138G>T");
        check(bsc, "C", "A", "chr7", 140_924_566);

        bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.608G>C");
        check(bsc, "C", "G", "chr7", 140_808_892);
    }

    @Test
    public void substitutionUpstreamOfStart()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-1C>A");
        check(bsc, "C", "A", "chr7", 143381572 - 1);
    }

    @Test
    public void substitutionUpstreamOfStartReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.-10C>T");
        check(bsc, "G", "A", "chr7", 140_924_713);
    }

    @Test
    public void substitutionDownstreamOfStop()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.*1G>A");
        check(bsc, "G", "A", "chr7", 143390682 + 1);
    }

    @Test
    public void substitutionAfterExonEndTest()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.208+1G>T");
        check(bsc, "G", "T", "chr7", 143381780);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.408+4C>A");
        check(bsc, "C", "A", "chr7", 143382451);
    }

    @Test
    public void substitutionBeforeExonStartTest()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.209-1G>A");
        check(bsc, "G", "A", "chr7", 143382247);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.409-3T>A");
        check(bsc, "T", "A", "chr7", 143382590);
    }

    @Test
    public void substitutionIn5PrimeUtrIntron()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-15-1G>A");
        check(bsc, "G", "A", "chr7", 143381556);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-16+1G>A");
        check(bsc, "G", "A", "chr7", 143381410);
    }

    @Test
    public void substitutionIn3PrimeUtrIntron()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(tatdn2, tatdn2Canonical, "c.*38+1G>A");
        check(bsc, "G", "A", "chr3", 10_279_064);

        bsc = reversePave.calculateDnaVariant(tatdn2, tatdn2Canonical, "c.*39-3C>T");
        check(bsc, "C", "T", "chr3", 10_279_218);
    }

    @Test
    public void deletionOfSingleBase()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(tatdn2, tatdn2Canonical, "c.977delA");
        check(bsc, "CA", "C", "chr3", 10_270_158);
    }

    @Test
    public void deletionOfRange()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(tatdn2, tatdn2Canonical, "c.974_977delAGCA");
        check(bsc, "GAGCA", "G", "chr3", 10_270_155);
    }

    @Test
    public void deletionOfRangeReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "602_606delTAGGA");
        check(bsc, "CATCCT", "C", "chr7", 140_808_893);
    }

    @Test
    public void duplication()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-15-2dupA");
        check(bsc, "C", "CA", "chr7", 143381554);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-15-5_-15-2dup");
        check(bsc, "C", "CCCCA", "chr7", 143381551);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-15-5_-14dup");
        check(bsc, "C", "CCCCAGCA", "chr7", 143381551);
    }

    @Test
    public void duplicationReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "606dupT");
        check(bsc, "C", "CA", "chr7", 140_808_893);

        bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "608+1dupG");
        check(bsc, "A", "AC", "chr7", 140_808_890);

        bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "608+3_608+4dupAT");
        check(bsc, "C", "CAT", "chr7", 140_808_887);
    }

    @Test
    public void insertion()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-15-3_-15-2insGGATCC");
        check(bsc, "C", "CGGATCC", "chr7", 143_381_554);
    }

    @Test
    public void insertionReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.*2_*3insGGATTC");
        check(bsc, "G", "GGAATCC", "chr7", 140_734_594);
    }

    @Test
    public void deletionsAreLeftAligned()
    {
        //               143381572
        //               |
        //               1     7
        // AGCCCGGCCCGGCCATGGCGGCCCCCCG...

        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.7del");
        check(bsc, "CG", "C", "chr7", 143381576);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.13del");
        check(bsc, "GC", "G", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.12_13del");
        check(bsc, "GCC", "G", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.11_13del");
        check(bsc, "GCCC", "G", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.10_13del");
        check(bsc, "GCCCC", "G", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.9_13del");
        check(bsc, "GCCCCC", "G", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-8_-4del");
        check(bsc, "AGCCCG", "A", "chr7", 143381558);
    }

    @Test
    public void deletionsOnReverseStrandAreAlreadyLeftAligned()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "592delT");
        check(bsc, "TA", "T", "chr7", 140_808_907);
    }

    @Test
    public void insertionsAreLeftAligned()
    {
        // 143381572
        // |
        // 1     7
        // ATGGCGGCCCCCCG...
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.7_8insGGG");
        check(bsc, "C", "CGGG", "chr7", 143381576);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.13_14insC");
        check(bsc, "G", "GC", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.13_14insCCC");
        check(bsc, "G", "GCCC", "chr7", 143381578);
    }

    @Test
    public void insertionsOnReverseStrandAreAlreadyLeftAligned()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "600_601insTTT");
        check(bsc, "G", "GAAA", "chr7", 140_808_899);
    }

    @Test
    public void duplicationsAreLeftAligned()
    {
        //               143381572
        //               |
        //               1     7
        // AGCCCGGCCCGGCCATGGCGGCCCCCCG...

        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.7dup");
        check(bsc, "C", "CG", "chr7", 143381576);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.12_13dup");
        check(bsc, "G", "GCC", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.11_13dup");
        check(bsc, "G", "GCCC", "chr7", 143381578);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-8_-4dup");
        check(bsc, "A", "AGCCCG", "chr7", 143381558);
    }

    @Test
    public void duplicationsOnReverseStrandAreAlreadyLeftAligned()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "585_586dupT");
        check(bsc, "G", "GCA", "chr7", 140_808_913);
    }
}