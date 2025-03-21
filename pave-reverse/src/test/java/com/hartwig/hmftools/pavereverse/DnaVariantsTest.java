package com.hartwig.hmftools.pavereverse;

import org.junit.Test;

public final class DnaVariantsTest extends ReversePaveTestBase
{

    @Test
    public void calculateVariantsForSuppliedTranscriptId()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.5C>A");
        check(bsc, "C", "A", "chr7", 143381572 + 5 - 1);
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

        bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.10G>C");
        check(bsc, "C", "G", "chr7", 140_924_694);
    }

    @Test
    public void substitutionAtEndOfExon()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.208A>T");
        check(bsc, "A", "T", "chr7", 143381779);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.408G>T");
        check(bsc, "G", "T", "chr7", 143382447);

        bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.1719A>G");
        check(bsc, "A", "G", "chr7", 143390682);
    }

    @Test
    public void substitutionAtEndOfExonReverseStrand()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "c.138C>T");
        check(bsc, "G", "A", "chr7", 140_924_566);

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

        bsc = reversePave.calculateDnaVariant(braf, brafCanonical, "608dupG");
        check(bsc, "C", "CC", "chr7", 140_808_891);

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
}