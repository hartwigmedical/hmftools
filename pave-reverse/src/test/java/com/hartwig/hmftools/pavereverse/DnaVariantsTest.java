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
    public void substitutionUpstreamOfStart()
    {
        BaseSequenceChange bsc = reversePave.calculateDnaVariant(zyx, zyxCanonical, "c.-1C>A");
        check(bsc, "C", "A", "chr7", 143381572 - 1);
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

    }
}