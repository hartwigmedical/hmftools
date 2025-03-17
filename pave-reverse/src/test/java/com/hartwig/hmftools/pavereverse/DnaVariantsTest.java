package com.hartwig.hmftools.pavereverse;

import org.junit.Test;

public final class DnaVariantsTest extends ReversePaveTestBase
{
    private final String zyx = "ZYX";
    private final String zyxCanonical = "ENST00000322764";

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
}