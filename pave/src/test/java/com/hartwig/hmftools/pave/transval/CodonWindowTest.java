package com.hartwig.hmftools.pave.transval;

import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class CodonWindowTest
{
    List<Integer> exonLengths = List.of(100, 100, 100, 120, 80, 31, 78);

    @Test
    public void windowInExon()
    {
        // |--- ....--|-
        CodonWindow bw = new CodonWindow(4,4);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(0, builder.ExonIndex);
        Assert.assertEquals(9, builder.ChangeStart);
        Assert.assertEquals(20, builder.ChangeEnd);
        Assert.assertEquals(0, builder.PaddingInPreviousExon);
        Assert.assertEquals(2, builder.PaddingInNextExon);
        Assert.assertEquals(1, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        //-|--...--|-
        bw = new CodonWindow(51,4);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.ExonIndex);
        Assert.assertEquals(50, builder.ChangeStart);
        Assert.assertEquals(61, builder.ChangeEnd);
        Assert.assertEquals(1, builder.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.PaddingInNextExon);
        Assert.assertEquals(35, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        //--|--...---|
        bw = new CodonWindow(95,5);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(2, builder.ExonIndex);
        Assert.assertEquals(82, builder.ChangeStart);
        Assert.assertEquals(96, builder.ChangeEnd);
        Assert.assertEquals(2, builder.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.PaddingInNextExon);
        Assert.assertEquals(68, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(175,2);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(5, builder.ExonIndex);
        Assert.assertEquals(22, builder.ChangeStart);
        Assert.assertEquals(27, builder.ChangeEnd);
        Assert.assertEquals(2, builder.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.PaddingInNextExon);
        Assert.assertEquals(168, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    @Test
    public void windowStartsNearStartOfExon()
    {
        CodonWindow bw = new CodonWindow(35,10);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.ExonIndex);
        Assert.assertEquals(2, builder.ChangeStart);
        Assert.assertEquals(31, builder.ChangeEnd);
        Assert.assertEquals(1, builder.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.PaddingInNextExon);
        Assert.assertEquals(35, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(101,8);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(3, builder.ExonIndex);
        Assert.assertEquals(0, builder.ChangeStart);
        Assert.assertEquals(23, builder.ChangeEnd);
        Assert.assertEquals(0, builder.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.PaddingInNextExon);
        Assert.assertEquals(101, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    @Test
    public void windowStartsAtEndOfExon()
    {
        CodonWindow bw = new CodonWindow(34,1);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.ExonIndex);
        Assert.assertEquals(0, builder.ChangeStart);
        Assert.assertEquals(1, builder.ChangeEnd);
        Assert.assertEquals(1, builder.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.PaddingInNextExon);
        Assert.assertEquals(35, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(67,2);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(2, builder.ExonIndex);
        Assert.assertEquals(0, builder.ChangeStart);
        Assert.assertEquals(3, builder.ChangeEnd);
        Assert.assertEquals(2, builder.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.PaddingInNextExon);
        Assert.assertEquals(68, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    @Test
    public void windowEndsAtStartOfExon()
    {
        CodonWindow bw = new CodonWindow(33,2);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(0, builder.ExonIndex);
        Assert.assertEquals(96, builder.ChangeStart);
        Assert.assertEquals(99, builder.ChangeEnd);
        Assert.assertEquals(0, builder.PaddingInPreviousExon);
        Assert.assertEquals(2, builder.PaddingInNextExon);
        Assert.assertEquals(1, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(61,7);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.ExonIndex);
        Assert.assertEquals(80, builder.ChangeStart);
        Assert.assertEquals(99, builder.ChangeEnd);
        Assert.assertEquals(1, builder.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.PaddingInNextExon);
        Assert.assertEquals(35, builder.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }
}
