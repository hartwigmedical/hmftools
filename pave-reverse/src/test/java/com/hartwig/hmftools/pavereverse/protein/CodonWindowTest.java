package com.hartwig.hmftools.pavereverse.protein;

import java.util.List;

import com.hartwig.hmftools.pavereverse.base.ChangeContextBuilder;

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
        Assert.assertEquals(0, builder.Data.ExonIndex);
        Assert.assertEquals(9, builder.Data.ChangeStart);
        Assert.assertEquals(20, builder.Data.ChangeEnd);
        Assert.assertEquals(0, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(2, builder.Data.PaddingInNextExon);
        Assert.assertEquals(1, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        //-|--...--|-
        bw = new CodonWindow(51,4);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.Data.ExonIndex);
        Assert.assertEquals(50, builder.Data.ChangeStart);
        Assert.assertEquals(61, builder.Data.ChangeEnd);
        Assert.assertEquals(1, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.Data.PaddingInNextExon);
        Assert.assertEquals(35, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        //--|--...---|
        bw = new CodonWindow(95,5);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(2, builder.Data.ExonIndex);
        Assert.assertEquals(82, builder.Data.ChangeStart);
        Assert.assertEquals(96, builder.Data.ChangeEnd);
        Assert.assertEquals(2, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.Data.PaddingInNextExon);
        Assert.assertEquals(68, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(175,2);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(5, builder.Data.ExonIndex);
        Assert.assertEquals(22, builder.Data.ChangeStart);
        Assert.assertEquals(27, builder.Data.ChangeEnd);
        Assert.assertEquals(2, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.Data.PaddingInNextExon);
        Assert.assertEquals(168, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    @Test
    public void windowStartsNearStartOfExon()
    {
        CodonWindow bw = new CodonWindow(35,10);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.Data.ExonIndex);
        Assert.assertEquals(2, builder.Data.ChangeStart);
        Assert.assertEquals(31, builder.Data.ChangeEnd);
        Assert.assertEquals(1, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.Data.PaddingInNextExon);
        Assert.assertEquals(35, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(101,8);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(3, builder.Data.ExonIndex);
        Assert.assertEquals(0, builder.Data.ChangeStart);
        Assert.assertEquals(23, builder.Data.ChangeEnd);
        Assert.assertEquals(0, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.Data.PaddingInNextExon);
        Assert.assertEquals(101, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    @Test
    public void windowStartsAtEndOfExon()
    {
        // |--- --- ... --- +|++ --- ... --- --|- --- ...   The window in this example is the + chars
        CodonWindow bw = new CodonWindow(34,1);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.Data.ExonIndex);
        Assert.assertEquals(0, builder.Data.ChangeStart);
        Assert.assertEquals(1, builder.Data.ChangeEnd);
        Assert.assertEquals(1, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.Data.PaddingInNextExon);
        Assert.assertEquals(35, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);
        Assert.assertNotNull(builder.CompanionData);
        Assert.assertEquals(0, builder.CompanionData.ExonIndex);
        Assert.assertEquals(99, builder.CompanionData.ChangeStart);
        Assert.assertEquals(99, builder.CompanionData.ChangeEnd);
        Assert.assertEquals(0, builder.CompanionData.PaddingInPreviousExon);
        Assert.assertEquals(2, builder.CompanionData.PaddingInNextExon);
        Assert.assertEquals(1, builder.CompanionData.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(67,2);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(2, builder.Data.ExonIndex);
        Assert.assertEquals(0, builder.Data.ChangeStart);
        Assert.assertEquals(3, builder.Data.ChangeEnd);
        Assert.assertEquals(2, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(0, builder.Data.PaddingInNextExon);
        Assert.assertEquals(68, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    @Test
    public void windowEndsAtStartOfExon()
    {
        CodonWindow bw = new CodonWindow(33,2);
        ChangeContextBuilder builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(0, builder.Data.ExonIndex);
        Assert.assertEquals(96, builder.Data.ChangeStart);
        Assert.assertEquals(99, builder.Data.ChangeEnd);
        Assert.assertEquals(0, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(2, builder.Data.PaddingInNextExon);
        Assert.assertEquals(1, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);

        bw = new CodonWindow(61,7);
        builder = bw.seekExonLocation(exonLengths);
        Assert.assertEquals(1, builder.Data.ExonIndex);
        Assert.assertEquals(80, builder.Data.ChangeStart);
        Assert.assertEquals(99, builder.Data.ChangeEnd);
        Assert.assertEquals(1, builder.Data.PaddingInPreviousExon);
        Assert.assertEquals(1, builder.Data.PaddingInNextExon);
        Assert.assertEquals(35, builder.Data.AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }
}
