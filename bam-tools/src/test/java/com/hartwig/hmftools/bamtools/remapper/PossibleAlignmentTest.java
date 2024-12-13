package com.hartwig.hmftools.bamtools.remapper;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Test;

public class PossibleAlignmentTest extends RemapperTestBase
{
    @Test
    public void mateMustBeChr6()
    {
        BwaMemAlignment mate = bwa("16,8,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        AlternativeAlignment alternativeAlignment = new AlternativeAlignment("chr6", 31354375, Orientation.FORWARD, "151M", 50);
        Assert.assertThrows(IllegalArgumentException.class, () ->
                new PossibleAlignment(mate, alternativeAlignment));
    }

    @Test
    public void alternativeAlignmentMustBeChr6()
    {
        BwaMemAlignment mate = bwa("16,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        AlternativeAlignment alternativeAlignment = new AlternativeAlignment("chr8", 31354375, Orientation.FORWARD, "151M", 50);
        Assert.assertThrows(IllegalArgumentException.class, () ->
                new PossibleAlignment(mate, alternativeAlignment));
    }

    @Test
    public void orientationsMustBeCompatible()
    {
        BwaMemAlignment mate = bwa("16,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        AlternativeAlignment alternativeAlignment = new AlternativeAlignment("chr6", 31354375, Orientation.REVERSE, "151M", 50);
        Assert.assertThrows(IllegalArgumentException.class, () ->
                new PossibleAlignment(mate, alternativeAlignment));

        BwaMemAlignment mate2 = bwa("35,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        AlternativeAlignment alternativeAlignment2 = new AlternativeAlignment("chr6", 31354375, Orientation.FORWARD, "151M", 50);
        Assert.assertThrows(IllegalArgumentException.class, () ->
                new PossibleAlignment(mate2, alternativeAlignment2));
    }

    @Test
    public void distanceFromMate()
    {
        BwaMemAlignment mate = bwa("97,5,29942927,29943078,0,151,60,7,116,78,151M,26G10G16T1G39C18G13G21,null,5,29943327,551");
        // chr6,-29888542,15S136M,7
        AlternativeAlignment alternativeAlignment = new AlternativeAlignment("chr6", 29888542, Orientation.REVERSE, "15S136M", 7);
        PossibleAlignment possibleAlignment = new PossibleAlignment(mate, alternativeAlignment);
        assertEquals(29942927 - 29888542, possibleAlignment.distanceFromMate());

        BwaMemAlignment mate2 = bwa("145,5,29945301,29945452,0,151,60,2,141,39,151M,19C38A92,null,5,29945035,-417");
        // chr6,+29830355,151M,7
        AlternativeAlignment alternativeAlignment2 = new AlternativeAlignment("chr6", 29830355, Orientation.FORWARD, "151M", 7);
        PossibleAlignment possibleAlignment2 = new PossibleAlignment(mate2, alternativeAlignment2);
        assertEquals(29945301 - 29830355, possibleAlignment2.distanceFromMate());
    }

    @Test
    public void compareTest()
    {
        BwaMemAlignment mate = bwa("145,5,29945301,29945452,0,151,60,2,141,39,151M,19C38A92,null,5,29945035,-417");
        // chr6,+29830355,151M,7;chr6,+29890241,151M,7;chr6,+30009122,151M,7;
        AlternativeAlignment aa1 = new AlternativeAlignment("chr6", 29830355, Orientation.FORWARD, "151M", 7);
        PossibleAlignment pa1 = new PossibleAlignment(mate, aa1);
        PossibleAlignment pa11 = new PossibleAlignment(mate, aa1);
        AlternativeAlignment aa2 = new AlternativeAlignment("chr6", 29890241, Orientation.FORWARD, "151M", 7);
        PossibleAlignment pa2 = new PossibleAlignment(mate, aa2);
        AlternativeAlignment aa3 = new AlternativeAlignment("chr6", 30009122, Orientation.FORWARD, "151M", 7);
        PossibleAlignment pa3 = new PossibleAlignment(mate, aa3);
        assertEquals(0, pa11.compareTo(pa1));
        assertTrue(pa1.compareTo(pa2) > 0);
        assertTrue(pa2.compareTo(pa1) < 0);
        assertTrue(pa3.compareTo(pa2) > 0);
    }
}
