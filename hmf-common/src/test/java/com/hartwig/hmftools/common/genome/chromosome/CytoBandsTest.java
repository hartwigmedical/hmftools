package com.hartwig.hmftools.common.genome.chromosome;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.junit.Test;

public class CytoBandsTest
{
    @Test
    public void testCytoBands37()
    {
        CytoBands cytoBands = new CytoBands(RefGenomeVersion.V37);
        assertEquals("p36.33", cytoBands.getCytoBandName("chr1", 1000));
    }

    @Test
    public void testCytoBands38()
    {
        CytoBands cytoBands = new CytoBands(RefGenomeVersion.V38);
        assertEquals("p36.33", cytoBands.getCytoBandName("chr1", 1000));
    }
}
