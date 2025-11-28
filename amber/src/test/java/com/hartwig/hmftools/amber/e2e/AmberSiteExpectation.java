package com.hartwig.hmftools.amber.e2e;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.amber.AmberBAF;

class AmberSiteExpectation
{
    private final int Position;
    private final int RefReadsCount;
    private final int AltReadsCount;
    private final int Depth;

    AmberSiteExpectation(int position, int refCount, int altCount, int depth)
    {
        Position = position;
        RefReadsCount = refCount;
        AltReadsCount = altCount;
        Depth = depth;
    }

    boolean aligns(AmberBAF baf)
    {
        return baf.Position == Position;
    }

    void checkTumorBaf(AmberBAF baf)
    {
        assertEquals((float) AltReadsCount / (float) (AltReadsCount + RefReadsCount), baf.tumorBAF(), 0.001);
        assertEquals(Depth, baf.TumorDepth);
        assertEquals(Position, baf.Position);
    }
}
