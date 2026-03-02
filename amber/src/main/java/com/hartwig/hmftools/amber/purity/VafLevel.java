package com.hartwig.hmftools.amber.purity;

import java.util.Set;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class VafLevel
{
    private final double Level;

    public VafLevel(double Level)
    {
        this.Level = Level;
    }

    public boolean hasSufficientDepthForEventDetection(PositionEvidence evidence)
    {
        return false;
    }

    public void test(PositionEvidence evidence)
    {

    }

    public double perArmConsistencyFactor(final ChrArmLocator chrArmLocator)
    {
        return 0;
    }

    public double perMutationTypeConsistencyFactor()
    {
        return 0;
    }

    public Set<PositionEvidence> homozygousEvidencePoints()
    {
        return null;
    }

    public Set<PositionEvidence> heterozygousEvidencePoints()
    {
        return null;
    }

    public int numberOfCapturedEvidencePoints()
    {
        return 0;
    }

    public double vaf()
    {
        return 0;
    }

    public double homozygousProportion()
    {
        return 0;
    }
}
