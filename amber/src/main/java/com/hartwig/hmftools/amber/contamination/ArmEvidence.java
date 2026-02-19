package com.hartwig.hmftools.amber.contamination;

import java.util.Objects;

import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.jetbrains.annotations.NotNull;

class ArmEvidence implements Comparable<ArmEvidence>
{
    private final ChrArm chrArm;
    private int TotalPoints = 0;
    private int EvidencePoints = 0;

    ArmEvidence(final ChrArm chrArm)
    {
        this.chrArm = chrArm;
    }

    public void register(boolean isEvidence)
    {
        TotalPoints++;
        if(isEvidence)
        {
            EvidencePoints++;
        }
    }

    public int totalPoints()
    {
        return TotalPoints;
    }

    public int evidencePoints()
    {
        return EvidencePoints;
    }

    @Override
    public int compareTo(@NotNull final ArmEvidence o)
    {
        int result = Double.compare(ratio(), o.ratio());
        if(result == 0)
        {
            result = chrArm.compareTo(o.chrArm);
        }
        return result;
    }

    public double ratio()
    {
        if(TotalPoints == 0)
        {
            return Double.NaN;
        }
        return (double) EvidencePoints / (double) TotalPoints;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ArmEvidence that = (ArmEvidence) o;
        return Objects.equals(chrArm, that.chrArm);
    }

    @Override
    public int hashCode()
    {
        return Objects.hashCode(chrArm);
    }
}
