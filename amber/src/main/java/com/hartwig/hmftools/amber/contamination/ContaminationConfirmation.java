package com.hartwig.hmftools.amber.contamination;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class ContaminationConfirmation
{
    private final ContaminationPredicate Classifier;
    private final ChrArmLocator mChrArmLocator;
    private final Map<ChrArm, Integer> ArmToEvidenceCount = new HashMap<>();
    private final Map<ChrArm, Integer> ArmToNumberOfPoints = new HashMap<>();

    public ContaminationConfirmation(final ContaminationPredicate classifier, final ChrArmLocator chrArmLocator)
    {
        Classifier = classifier;
        mChrArmLocator = chrArmLocator;
    }

    public void offer(TumorContamination contamination)
    {
        ChrArm chrArm = mChrArmLocator.map(contamination.Chromosome, contamination.Position);
        ArmToNumberOfPoints.put(chrArm, ArmToNumberOfPoints.getOrDefault(chrArm, 0) + 1);
        if(Classifier.test(contamination))
        {
            ArmToEvidenceCount.put(chrArm, ArmToEvidenceCount.getOrDefault(chrArm, 0) + 1);
        }
    }

    public double gini()
    {
        return 0.0; // todo
    }
}

class ArmEvidence
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
