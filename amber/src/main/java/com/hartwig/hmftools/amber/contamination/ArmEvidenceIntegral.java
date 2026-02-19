package com.hartwig.hmftools.amber.contamination;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.Set;

class ArmEvidenceIntegral
{
    private final List<ArmEvidence> mArmEvidence;

    public ArmEvidenceIntegral(Set<ArmEvidence> armEvidence)
    {
        mArmEvidence = armEvidence.stream().sorted().collect(toList());
    }

    public int totalHits()
    {
        return mArmEvidence.stream().mapToInt(ArmEvidence::evidencePoints).sum();
    }

    public int totalPoints()
    {
        return mArmEvidence.stream().map(ArmEvidence::totalPoints).reduce(0, Integer::sum);
    }

    public double value()
    {
        double heightSoFar = 0.0;
        double area = 0.0;
        for(ArmEvidence armEvidence : mArmEvidence)
        {
            double rectangleSize = armEvidence.totalPoints() * heightSoFar;
            area += rectangleSize;
            double triangleSize = armEvidence.evidencePoints() * armEvidence.totalPoints() / 2.0;
            area += triangleSize;
            heightSoFar += armEvidence.evidencePoints();
        }
        return area;
    }
}
