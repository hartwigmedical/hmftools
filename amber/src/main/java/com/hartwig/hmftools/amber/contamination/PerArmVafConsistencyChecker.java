package com.hartwig.hmftools.amber.contamination;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class PerArmVafConsistencyChecker
{
    public static VafConsistencyCheckResult calculateConfirmationFactor(final ChrArmLocator chrArmLocator,
            double vaf,
            List<TumorContamination> data)
    {
        VafPredicate classifier = new BinomialVafPredicate(vaf);
        PerArmVafConsistencyChecker confirmation = new PerArmVafConsistencyChecker(classifier, chrArmLocator);
        data.forEach(confirmation::offer);
        return confirmation.gini();
    }

    private final VafPredicate Classifier;
    private final ChrArmLocator mChrArmLocator;
    private final Map<ChrArm, ArmEvidence> ArmToEvidence = new HashMap<>();

    PerArmVafConsistencyChecker(final VafPredicate classifier, final ChrArmLocator chrArmLocator)
    {
        Classifier = classifier;
        mChrArmLocator = chrArmLocator;
    }

    void offer(TumorContamination contamination)
    {
        ChrArm chrArm = mChrArmLocator.map(contamination.Chromosome, contamination.Position);
        boolean isEvidenceOfContamination = Classifier.test(contamination);
        ArmToEvidence.computeIfAbsent(chrArm, ArmEvidence::new).register(isEvidenceOfContamination);
    }

    VafConsistencyCheckResult gini()
    {
        Set<ArmEvidence> armEvidenceValues = new HashSet<>(ArmToEvidence.values());
        ArmEvidenceIntegral integral = new ArmEvidenceIntegral(armEvidenceValues);
        double integralValue = integral.value();
        final int totalHits = integral.totalHits();
        final int totalPoints = integral.totalPoints();
        double possibleMax = 0.5 * totalHits * totalPoints;
        double gini = 1.0 - (integralValue / possibleMax);
        return new VafConsistencyCheckResult(gini, totalHits, totalPoints);
    }
}
