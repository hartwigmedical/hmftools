package com.hartwig.hmftools.amber.contamination;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.amber.VafReading;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class PerArmVafConsistencyChecker
{
    public static VafConsistencyCheckResult<ChrArm> calculateConfirmationFactor(final ChrArmLocator chrArmLocator,
            double vaf,
            List<VafReading> data)
    {
        VafPredicate classifier = new BinomialVafPredicate(vaf);
        PerArmVafConsistencyChecker confirmation = new PerArmVafConsistencyChecker(classifier, chrArmLocator);
        data.forEach(confirmation::offer);
        return confirmation.unevenDistributionCost();
    }

    private final VafPredicate Classifier;
    private final ChrArmLocator mChrArmLocator;
    private final Map<ChrArm, CategoryEvidence<ChrArm>> ArmToEvidence = new HashMap<>();

    public PerArmVafConsistencyChecker(final VafPredicate classifier, final ChrArmLocator chrArmLocator)
    {
        Classifier = classifier;
        mChrArmLocator = chrArmLocator;
    }

    public void offer(VafReading contamination)
    {
        ChrArm chrArm = mChrArmLocator.map(contamination.chromosome(), contamination.position());
        boolean isEvidenceOfContamination = Classifier.test(contamination);
        ArmToEvidence.computeIfAbsent(chrArm, CategoryEvidence::new).register(isEvidenceOfContamination);
    }

    public VafConsistencyCheckResult<ChrArm> unevenDistributionCost()
    {
        Set<CategoryEvidence<ChrArm>> categoryEvidenceValues = new HashSet<>(ArmToEvidence.values());
        CategoryEvidenceIntegral<ChrArm> integral = new CategoryEvidenceIntegral<>(categoryEvidenceValues);
        double integralValue = integral.value();
        final int totalHits = integral.totalHits();
        final int totalPoints = integral.totalPoints();
        double possibleMax = 0.5 * totalHits * totalPoints;
        double cost = integralValue / possibleMax;
        return new VafConsistencyCheckResult<>(cost, totalHits, totalPoints, integral.mCategoryEvidence);
    }
}
