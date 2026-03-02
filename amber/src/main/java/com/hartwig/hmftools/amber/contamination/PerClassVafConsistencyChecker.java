package com.hartwig.hmftools.amber.contamination;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.amber.VafReading;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class PerClassVafConsistencyChecker<S, T extends Comparable<T>>
{
    public static VafConsistencyCheckResult<ChrArm> calculateConfirmationFactor(final ChrArmLocator chrArmLocator,
            double vaf,
            List<VafReading> data)
    {
        VafPredicate<VafReading> classifier = new BinomialVafPredicate(vaf);
        VafClassifier<VafReading, ChrArm> chrArmClassifier = VafClassifier.chrArmClassifier(chrArmLocator);
        PerClassVafConsistencyChecker<VafReading, ChrArm> confirmation = new PerClassVafConsistencyChecker<>(classifier, chrArmClassifier);
        data.forEach(confirmation::offer);
        return confirmation.unevenDistributionCost();
    }

    private final VafPredicate<S> Classifier;
    private final VafClassifier<S, T> mChrArmLocator;
    private final Map<T, CategoryEvidence<T>> ArmToEvidence = new HashMap<>();

    public PerClassVafConsistencyChecker(final VafPredicate<S> classifier, final VafClassifier<S, T> chrArmLocator)
    {
        Classifier = classifier;
        mChrArmLocator = chrArmLocator;
    }

    public void offer(S contamination)
    {
        T chrArm = mChrArmLocator.classify(contamination);
        Preconditions.checkArgument(chrArm != null);
        boolean isEvidenceOfContamination = Classifier.test(contamination);
        ArmToEvidence.computeIfAbsent(chrArm, CategoryEvidence::new).register(isEvidenceOfContamination);
    }

    public VafConsistencyCheckResult<T> unevenDistributionCost()
    {
        Set<CategoryEvidence<T>> categoryEvidenceValues = new HashSet<>(ArmToEvidence.values());
        CategoryEvidenceIntegral<T> integral = new CategoryEvidenceIntegral<>(categoryEvidenceValues);
        double integralValue = integral.value();
        final int totalHits = integral.totalHits();
        final int totalPoints = integral.totalPoints();
        double possibleMax = 0.5 * totalHits * totalPoints;
        double cost = integralValue / possibleMax;
        return new VafConsistencyCheckResult<>(cost, totalHits, totalPoints, integral.mCategoryEvidence);
    }
}
