package com.hartwig.hmftools.orange.algo.cuppa;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa2.Cuppa2Data;
import com.hartwig.hmftools.datamodel.cuppa2.FeatureContributionEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableCuppa2Data;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableFeatureContributionEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableProbabilityEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ImmutableSignatureQuantileEntry;
import com.hartwig.hmftools.datamodel.cuppa2.ProbabilityEntry;
import com.hartwig.hmftools.datamodel.cuppa2.SignatureQuantileEntry;

import org.jetbrains.annotations.NotNull;

public final class TestCuppaFactory
{
    @NotNull
    public static CuppaData createMinimalCuppaData()
    {
        // Some downstream algo's expect at least one prediction, so that is considered "minimal"
        List<CuppaPrediction> predictions = Lists.newArrayList();
        predictions.add(ImmutableCuppaPrediction.builder().cancerType("cancer").likelihood(1D).build());

        return ImmutableCuppaData.builder()
                .predictions(predictions)
                .simpleDups32To200B(0)
                .maxComplexSize(0)
                .telomericSGLs(0)
                .lineCount(0)
                .build();
    }

    @NotNull
    public static Cuppa2Data createMinimalCuppa2Data()
    {
        ProbabilityEntry minimalTopPrediction = ImmutableProbabilityEntry.builder().
                classifierName("dummy_classifier_name").
                cancerType("dummy_cancer_type_1").
                probability(0.95).
                rank(1).
                rankGroup(0).
                build();

        List<ProbabilityEntry> minimalProbs = Arrays.asList(minimalTopPrediction);

        FeatureContributionEntry minimalFeatureContributionEntry = ImmutableFeatureContributionEntry.builder().
                classifierName("dummy_classifier_name").
                featureName("dummy_feat_name").
                featureValue(0.0).
                cancerType("dummy_cancer_type_1").
                featureContribution(0.0).
                rank(1).
                rankGroup(1).
                build();

        List<FeatureContributionEntry> minimalFeatContribs = Arrays.asList(minimalFeatureContributionEntry);

        SignatureQuantileEntry minimalSignatureQuantileEntry = ImmutableSignatureQuantileEntry.builder().
                signatureName("dummy_feat_name").
                signatureCount(0.0).
                cancerType("dummy_cancer_type_1").
                signatureQuantile(0.0).
                rank(1).
                rankGroup(2).
                build();

        List<SignatureQuantileEntry> minimalSigQuantiles = Arrays.asList(minimalSignatureQuantileEntry);

        return ImmutableCuppa2Data.builder().
                topPrediction(minimalTopPrediction).
                probabilities(minimalProbs).
                featureContributions(minimalFeatContribs).
                signatureQuantiles(minimalSigQuantiles).
                build();
    }
}
