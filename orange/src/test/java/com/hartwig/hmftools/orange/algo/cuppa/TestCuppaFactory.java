package com.hartwig.hmftools.orange.algo.cuppa;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaMode;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;

import org.jetbrains.annotations.NotNull;

public final class TestCuppaFactory
{
    @NotNull
    public static CuppaData createMinimalCuppaData()
    {
        // Some downstream algo's expect at least one prediction, so that is considered "minimal"
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder()
                .cancerType("cancer")
                .likelihood(1D)
                .build();

        List<CuppaPrediction> predictions = Lists.newArrayList();
        predictions.add(prediction);

        return ImmutableCuppaData.builder()
                .bestPrediction(prediction)
                .predictions(predictions)
                .simpleDups32To200B(0)
                .maxComplexSize(0)
                .telomericSGLs(0)
                .cuppaMode(CuppaMode.WGTS)
                .lineCount(0)
                .build();
    }
}
