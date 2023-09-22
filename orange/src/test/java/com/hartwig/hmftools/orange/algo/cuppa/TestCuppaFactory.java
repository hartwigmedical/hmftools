package com.hartwig.hmftools.orange.algo.cuppa;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
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
}
