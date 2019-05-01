package com.hartwig.hmftools.patientreporter.variants.germline;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GermlineGenesReporting {

    @NotNull
    public abstract Map<String, Boolean> germlineGenes();

    @NotNull
    public Set<String> reportableGermlineGenes() {
        return germlineGenes().keySet();
    }

    public boolean genesToNotifyClinicalGeneticist(@NotNull String germlineVariantGene) {
        return germlineGenes().get(germlineVariantGene);
    }
}
