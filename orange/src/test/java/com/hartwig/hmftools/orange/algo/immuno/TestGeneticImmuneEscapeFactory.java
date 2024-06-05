package com.hartwig.hmftools.orange.algo.immuno;

import com.hartwig.hmftools.datamodel.immuno.ImmutableGeneticImmuneEscapeRecord;

import org.jetbrains.annotations.NotNull;

public final class TestGeneticImmuneEscapeFactory
{
    @NotNull
    public static ImmutableGeneticImmuneEscapeRecord.Builder builder() {
        return ImmutableGeneticImmuneEscapeRecord.builder()
                .hasHlaEscapePresent(false)
                .hasAntigenPresentationPathwayEscape(false)
                .hasIFNGammaInactivationEscape(false)
                .hasPDL1OverexpressionEscape(false)
                .hasCD58InactivationEscape(false)
                .hasEpigeneticSETDB1Escape(false);
    }
}
