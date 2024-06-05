package com.hartwig.hmftools.orange.algo.immuno;

import com.hartwig.hmftools.datamodel.immuno.ImmutableImmuneEscapeRecord;

import org.jetbrains.annotations.NotNull;

public final class TestImmuneEscapeFactory
{
    @NotNull
    public static ImmutableImmuneEscapeRecord.Builder builder() {
        return ImmutableImmuneEscapeRecord.builder()
                .hasHlaEscape(false)
                .hasAntigenPresentationPathwayEscape(false)
                .hasIFNGammaPathwayEscape(false)
                .hasPDL1OverexpressionEscape(false)
                .hasCD58InactivationEscape(false)
                .hasEpigeneticSETDB1Escape(false);
    }
}
