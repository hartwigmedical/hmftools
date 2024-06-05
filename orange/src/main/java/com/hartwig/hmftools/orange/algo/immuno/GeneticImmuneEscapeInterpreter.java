package com.hartwig.hmftools.orange.algo.immuno;

import com.hartwig.hmftools.datamodel.immuno.GeneticImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.immuno.ImmutableGeneticImmuneEscapeRecord;

import org.jetbrains.annotations.NotNull;

public class GeneticImmuneEscapeInterpreter
{
    @NotNull
    public GeneticImmuneEscapeRecord interpret() {
        return ImmutableGeneticImmuneEscapeRecord.builder()
                .hasHlaEscapePresent(false)
                .hasAntigenPresentationPathwayEscape(false)
                .hasIFNGammaInactivationEscape(false)
                .hasPDL1OverexpressionEscape(false)
                .hasCD58InactivationEscape(false)
                .hasEpigeneticSETDB1Escape(false)
                .build();
    }
}
