package com.hartwig.hmftools.orange.isofox;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxData;

import org.jetbrains.annotations.NotNull;

public final class IsofoxInterpreter {

    private IsofoxInterpreter() {
    }

    @NotNull
    public static IsofoxInterpretedData interpret(@NotNull IsofoxData isofox, @NotNull List<DriverGene> driverGenes,
            @NotNull KnownFusionCache knownFusionCache) {
        return ImmutableIsofoxInterpretedData.builder()
                .summary(isofox.summary())
                .allGeneExpressions(isofox.geneExpressions())
                .allFusions(isofox.fusions())
                .allNovelSpliceJunctions(isofox.novelSpliceJunctions())
                .build();
    }
}
