package com.hartwig.hmftools.serve.sources.vicc.check;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.jetbrains.annotations.NotNull;

public final class GeneCheckerTestFactory {

    private GeneCheckerTestFactory() {
    }

    @NotNull
    public static GeneChecker buildForHG19() {
        return new GeneChecker(HmfGenePanelSupplier.allGenesMap37().keySet());
    }
}
