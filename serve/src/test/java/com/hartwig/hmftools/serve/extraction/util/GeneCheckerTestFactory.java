package com.hartwig.hmftools.serve.extraction.util;

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
