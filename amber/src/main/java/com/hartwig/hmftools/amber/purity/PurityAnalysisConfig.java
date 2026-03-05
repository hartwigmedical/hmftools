package com.hartwig.hmftools.amber.purity;

import com.hartwig.hmftools.amber.AmberConfig;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public record PurityAnalysisConfig(String tumorId, RefGenomeVersion refGenomeVersion, String outputDir, int threads)
{
    public PurityAnalysisConfig(AmberConfig config)
    {
        this(config.TumorId, config.RefGenVersion, config.OutputDir, config.Threads);
    }
}
