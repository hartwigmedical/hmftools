package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._7;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.EnsemblMini;

public class TargetRegionsTestBase
{
    final RefGenomeVersion refGenomeVersion = V38;
    final String chr1 = V38.versionedChromosome(_1.toString());
    final String chr7 = V38.versionedChromosome(_7.toString());
    private final String indelsFilePath = Resources.getResource("reference/msi_indels_1.tsv").getPath();
    private final String ratiosFilePath = Resources.getResource("reference/target_regions_ratios_1.tsv").getPath();
    private final String panelFilePath = Resources.getResource("panel/panel1.bed").getPath();
    EnsemblDataCache ensemblDataCache = EnsemblMini.ensemblMiniDataCache();
    TargetRegionsData targetRegionsData;

    public void setup()
    {
        targetRegionsData = new TargetRegionsData(ratiosFilePath, indelsFilePath);
        targetRegionsData.loadTargetRegionsBed(panelFilePath, ensemblDataCache);
    }
}
