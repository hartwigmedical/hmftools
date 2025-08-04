package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._7;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.nio.file.Paths;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class TargetRegionsTestBase
{
    final RefGenomeVersion refGenomeVersion = V38;
    final String chr1 = V38.versionedChromosome(_1.toString());
    final String chr7 = V38.versionedChromosome(_7.toString());
    private final String ensemblPath = Paths.get("src", "test", "resources", "ensembl").toAbsolutePath().toString();
    private final String indelsFilePath =
            Paths.get("src", "test", "resources", "reference", "msi_indels_1.tsv").toAbsolutePath().toString();
    private final String ratiosFilePath =
            Paths.get("src", "test", "resources", "reference", "target_regions_ratios_1.tsv").toAbsolutePath().toString();
    private final String panelFilePath =
            Paths.get("src", "test", "resources", "panel", "panel_1.bed").toAbsolutePath().toString();
    EnsemblDataCache ensemblDataCache;
    TargetRegionsData targetRegionsData;

    public void setup()
    {
        ensemblDataCache = new EnsemblDataCache(ensemblPath, RefGenomeVersion.V38);
        ensemblDataCache.setRequiredData(true, true, true, false);
        ensemblDataCache.load(false);
        ensemblDataCache.createTranscriptIdMap();

        targetRegionsData = new TargetRegionsData(ratiosFilePath, indelsFilePath);
        targetRegionsData.loadTargetRegionsBed(panelFilePath, ensemblDataCache);
    }
}
