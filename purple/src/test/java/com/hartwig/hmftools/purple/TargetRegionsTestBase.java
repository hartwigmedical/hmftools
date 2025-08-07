package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._7;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.nio.file.Paths;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.purple.SegmentSupport;

public class TargetRegionsTestBase
{
    final RefGenomeVersion refGenomeVersion = V38;
    final Chromosome chromosome1 = _1;
    final Chromosome chromosome7 = _7;
    final String chr1 = V38.versionedChromosome(chromosome1.toString());
    final String chr7 = V38.versionedChromosome(chromosome7.toString());
    private final String ensemblPath = resourceDirectoryPath("ensembl");
    private final String indelsFilePath = referenceFilePath("msi_indels_1.tsv");
    private final String ratiosFilePath = referenceFilePath("target_regions_ratios_1.tsv");
    EnsemblDataCache ensemblDataCache;
    TargetRegionsData targetRegionsData;

    public void setup()
    {
        ensemblDataCache = new EnsemblDataCache(ensemblPath, RefGenomeVersion.V38);
        ensemblDataCache.setRequiredData(true, true, true, false);
        ensemblDataCache.load(false);
        ensemblDataCache.createTranscriptIdMap();

        targetRegionsData = new TargetRegionsData(ratiosFilePath, indelsFilePath);
    }

    String panelFilePath(String fileName)
    {
        return resourceFilePath("panel", fileName);
    }

    String cobaltFilePath(String fileName)
    {
        return resourceFilePath("cobalt", fileName);
    }

    String purpleFilePath(String fileName)
    {
        return resourceFilePath("purple", fileName);
    }

    String referenceFilePath(String fileName)
    {
        return resourceFilePath("reference", fileName);
    }

    private String resourceFilePath(String directory, String fileName)
    {
        return Paths.get("src", "test", "resources", directory, fileName).toAbsolutePath().toString();
    }

    private String resourceDirectoryPath(String directory)
    {
        return Paths.get("src", "test", "resources", directory).toAbsolutePath().toString();
    }

    protected PurpleSegment ps(String chromosome, int start, int end, GermlineStatus status)
    {
        return new PurpleSegment(chromosome, start, end, false, SegmentSupport.TELOMERE,
                123, 0.33, 123, 0.34, 0.45, 0.34,
                status, true, 0.44, 123, 324, 0.45, 0.45,
                0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45);
    }
}
