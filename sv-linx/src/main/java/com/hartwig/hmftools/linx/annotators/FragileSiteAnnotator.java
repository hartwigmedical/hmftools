package com.hartwig.hmftools.linx.annotators;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.refGenomeChromosome;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.linx.types.SvVarData;

public class FragileSiteAnnotator
{
    private List<GenomeRegion> mFragileSites;

    private final static int CSV_REQUIRED_FIELDS = 3;

    public FragileSiteAnnotator()
    {
        mFragileSites = Lists.newArrayList();
    }

    private static final int FS_COL_CHR = 0;
    private static final int FS_COL_POS_START = 1;
    private static final int FS_COL_POS_END = 2;

    public void loadFragileSitesFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            for(final String line : fileContents)
            {
                if(line.contains("Chromosome"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < CSV_REQUIRED_FIELDS)
                    continue;

                final GenomeRegion genomeRegion = GenomeRegions.create(
                        refGenomeChromosome(items[FS_COL_CHR], RG_VERSION),
                        Long.parseLong(items[FS_COL_POS_START]),
                        Long.parseLong(items[FS_COL_POS_END]));

                mFragileSites.add(genomeRegion);
            }

            LNX_LOGGER.info("loaded {} known fragile sites from file: {}", mFragileSites.size(), filename);
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read fragile site CSV file({})", filename);
        }
    }

    public boolean isFragileSite(final SvVarData svData, final boolean useStart)
    {
        if(mFragileSites.isEmpty())
            return false;

        for(final GenomeRegion genomeRegion : mFragileSites)
        {
            if(genomeRegion.chromosome().equals(svData.chromosome(useStart))
            && genomeRegion.start() <= svData.position(useStart) && genomeRegion.end() >= svData.position(useStart))
            {
                LNX_LOGGER.debug("var({}) found in known fragile site",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return true;
            }
        }

        return false;
    }

}
