package com.hartwig.hmftools.linx.annotators;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.RG_VERSION;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.linx.types.SvVarData;

public class FragileSiteAnnotator
{
    private List<BaseRegion> mFragileSites;

    private static final int CSV_REQUIRED_FIELDS = 3;

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

                final BaseRegion genomeRegion = new BaseRegion(
                        RG_VERSION.versionedChromosome(items[FS_COL_CHR]),
                        Integer.parseInt(items[FS_COL_POS_START]),
                        Integer.parseInt(items[FS_COL_POS_END]));

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

        for(final BaseRegion fsRegion : mFragileSites)
        {
            if(fsRegion.containsPosition(svData.chromosome(useStart), svData.position(useStart)))
            {
                LNX_LOGGER.debug("var({}) found in known fragile site({})", svData.posId(), fsRegion);
                return true;
            }
        }

        return false;
    }

}
