package com.hartwig.hmftools.linx.annotators;

import static java.lang.String.format;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxConfig.REF_GENOME_VERSION;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.loadConfigFile;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.jetbrains.annotations.Nullable;

public class FragileSiteAnnotator
{
    private List<ChrBaseRegion> mFragileSites;

    public FragileSiteAnnotator()
    {
        mFragileSites = Lists.newArrayList();
    }

    public static String fragileSitesResourceFile(final RefGenomeVersion refGenomeVersion)
    {
        return format("/sites/fragile_sites.%s.csv", refGenomeVersion.identifier());
    }

    public void loadFragileSitesFile(@Nullable final String filename)
    {
        try
        {
            List<String> fileContents = Lists.newArrayList();

            if(filename != null)
            {
                fileContents.addAll(Files.readAllLines(new File(filename).toPath()));
            }
            else
            {
                fileContents.addAll(new BufferedReader(new InputStreamReader(
                        FragileSiteAnnotator.class.getResourceAsStream(fragileSitesResourceFile(REF_GENOME_VERSION))))
                        .lines().collect(Collectors.toList()));
            }

            mFragileSites.addAll(loadConfigFile(fileContents, REF_GENOME_VERSION));
            LNX_LOGGER.info("loaded {} known fragile sites", mFragileSites.size());
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to read fragile site file: {}", e.toString());
        }
    }

    public boolean isFragileSite(final SvVarData svData, final boolean useStart)
    {
        if(mFragileSites.isEmpty())
            return false;

        for(final ChrBaseRegion fsRegion : mFragileSites)
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
