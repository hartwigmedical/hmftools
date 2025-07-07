package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.PANEL_DEFINITION_FILE_EXTENSION;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.writePanelDefinition;

import java.io.FileNotFoundException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class PanelBuilder
{
    private final PanelConfig mConfig;

    private final PanelCache mPanelCache;

    private final ProbeQualityProfile mProbeQualityProfile;

    public PanelBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelConfig(configBuilder);

        mPanelCache = new PanelCache();

        mProbeQualityProfile = new ProbeQualityProfile(mConfig.ProbeQualityProfileFile);
    }

    public void run()
    {
        GU_LOGGER.info("starting panel builder");

        long startTimeMs = System.currentTimeMillis();

        // first load custom regions
        CustomRegions customRegions = new CustomRegions(mConfig, mPanelCache);
        customRegions.run();

        // generate gene probes
        GeneProbesGenerator geneProbesGenerator = new GeneProbesGenerator(mConfig, mPanelCache, mProbeQualityProfile);
        geneProbesGenerator.run();

        // build copy-number backbone to fill in gaps
        CopyNumberBackbone copyNumberBackbone = new CopyNumberBackbone(mConfig, mPanelCache, mProbeQualityProfile);
        copyNumberBackbone.run();

        // write final regions and probes
        writeFinalPanelRegions();

        GU_LOGGER.info("panel builder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void writeFinalPanelRegions()
    {
        List<PanelRegion> panelRegions = Lists.newArrayList();

        mPanelCache.chrRegionsMap().values().forEach(panelRegions::addAll);

        // sort and merge - for now keep the first region's source info and type
        Collections.sort(panelRegions);

        // merge any adjacent regions
        int index = 0;
        while(index < panelRegions.size() - 1)
        {
            PanelRegion region = panelRegions.get(index);
            PanelRegion nextRegion = panelRegions.get(index + 1);

            if(region.Chromosome.equals(nextRegion.Chromosome) && region.end() >= nextRegion.start() - 2)
            {
                ChrBaseRegion newRegion = new ChrBaseRegion(
                        region.Chromosome, min(region.start(), nextRegion.start()), max(region.end(), nextRegion.end()));

                String sourceInfo = format("%s;%s", region.SourceInfo, nextRegion.SourceInfo);

                PanelRegion newPanelRegion = new PanelRegion(newRegion, RegionType.MIXED, sourceInfo);

                panelRegions.set(index, newPanelRegion);
                panelRegions.remove(index + 1);
            }
            else
            {
                ++index;
            }
        }

        String panelDefinitionFilename = mConfig.formOutputFilename(PANEL_DEFINITION_FILE_EXTENSION);
        writePanelDefinition(panelDefinitionFilename, panelRegions);
    }

    public static void main(@NotNull final String[] args) throws FileNotFoundException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PanelConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PanelBuilder panelBuilder = new PanelBuilder(configBuilder);
        panelBuilder.run();
    }
}
