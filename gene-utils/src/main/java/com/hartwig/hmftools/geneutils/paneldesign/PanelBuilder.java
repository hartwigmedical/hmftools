package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.writePanelProbes;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.writeRejectedRegions;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PANEL_PROBES_FILE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_GC_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_GC_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_QUALITY_SCORE_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.REJECTED_REGIONS_FILE;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PanelBuilder
{
    private final PanelBuilderConfig mConfig;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeEvaluator mProbeEvaluator;

    private static final Logger LOGGER = LogManager.getLogger(PanelBuilder.class);

    public PanelBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelBuilderConfig(configBuilder);
        RefGenomeSource mRefGenome = loadRefGenome(mConfig.RefGenomeFile);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        ProbeQualityProfile probeQualityProfile = new ProbeQualityProfile(mConfig.ProbeQualityProfileFile);
        mProbeEvaluator = new ProbeEvaluator(mRefGenome, probeQualityProfile, PROBE_QUALITY_SCORE_MIN, PROBE_GC_MIN, PROBE_GC_MAX);
    }

    public void run()
    {
        LOGGER.info("Starting panel builder");

        long startTimeMs = System.currentTimeMillis();

        LOGGER.info("Starting probe generation");
        ProbeGenerationResult customRegionProbes = generateCustomRegionProbes();
        ProbeGenerationResult geneProbes = generateTargetGeneProbes();
        ProbeGenerationResult cnBackboneProbes = generateCopyNumberBackboneProbes();
        LOGGER.info("Probe generation done");

        LOGGER.debug("Writing output");
        {
            ProbeGenerationResult aggregate = customRegionProbes.add(geneProbes).add(cnBackboneProbes);
            writePanelProbes(mConfig.outputFilePath(PANEL_PROBES_FILE), aggregate.probes().stream());
            writeRejectedRegions(mConfig.outputFilePath(REJECTED_REGIONS_FILE), aggregate.rejectedRegions().stream());
        }

        // TODO: other output to write?

        // TODO: remove duplicate/overlapping probes?

        LOGGER.info("Panel builder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private ProbeGenerationResult generateCustomRegionProbes()
    {
        if(mConfig.CustomRegionsFile == null)
        {
            LOGGER.info("Custom regions not provided; skipping custom region probes");
            return new ProbeGenerationResult();
        }
        else
        {
            return CustomRegions.generateProbes(mConfig.CustomRegionsFile, mProbeEvaluator);
        }
    }

    private ProbeGenerationResult generateTargetGeneProbes()
    {
        if(mConfig.TargetGenesFile == null)
        {
            LOGGER.info("Target genes not provided; skipping gene probes");
            return new ProbeGenerationResult();
        }
        else
        {
            EnsemblDataCache ensemblData = new EnsemblDataCache(mConfig.EnsemblDir, mRefGenomeVersion);
            return TargetGenes.generateProbes(mConfig.TargetGenesFile, ensemblData, mProbeEvaluator);
        }
    }

    private ProbeGenerationResult generateCopyNumberBackboneProbes()
    {
        if(mConfig.AmberSitesFile == null)
        {
            LOGGER.info("Amber sites not provided; skipping copy number backbone probes");
            return new ProbeGenerationResult();
        }
        else
        {
            return CopyNumberBackbone.generateProbes(mConfig.AmberSitesFile, mRefGenomeVersion, mProbeEvaluator);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PanelBuilderConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PanelBuilder panelBuilder = new PanelBuilder(configBuilder);
        panelBuilder.run();
    }
}
