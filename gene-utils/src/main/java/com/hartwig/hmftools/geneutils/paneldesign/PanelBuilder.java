package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CANDIDATE_PROBES_FILE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PANEL_PROBES_FILE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.REJECTED_REGIONS_FILE;

import java.util.List;

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
    private final ProbeGenerator mProbeGenerator;
    private OutputWriter mOutputWriter;

    private static final Logger LOGGER = LogManager.getLogger(PanelBuilder.class);

    public PanelBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelBuilderConfig(configBuilder);
        RefGenomeSource mRefGenome = loadRefGenome(mConfig.RefGenomeFile);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        ProbeQualityProfile probeQualityProfile = new ProbeQualityProfile(mConfig.ProbeQualityProfileFile);
        ProbeEvaluator probeEvaluator = new ProbeEvaluator(mRefGenome, probeQualityProfile, this::writeCandidateProbe);
        mProbeGenerator = new ProbeGenerator(probeEvaluator);
    }

    public void run()
    {
        LOGGER.info("Starting panel builder");

        long startTimeMs = System.currentTimeMillis();

        mOutputWriter = new OutputWriter(
                mConfig.outputFilePath(PANEL_PROBES_FILE),
                mConfig.outputFilePath(REJECTED_REGIONS_FILE),
                mConfig.outputFilePath(CANDIDATE_PROBES_FILE));

        LOGGER.info("Generating probes");
        ProbeGenerationResult customRegionProbes = generateCustomRegionProbes();
        ProbeGenerationResult geneProbes = generateTargetGeneProbes();
        ProbeGenerationResult cnBackboneProbes = generateCopyNumberBackboneProbes();
        LOGGER.info("Probe generation done");

        LOGGER.debug("Writing output");
        {
            checkCreateOutputDir(mConfig.OutputDir);
            ProbeGenerationResult aggregate = customRegionProbes.add(geneProbes).add(cnBackboneProbes);
            mOutputWriter.writePanelProbes(aggregate.probes());
            mOutputWriter.writeRejectedRegions(aggregate.rejectedRegions());
        }

        // TODO: profile code (it's kinda slow)

        // TODO: other output to write?

        // TODO: remove duplicate/overlapping probes?

        // TODO: check probes are within chromosome bounds?

        mOutputWriter.close();
        mOutputWriter = null;

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
            return CustomRegions.generateProbes(mConfig.CustomRegionsFile, mProbeGenerator);
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
            EnsemblDataCache ensemblData = loadEnsemblData();
            return TargetGenes.generateProbes(mConfig.TargetGenesFile, ensemblData, mProbeGenerator);
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
            return CopyNumberBackbone.generateProbes(mConfig.AmberSitesFile, mRefGenomeVersion, mProbeGenerator);
        }
    }

    private EnsemblDataCache loadEnsemblData()
    {
        EnsemblDataCache ensemblData = new EnsemblDataCache(mConfig.EnsemblDir, mRefGenomeVersion);
        ensemblData.setRequiredData(true, false, false, false);
        ensemblData.load(false);
        return ensemblData;
    }

    private void writeCandidateProbe(final EvaluatedProbe probe)
    {
        mOutputWriter.writeCandidateProbes(List.of(probe));
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
