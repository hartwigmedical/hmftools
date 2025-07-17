package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.CANDIDATE_PROBES_FILE_NAME;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.REJECTED_REGIONS_FILE_STEM;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.TARGET_REGIONS_FILE_NAME;

import java.io.IOException;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PanelBuilder
{
    private final PanelBuilderConfig mConfig;
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeGenerator mProbeGenerator;
    // TODO: is this the best way to scope this?
    private final CoveredRegions mCoveredRegions;
    @Nullable
    private OutputWriter mOutputWriter;

    private static final Logger LOGGER = LogManager.getLogger(PanelBuilder.class);

    public PanelBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelBuilderConfig(configBuilder);
        mRefGenome = loadRefGenome(mConfig.RefGenomeFile);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        ProbeQualityProfile probeQualityProfile = new ProbeQualityProfile(mConfig.ProbeQualityProfileFile);
        ProbeEvaluator probeEvaluator = new ProbeEvaluator(mRefGenome, probeQualityProfile, this::writeCandidateProbe);
        ProbeSelector probeSelector = new ProbeSelector(probeEvaluator);
        CandidateProbeGenerator candidateGenerator = new CandidateProbeGenerator(mRefGenome.chromosomeLengths());
        mProbeGenerator = new ProbeGenerator(candidateGenerator, probeSelector);
        mCoveredRegions = new CoveredRegions();
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting panel builder");

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.OutputDir);
        mOutputWriter = new OutputWriter(
                mConfig.outputFilePath(PANEL_PROBES_FILE_STEM),
                mConfig.outputFilePath(TARGET_REGIONS_FILE_NAME),
                mConfig.outputFilePath(REJECTED_REGIONS_FILE_STEM),
                mConfig.VerboseOutput ? mConfig.outputFilePath(CANDIDATE_PROBES_FILE_NAME) : null);

        LOGGER.info("Generating probes");
        ProbeGenerationResult geneProbes = generateTargetGeneProbes();
        ProbeGenerationResult cnBackboneProbes = generateCopyNumberBackboneProbes();
        ProbeGenerationResult customRegionProbes = generateCustomRegionProbes();
        LOGGER.info("Probe generation done");

        LOGGER.info("Writing output");
        {
            ProbeGenerationResult aggregate = customRegionProbes.add(geneProbes).add(cnBackboneProbes);
            mOutputWriter.writePanelProbes(aggregate.probes());
            mOutputWriter.writeTargetRegions(aggregate.targetRegions());
            mOutputWriter.writeRejectedRegions(aggregate.rejectedRegions());
        }

        // TODO: strategy for handling overlapping probes

        mOutputWriter.close();
        mOutputWriter = null;

        mCoveredRegions.clear();

        LOGGER.info("Panel builder complete, mins({})", runTimeMinsStr(startTimeMs));
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
            ProbeGenerationResult result = TargetGenes.generateProbes(mConfig.TargetGenesFile, ensemblData, mProbeGenerator);
            mCoveredRegions.addFromProbes(result.probes());
            return result;
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
            ProbeGenerationResult result =
                    CopyNumberBackbone.generateProbes(mConfig.AmberSitesFile, mRefGenomeVersion, mProbeGenerator, mCoveredRegions);
            mCoveredRegions.addFromProbes(result.probes());
            return result;
        }
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
            ProbeGenerationResult result =
                    CustomRegions.generateProbes(mConfig.CustomRegionsFile, mRefGenome.chromosomeLengths(), mProbeGenerator);
            mCoveredRegions.addFromProbes(result.probes());
            return result;
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
        mOutputWriter.writeCandidateProbe(probe);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PanelBuilderConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        try
        {
            PanelBuilder panelBuilder = new PanelBuilder(configBuilder);
            panelBuilder.run();
        }
        catch(IOException e)
        {
            LOGGER.error("IO error", e);
            exit(1);
        }
        catch(UserInputError e)
        {
            LOGGER.error("Bad input data");
            LOGGER.error(e.getMessage());
            exit(1);
        }
    }
}
