package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.System.exit;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;

import java.io.IOException;
import java.util.Optional;

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
    private PanelData mPanelData;
    @Nullable
    private OutputWriter mOutputWriter;

    private static final Logger LOGGER = LogManager.getLogger(PanelBuilder.class);

    public PanelBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelBuilderConfig(configBuilder);
        mRefGenome = loadRefGenome(mConfig.RefGenomeFile);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        ProbeQualityProfile probeQualityProfile = ProbeQualityProfile.loadFromResourceFile(mConfig.ProbeQualityProfileFile);
        ProbeEvaluator probeEvaluator = new ProbeEvaluator(mRefGenome, probeQualityProfile, this::writeCandidateProbe);
        CandidateProbeGenerator candidateGenerator = new CandidateProbeGenerator(mRefGenome.chromosomeLengths());
        mProbeGenerator = new ProbeGenerator(candidateGenerator, probeEvaluator);
        mPanelData = new PanelData();
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting panel builder");

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.OutputDir);
        mOutputWriter = new OutputWriter(mConfig.OutputDir, mConfig.OutputPrefix, mConfig.VerboseOutput);

        LOGGER.info("Generating probes");
        mPanelData = new PanelData();
        // Note the order of generation here determines the priority of probe overlap resolution.
        // Probes generated first will exclude overlapping probes generated afterward.
        Optional<TargetGenes.Stats> geneStats = generateTargetGeneProbes();
        generateCustomRegionProbes();
        generateCopyNumberBackboneProbes();

        LOGGER.info("Writing output");
        {
            mOutputWriter.writePanelProbes(mPanelData.probes());
            mOutputWriter.writeTargetRegions(mPanelData.coveredTargetRegions());
            mOutputWriter.writeCandidateRegions(mPanelData.candidateTargetRegions());
            mOutputWriter.writeRejectedRegions(mPanelData.rejectedRegions());
            geneStats.ifPresent(stats -> mOutputWriter.writeGeneStats(stats.perGene()));
        }

        mOutputWriter.close();
        mOutputWriter = null;

        LOGGER.info("Panel builder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private Optional<TargetGenes.Stats> generateTargetGeneProbes()
    {
        if(mConfig.TargetGenesFile == null)
        {
            LOGGER.info("Target genes not provided; skipping gene probes");
            return Optional.empty();
        }
        else
        {
            EnsemblDataCache ensemblData = loadEnsemblData();
            TargetGenes.Stats stats = TargetGenes.generateProbes(mConfig.TargetGenesFile, ensemblData, mProbeGenerator, mPanelData);
            // Result is stored into mPanelData
            return Optional.of(stats);
        }
    }

    private void generateCopyNumberBackboneProbes()
    {
        if(mConfig.AmberSitesFile == null)
        {
            LOGGER.info("Amber sites not provided; skipping copy number backbone probes");
        }
        else
        {
            CopyNumberBackbone copyNumberBackbone =
                    new CopyNumberBackbone(mConfig.AmberSitesFile, mRefGenomeVersion, mProbeGenerator, mPanelData);
            copyNumberBackbone.generateProbes();
            // Result is stored into mPanelData
        }
    }

    private void generateCustomRegionProbes()
    {
        if(mConfig.CustomRegionsFile == null)
        {
            LOGGER.info("Custom regions not provided; skipping custom region probes");
        }
        else
        {
            CustomRegions.generateProbes(mConfig.CustomRegionsFile, mRefGenome.chromosomeLengths(), mProbeGenerator, mPanelData);
            // Result is stored into mPanelData
        }
    }

    private EnsemblDataCache loadEnsemblData()
    {
        EnsemblDataCache ensemblData = new EnsemblDataCache(mConfig.EnsemblDir, mRefGenomeVersion);
        ensemblData.setRequiredData(true, false, false, false);
        ensemblData.load(false);
        return ensemblData;
    }

    private void writeCandidateProbe(final Probe probe)
    {
        requireNonNull(mOutputWriter).writeCandidateProbe(probe);
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
