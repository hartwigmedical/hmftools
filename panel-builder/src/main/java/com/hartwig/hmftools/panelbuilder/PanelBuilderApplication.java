package com.hartwig.hmftools.panelbuilder;

import static java.lang.System.exit;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.APP_NAME;

import java.io.IOException;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PanelBuilderApplication
{
    private final PanelBuilderConfig mConfig;
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeGenerator mProbeGenerator;
    private PanelData mPanelData;
    @Nullable
    private OutputWriter mOutputWriter;

    private static final Logger LOGGER = LogManager.getLogger(PanelBuilderApplication.class);

    public PanelBuilderApplication(final PanelBuilderConfig config)
    {
        mConfig = config;
        mRefGenome = loadRefGenome(mConfig.refGenomeFile());
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        ProbeQualityProfile probeQualityProfile = ProbeQualityProfile.loadFromResourceFile(mConfig.probeQualityProfileFile());
        mProbeGenerator = ProbeGenerator.create(mRefGenome, probeQualityProfile, this::writeCandidateProbe);
        mPanelData = new PanelData();
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting panel builder");

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());
        mOutputWriter = new OutputWriter(mConfig.outputDir(), mConfig.outputId(), mConfig.verboseOutput());

        LOGGER.info("Generating probes");
        mPanelData = new PanelData();
        // Note the order of generation here determines the priority of probe overlap resolution.
        // Probes generated first will exclude overlapping probes generated afterward.
        Genes.ExtraOutput geneExtraOutput = generateTargetGeneProbes();
        generateCustomRegionProbes();
        generateCopyNumberBackboneProbes();
        SampleVariants.ExtraOutput sampleExtraOutput = generateSampleVariantProbes();

        LOGGER.info("Writing output");
        mOutputWriter.writePanelProbes(mPanelData.probes());
        mOutputWriter.writeTargetRegions(mPanelData.coveredTargetRegions());
        mOutputWriter.writeCandidateRegions(mPanelData.candidateTargetRegions());
        mOutputWriter.writeRejectedRegions(mPanelData.rejectedRegions());
        if(geneExtraOutput != null)
        {
            mOutputWriter.writeGeneStats(geneExtraOutput.geneStats());
        }
        if(sampleExtraOutput != null)
        {
            mOutputWriter.writeSampleVariants(sampleExtraOutput.variants());
        }

        mOutputWriter.close();
        mOutputWriter = null;

        LOGGER.info("Panel builder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    @Nullable
    private Genes.ExtraOutput generateTargetGeneProbes()
    {
        if(mConfig.genesFile() == null)
        {
            LOGGER.info("Genes not provided; skipping gene probes");
            return null;
        }
        else
        {
            EnsemblDataCache ensemblData = loadEnsemblData();
            Genes.ExtraOutput extraOutput =
                    Genes.generateProbes(mConfig.genesFile(), ensemblData, mProbeGenerator, mPanelData);
            // Result is stored into mPanelData
            return extraOutput;
        }
    }

    private void generateCopyNumberBackboneProbes()
    {
        if(mConfig.amberSitesFile() == null)
        {
            LOGGER.warn("Amber sites not provided; skipping copy number backbone probes");
        }
        else
        {
            CopyNumberBackbone copyNumberBackbone =
                    new CopyNumberBackbone(mConfig.amberSitesFile(), mRefGenomeVersion, mProbeGenerator, mPanelData);
            copyNumberBackbone.generateProbes();
            // Result is stored into mPanelData
        }
    }

    private void generateCustomRegionProbes()
    {
        if(mConfig.customRegionsFile() == null)
        {
            LOGGER.info("Custom regions not provided; skipping custom region probes");
        }
        else
        {
            CustomRegions.generateProbes(mConfig.customRegionsFile(), mRefGenome.chromosomeLengths(), mProbeGenerator, mPanelData);
            // Result is stored into mPanelData
        }
    }

    private SampleVariants.ExtraOutput generateSampleVariantProbes()
    {
        if(mConfig.sampleVariants() == null)
        {
            LOGGER.info("Sample data not provided; skipping sample variants probes");
            return null;
        }
        else
        {
            SampleVariants.ExtraOutput extraOutput =
                    SampleVariants.generateProbes(mConfig.sampleVariants(), mRefGenome, mProbeGenerator, mPanelData);
            // Result is stored into mPanelData
            return extraOutput;
        }
    }

    private EnsemblDataCache loadEnsemblData()
    {
        EnsemblDataCache ensemblData = new EnsemblDataCache(mConfig.ensemblDir(), mRefGenomeVersion);
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
            PanelBuilderConfig config = PanelBuilderConfig.fromConfigBuilder(configBuilder);
            PanelBuilderApplication panelBuilder = new PanelBuilderApplication(config);
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
