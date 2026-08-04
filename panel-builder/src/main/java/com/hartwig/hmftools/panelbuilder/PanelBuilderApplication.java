package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.round;
import static java.lang.System.exit;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.APP_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.Utils.estimatePanelOnTargetRate;
import static com.hartwig.hmftools.panelbuilder.probequality.Utils.createBwaMemAligner;

import java.io.IOException;
import java.util.List;
import java.util.function.Supplier;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.CachedRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PanelBuilderApplication
{
    private final PanelBuilderConfig mConfig;
    private final RefGenomeInterface mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
    private final ProbeGenerator mProbeGenerator;
    private final ProbeGenerator mRnaProbeGenerator;
    private PanelData mPanelData;
    // RNA probes are kept in a separate panel so RNA and DNA coverage/overlap never interact.
    private PanelData mRnaPanelData;
    // Loaded on first use and shared between DNA and RNA gene probe generation.
    @Nullable
    private EnsemblDataCache mEnsemblData;
    @Nullable
    private OutputWriter mOutputWriter;

    private static final Logger LOGGER = LogManager.getLogger(PanelBuilderApplication.class);

    public PanelBuilderApplication(final PanelBuilderConfig config)
    {
        mConfig = config;

        LOGGER.info("Loading prerequisite data");

        RefGenomeSource refGenomeSource = loadRefGenome(mConfig.refGenomeFile());
        mRefGenome = new CachedRefGenome(refGenomeSource, 2000, 10);
        mRefGenomeVersion = deriveRefGenomeVersion(refGenomeSource);

        ProbeQualityProfile probeQualityProfile = ProbeQualityProfile.loadFromResourceFile(mConfig.probeQualityProfileFile());
        loadAlignerLibrary(mConfig.bwaLibPath());
        Supplier<BwaMemAligner> alignerFactory = () -> createBwaMemAligner(mConfig.bwaIndexImageFile(), mConfig.threads());
        ProbeQualityModel probeQualityModel = new ProbeQualityModel(
                alignerFactory, PROBE_LENGTH,
                probeQualityProfile.matchScoreThreshold(), probeQualityProfile.matchScoreOffset());

        mProbeGenerator = ProbeGenerator.construct(mRefGenome, probeQualityProfile, probeQualityModel, this::writeCandidateProbe);
        mRnaProbeGenerator = ProbeGenerator.construct(mRefGenome, probeQualityProfile, probeQualityModel, this::writeRnaCandidateProbe);
        mPanelData = new PanelData();
        mRnaPanelData = new PanelData();
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting PanelBuilder");

        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());
        mOutputWriter = new OutputWriter(
                mConfig.outputDir(), mConfig.outputId(), mConfig.verboseOutput(), mConfig.rnaGenesFile() != null);

        LOGGER.info("Generating probes");
        mPanelData = new PanelData();
        mRnaPanelData = new PanelData();
        // Note the order of generation here determines the priority of probe overlap resolution.
        // Probes generated first will exclude overlapping probes generated afterward.
        // RNA probes are a separate panel (separate PanelData) and don't interact with the DNA panel overlap resolution.
        Genes.ExtraOutput genesExtraOutput = generateTargetGeneProbes();
        GenesRna.ExtraOutput rnaGenesExtraOutput = generateRnaGeneProbes();
        generateCustomRegionProbes();
        generateCustomStructuralVariantProbes();
        generateCustomSmallVariantProbes();
        generateCopyNumberBackboneProbes();
        generateCdr3Probes();
        SampleVariants.ExtraOutput sampleVariantsExtraOutput = generateSampleVariantProbes();

        LOGGER.info("Writing output");
        ProbeOutputWriter dnaOutput = mOutputWriter.panelOutput();
        dnaOutput.writeProbes(mPanelData.probes());
        dnaOutput.writeCoveredTargetRegions(mPanelData.coveredTargetRegions());
        dnaOutput.writeCoveredRegions(mPanelData.probes());
        dnaOutput.writeCandidateTargetRegions(mPanelData.candidateTargetRegions());
        dnaOutput.writeRejectedFeatures(mPanelData.rejectedFeatures());
        if(genesExtraOutput != null)
        {
            dnaOutput.writeGeneStats(genesExtraOutput.geneStats());
        }
        if(rnaGenesExtraOutput != null)
        {
            ProbeOutputWriter rnaOutput = requireNonNull(mOutputWriter.rnaPanelOutput());
            rnaOutput.writeProbes(mRnaPanelData.probes());
            rnaOutput.writeCoveredTargetRegions(mRnaPanelData.coveredTargetRegions());
            rnaOutput.writeCoveredRegions(mRnaPanelData.probes());
            rnaOutput.writeCandidateTargetRegions(mRnaPanelData.candidateTargetRegions());
            rnaOutput.writeRejectedFeatures(mRnaPanelData.rejectedFeatures());
            rnaOutput.writeGeneStats(rnaGenesExtraOutput.geneStats());
        }
        if(sampleVariantsExtraOutput != null)
        {
            mOutputWriter.writeSampleVariantInfos(sampleVariantsExtraOutput.variantInfos());
        }
        mOutputWriter.close();
        mOutputWriter = null;

        printPanelStats();

        LOGGER.info("PanelBuilder complete, mins({})", runTimeMinsStr(startTimeMs));
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
            if(mConfig.ensemblDir() == null)
            {
                throw new UserInputError("Genes requested but Ensembl data directory not provided");
            }
            {
                EnsemblDataCache ensemblData = ensemblData();
                Genes.ExtraOutput extraOutput =
                        Genes.generateProbes(mConfig.genesFile(), ensemblData, mProbeGenerator, mPanelData);
                // Result is stored into mPanelData.
                return extraOutput;
            }
        }
    }

    @Nullable
    private GenesRna.ExtraOutput generateRnaGeneProbes()
    {
        if(mConfig.rnaGenesFile() == null)
        {
            LOGGER.info("Genes RNA not provided; skipping gene RNA probes");
            return null;
        }
        else
        {
            if(mConfig.ensemblDir() == null)
            {
                throw new UserInputError("Genes RNA requested but Ensembl data directory not provided");
            }
            EnsemblDataCache ensemblData = ensemblData();
            GenesRna.ExtraOutput extraOutput =
                    GenesRna.generateProbes(mConfig.rnaGenesFile(), ensemblData, mRnaProbeGenerator, mRnaPanelData);
            // Result is stored into mRnaPanelData.
            return extraOutput;
        }
    }

    private void generateCopyNumberBackboneProbes()
    {
        if(mConfig.includeCnBackbone())
        {
            if(mConfig.hetSitesFile() == null)
            {
                throw new UserInputError("Copy number backbone requested but heterozygous sites file not provided");
            }
            else
            {
                new CopyNumberBackbone(
                        mConfig.hetSitesFile(), mConfig.cnBackboneResolution(), mRefGenomeVersion, mProbeGenerator, mPanelData)
                        .generateProbes();
                // Result is stored into mPanelData.
            }
        }
        else
        {
            LOGGER.info("Copy number backbone not requested; skipping copy number backbone probes");
        }
    }

    private void generateCdr3Probes()
    {
        if(mConfig.includeCdr3())
        {
            Cdr3Regions.generateProbes(mRefGenomeVersion, mProbeGenerator, mPanelData);
            // Result is stored into mPanelData.
        }
        else
        {
            LOGGER.info("CDR3 not requested; skipping CDR3 probe generation");
        }
    }

    private void generateCustomRegionProbes()
    {
        if(mConfig.customRegionsFiles() == null)
        {
            LOGGER.info("Custom regions not provided; skipping custom region probes");
        }
        else
        {
            CustomRegions.generateProbes(mConfig.customRegionsFiles(), mRefGenome.chromosomeLengths(), mProbeGenerator, mPanelData);
            // Result is stored into mPanelData.
        }
    }

    private void generateCustomSmallVariantProbes()
    {
        if(mConfig.customSmallVariantsFile() == null)
        {
            LOGGER.info("Custom small variants not provided; skipping custom small variant probes");
        }
        else
        {
            CustomSmallVariants.generateProbes(mConfig.customSmallVariantsFile(), mRefGenome, mProbeGenerator, mPanelData);
            // Result is stored into mPanelData.
        }
    }

    private void generateCustomStructuralVariantProbes()
    {
        if(mConfig.customStructuralVariantsFile() == null)
        {
            LOGGER.info("Custom structural variants not provided; skipping custom structural variant probes");
        }
        else
        {
            CustomStructuralVariants.generateProbes(
                    mConfig.customStructuralVariantsFile(), mRefGenome.chromosomeLengths(), mProbeGenerator, mPanelData);
            // Result is stored into mPanelData.
        }
    }

    @Nullable
    private SampleVariants.ExtraOutput generateSampleVariantProbes()
    {
        if(mConfig.sampleVariants() == null)
        {
            LOGGER.info("Sample data not provided; skipping sample variants probes");
            return null;
        }
        else
        {
            SampleVariants.ExtraOutput extraOutput = SampleVariants.generateProbes(mConfig.sampleVariants(), mProbeGenerator, mPanelData);
            // Result is stored into mPanelData.
            return extraOutput;
        }
    }

    // Loads the Ensembl cache once and shares it between DNA and RNA gene probe generation.
    private EnsemblDataCache ensemblData()
    {
        if(mEnsemblData == null)
        {
            EnsemblDataCache ensemblData = new EnsemblDataCache(mConfig.ensemblDir(), mRefGenomeVersion);
            ensemblData.setRequiredData(true, false, false, false);
            ensemblData.load(false);
            mEnsemblData = ensemblData;
        }
        return mEnsemblData;
    }

    private void writeCandidateProbe(final Probe probe)
    {
        requireNonNull(mOutputWriter).writeCandidateProbe(probe);
    }

    private void writeRnaCandidateProbe(final Probe probe)
    {
        requireNonNull(mOutputWriter).writeRnaCandidateProbe(probe);
    }

    private void printPanelStats()
    {
        List<Probe> probes = mPanelData.probes();
        long probeBases = probes.stream().mapToLong(probe -> probe.definition().baseLength()).sum();
        double estimatedOnTargetRate = estimatePanelOnTargetRate(probes) * 100;
        LOGGER.info("Panel stats:");
        LOGGER.info("  Probe count: {}", probes.size());
        LOGGER.info("  Probe bases: {}", probeBases);
        LOGGER.info("  Estimated on-target rate: {}%", round(estimatedOnTargetRate));

        List<Probe> rnaProbes = mRnaPanelData.probes();
        if(!rnaProbes.isEmpty())
        {
            long rnaProbeBases = rnaProbes.stream().mapToLong(probe -> probe.definition().baseLength()).sum();
            LOGGER.info("  RNA probe count: {}", rnaProbes.size());
            LOGGER.info("  RNA probe bases: {}", rnaProbeBases);
        }
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
        catch(UserInputError e)
        {
            LOGGER.error("Bad input data: {}", e.getMessage());
            exit(1);
        }
        catch(IOException e)
        {
            LOGGER.error("IO error", e);
            exit(1);
        }
        catch(RuntimeException e)
        {
            LOGGER.error("Runtime error", e);
            exit(1);
        }
    }
}
