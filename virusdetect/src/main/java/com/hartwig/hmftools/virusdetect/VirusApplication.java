package com.hartwig.hmftools.virusdetect;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.virusdetect.VirusConstants.ALIGNED_BAM_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.APP_NAME;
import static com.hartwig.hmftools.virusdetect.VirusConstants.CANDIDATE_FASTA_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.CONTIG_STATS_TSV_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.DECOY_CONTIGS;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_SOFT_CLIP_BASES_DEFAULT;

import java.io.IOException;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class VirusApplication
{
    private final VirusConfig mConfig;
    private final ViralReference mViralReference;
    private final CandidateReadExtractor mCandidateExtractor;
    private final ViralReadAligner mAligner;

    private static final Logger LOGGER = LogManager.getLogger(VirusApplication.class);

    public VirusApplication(VirusConfig config)
    {
        mConfig = config;

        LOGGER.info("loading viral reference model");
        mViralReference = ViralReference.load(config.viralRefFile(), config.viralRefInfoFile());
        LOGGER.info("viral reference model loaded");

        CandidateReadFilter candidateFilter = new CandidateReadFilter(MIN_SOFT_CLIP_BASES_DEFAULT, DECOY_CONTIGS);
        mCandidateExtractor = new CandidateReadExtractor(config.refGenomeFile(), candidateFilter);
        mAligner = ViralReadAligner.create(config, mViralReference);
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting VirusDetect");
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        LOGGER.info("Extracting candidate viral reads from tumor BAM");
        String candidateFastaFile = candidateFastaFile();
        mCandidateExtractor.extractToFasta(mConfig.tumorBam(), candidateFastaFile);
        LOGGER.info("Candidate read extraction complete");

        LOGGER.info("Aligning candidate reads to viral reference");
        String alignedBamFile = alignedBamFile();
        mAligner.align(candidateFastaFile, alignedBamFile);
        LOGGER.info("Alignment complete");

        LOGGER.info("Computing per-contig statistics");
        Map<String, ContigStats> contigStats = new ContigStatsCalculator().compute(alignedBamFile, mViralReference);
        VirusOutputWriter.writeContigStats(contigStatsFile(), contigStats.values(), mViralReference);
        LOGGER.info("Per-contig statistics complete");

        // TODO: placeholder pipeline; each step is replaced by its implementation as it lands.
        LOGGER.info("Selecting representative contig per oncology group (stub)");
        LOGGER.info("Filtering aligned BAM to representatives -> BAM (stub)");
        LOGGER.info("Computing per-contig stats over representative BAM (stub)");
        LOGGER.info("Annotating QC and writing detected TSV (stub)");

        LOGGER.info("VirusDetect complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private String candidateFastaFile()
    {
        return mConfig.outputDir() + mConfig.sampleId() + CANDIDATE_FASTA_SUFFIX;
    }

    private String alignedBamFile()
    {
        return mConfig.outputDir() + mConfig.sampleId() + ALIGNED_BAM_SUFFIX;
    }

    private String contigStatsFile()
    {
        return mConfig.outputDir() + mConfig.sampleId() + CONTIG_STATS_TSV_SUFFIX;
    }

    public static void main(@NotNull String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        VirusConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        try
        {
            VirusConfig config = VirusConfig.fromConfigBuilder(configBuilder);
            VirusApplication virusDetect = new VirusApplication(config);
            virusDetect.run();
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
