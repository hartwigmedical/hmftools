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
import static com.hartwig.hmftools.virusdetect.VirusConstants.PAIRWISE_MARGINS_TSV_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.util.List;
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
        mCandidateExtractor = new CandidateReadExtractor(config.refGenomeFile(), candidateFilter, config.threads());
        mAligner = ViralReadAligner.create(config, mViralReference);
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting VirusDetect");
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        String candidateFastaFile = candidateFastaFile();
        boolean reuseExisting = mConfig.reuseCandidateFasta() && new File(candidateFastaFile).exists();
        if(reuseExisting)
        {
            LOGGER.info("Reusing existing candidate FASTA: {}", candidateFastaFile);
        }
        else
        {
            if(mConfig.reuseCandidateFasta())
            {
                LOGGER.info("Candidate FASTA not found, extracting: {}", candidateFastaFile);
            }
            LOGGER.info("Extracting candidate viral reads from tumor BAM");
            mCandidateExtractor.extractToFasta(mConfig.tumorBam(), candidateFastaFile);
            LOGGER.info("Candidate read extraction complete");
        }

        LOGGER.info("Aligning candidate reads to viral reference");
        String alignedBamFile = alignedBamFile();
        mAligner.align(candidateFastaFile, alignedBamFile);
        LOGGER.info("Alignment complete");

        ViralAlignments viralAlignments = ViralAlignments.load(alignedBamFile, mViralReference);

        LOGGER.info("Computing per-contig statistics");
        Map<ViralContig, ContigSupport> contigStats = new ContigSupportCalculator().compute(viralAlignments);
        LOGGER.info("Per-contig statistics complete");

        LOGGER.info("Selecting representative contig per oncology group");
        PairwiseMargins pairwiseMargins = PairwiseMargins.from(viralAlignments);
        List<OncologyGroupSelection> selections = new RepresentativeSelector().select(
                contigStats.values(), pairwiseMargins, viralAlignments.readCountsByOncologyGroup());
        logSelections(selections);

        VirusOutputWriter.writeContigStats(contigStatsFile(), selections);
        if(mConfig.verboseOutput())
        {
            VirusOutputWriter.writePairwiseMargins(pairwiseMarginsFile(), pairwiseMargins, selections);
        }

        // TODO: placeholder pipeline; each step is replaced by its implementation as it lands.
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

    private String pairwiseMarginsFile()
    {
        return mConfig.outputDir() + mConfig.sampleId() + PAIRWISE_MARGINS_TSV_SUFFIX;
    }

    private static void logSelections(List<OncologyGroupSelection> selections)
    {
        for(OncologyGroupSelection selection : selections)
        {
            ViralContig representative = selection.representative();
            if(representative != null)
            {
                LOGGER.info("oncologyGroup({}) representative({})", selection.oncologyGroup(), representative.name());
            }
            else if(selection.resolution() == OncologyGroupResolution.UNRESOLVED)
            {
                LOGGER.warn("oncologyGroup({}) unresolved({})", selection.oncologyGroup(), selection.outcome());
            }
        }
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
