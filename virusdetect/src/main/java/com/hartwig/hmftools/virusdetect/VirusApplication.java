package com.hartwig.hmftools.virusdetect;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.virusdetect.VirusConstants.ALIGNED_BAM_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.APP_NAME;
import static com.hartwig.hmftools.virusdetect.VirusConstants.CANDIDATE_FASTA_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.CONTIG_INFO_TSV_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.PAIRWISE_MARGINS_TSV_SUFFIX;
import static com.hartwig.hmftools.virusdetect.VirusConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;
import static com.hartwig.hmftools.virusdetect.VirusConstants.VIRAL_HOST_CONTIGS;
import static com.hartwig.hmftools.virusdetect.VirusConstants.VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.virusdetect.integration.CandidateIntegration;
import com.hartwig.hmftools.virusdetect.integration.CandidateIntegrationExtractor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// TODO: rename when app name is decided
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

        LOGGER.info("Loading viral reference model");
        mViralReference = ViralReference.load(config.viralRefFile(), config.viralRefInfoFile());
        LOGGER.info("Viral reference model loaded");

        CandidateReadFilter candidateFilter = new CandidateReadFilter(VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT, VIRAL_HOST_CONTIGS);
        mCandidateExtractor = new CandidateReadExtractor(config.refGenomeFile(), candidateFilter, config.threads());
        mAligner = ViralReadAligner.create(config, mViralReference);
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting VirusDetect");
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        String candidateFastaFile = outputFile(CANDIDATE_FASTA_SUFFIX);
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
        String alignedBamFile = outputFile(ALIGNED_BAM_SUFFIX);
        mAligner.align(candidateFastaFile, alignedBamFile);
        LOGGER.info("Alignment complete");

        ViralAlignments viralAlignments = ViralAlignments.load(alignedBamFile, mViralReference);

        LOGGER.info("Computing per-contig support");
        List<ContigSupport> contigSupports = new ContigSupportCalculator(READ_VOTE_CORRECT_BASE_PROBABILITY).compute(viralAlignments);
        LOGGER.info("Per-contig support complete");

        LOGGER.info("Selecting representative contig per oncology group");
        PairwiseMargins pairwiseMargins = PairwiseMargins.from(viralAlignments);
        List<OncologyGroupRepresentativeSelection> selections = RepresentativeSelector.select(
                contigSupports, pairwiseMargins, viralAlignments.readCountsByOncologyGroup());
        logSelections(selections);

        OutputWriter.writeContigInfo(outputFile(CONTIG_INFO_TSV_SUFFIX), selections);
        if(mConfig.verboseOutput())
        {
            OutputWriter.writePairwiseMargins(outputFile(PAIRWISE_MARGINS_TSV_SUFFIX), pairwiseMargins, selections);
        }

        callIntegrations();

        LOGGER.info("VirusDetect complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // TODO: placeholder stages, each replaced by its implementation as it lands.
    private void callIntegrations()
    {
        String esveeVcf = mConfig.esveeUnfilteredVcf();
        if(esveeVcf == null)
        {
            LOGGER.info("ESVEE VCF not specified; skipping integration site calling");
            return;
        }

        LOGGER.info("Extracting integration candidates from ESVEE VCF: {}", esveeVcf);
        List<CandidateIntegration> candidates = new CandidateIntegrationExtractor(mConfig.sampleId()).extract(esveeVcf);

        LOGGER.info("Aligning {} candidate viral integration sequences to the viral reference", candidates.size());
        Map<CandidateIntegration, ViralSequenceAlignment> alignments = alignInserts(candidates);
        LOGGER.info("Aligned {} of {} inserted sequences to a viral contig", alignments.size(), candidates.size());

        LOGGER.info("Annotating {} integrations and writing TSV (stub)", candidates.size());
    }

    // A candidate whose insert aligned nowhere is simply absent from the result, and is still reported.
    private Map<CandidateIntegration, ViralSequenceAlignment> alignInserts(List<CandidateIntegration> candidates)
    {
        List<String> insertSequences = candidates.stream().map(CandidateIntegration::insertSequence).toList();
        List<ViralSequenceAlignment> alignments = ViralSequenceAligner.create(mConfig, mViralReference).alignAll(insertSequences);

        Map<CandidateIntegration, ViralSequenceAlignment> byCandidate = new LinkedHashMap<>();
        for(int i = 0; i < candidates.size(); ++i)
        {
            if(alignments.get(i) != null)
            {
                byCandidate.put(candidates.get(i), alignments.get(i));
            }
        }
        return byCandidate;
    }

    private String outputFile(String suffix)
    {
        String outputId = mConfig.outputId();
        String f = mConfig.outputDir() + mConfig.sampleId();
        if(outputId != null)
        {
            f += "." + outputId;
        }
        f += suffix;
        return f;
    }

    private static void logSelections(List<OncologyGroupRepresentativeSelection> selections)
    {
        for(OncologyGroupRepresentativeSelection selection : selections)
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
