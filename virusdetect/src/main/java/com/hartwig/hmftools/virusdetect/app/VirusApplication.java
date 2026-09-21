package com.hartwig.hmftools.virusdetect.app;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.ALIGNED_READ_BAM_SUFFIX;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.APP_NAME;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.CANDIDATE_READ_FASTA_SUFFIX;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.CONTIG_INFO_TSV_SUFFIX;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.PAIRWISE_MARGINS_TSV_SUFFIX;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT;
import static com.hartwig.hmftools.virusdetect.common.VirusConstants.VIRAL_REF_CONTIGS;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.virusdetect.common.UserInputError;
import com.hartwig.hmftools.virusdetect.detection.contig_support.ContigSupport;
import com.hartwig.hmftools.virusdetect.detection.contig_support.ContigSupportCalculator;
import com.hartwig.hmftools.virusdetect.detection.read_align.ViralReadAligner;
import com.hartwig.hmftools.virusdetect.detection.read_align.ViralReadAlignments;
import com.hartwig.hmftools.virusdetect.detection.read_extract.CandidateReadExtractor;
import com.hartwig.hmftools.virusdetect.detection.read_extract.CandidateReadFilter;
import com.hartwig.hmftools.virusdetect.integration.seq_align.ViralSequenceAligner;
import com.hartwig.hmftools.virusdetect.integration.seq_align.ViralSequenceAlignment;
import com.hartwig.hmftools.virusdetect.integration.variant_extract.CandidateIntegration;
import com.hartwig.hmftools.virusdetect.integration.variant_extract.CandidateIntegrationExtractor;
import com.hartwig.hmftools.virusdetect.reference.ViralContig;
import com.hartwig.hmftools.virusdetect.reference.ViralReference;
import com.hartwig.hmftools.virusdetect.selection.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.virusdetect.selection.OncologyGroupResolution;
import com.hartwig.hmftools.virusdetect.selection.PairwiseMargins;
import com.hartwig.hmftools.virusdetect.selection.RepresentativeContigSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// TODO: rename when app name is decided
public class VirusApplication
{
    private final VirusConfig mConfig;
    private final ViralReference mViralReference;

    private static final Logger LOGGER = LogManager.getLogger(VirusApplication.class);

    public VirusApplication(VirusConfig config)
    {
        mConfig = config;

        LOGGER.info("Loading viral reference data");
        mViralReference = ViralReference.load(config.viralRefFile(), config.viralRefInfoFile());
        LOGGER.info("Viral reference data loaded");
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting VirusDetect");
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        BwaMemAligner.initLibrary(mConfig.bwaLibPath());

        String candidateReadFasta = extractCandidateReads();

        String viralReadBam = alignCandidateReadsToViralContigs(candidateReadFasta);

        ViralReadAlignments viralReadAlignments = ViralReadAlignments.load(viralReadBam, mViralReference);

        List<ContigSupport> viralContigSupports = computeViralContigSupport(viralReadAlignments);

        selectRepresentativeViralContigs(viralReadAlignments, viralContigSupports);

        callHostIntegrations();

        LOGGER.info("VirusDetect complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // Extract reads which may be viral into a FASTA.
    private String extractCandidateReads()
    {
        String candidateFastaFile = outputFile(CANDIDATE_READ_FASTA_SUFFIX);
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
            CandidateReadFilter candidateFilter = new CandidateReadFilter(VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT, VIRAL_REF_CONTIGS);
            CandidateReadExtractor mCandidateExtractor = new CandidateReadExtractor(
                    mConfig.refGenomeFile(), candidateFilter, mConfig.threads());
            mCandidateExtractor.extractToFasta(mConfig.tumorBam(), candidateFastaFile);
            LOGGER.info("Candidate read extraction complete");
        }
        return candidateFastaFile;
    }

    // Align potentially viral reads to all virus genomes, so we can decide which viruses are present.
    private String alignCandidateReadsToViralContigs(String candidateReadFasta)
    {
        LOGGER.info("Aligning candidate reads to viral genomes");
        ViralReadAligner viralReadAligner = ViralReadAligner.create(mConfig, mViralReference);
        String viralReadBamFile = outputFile(ALIGNED_READ_BAM_SUFFIX);
        viralReadAligner.align(candidateReadFasta, viralReadBamFile);
        LOGGER.info("Candidate read alignment complete");
        return viralReadBamFile;
    }

    // Compute support information for each virus genome and decide which genomes may be present.
    private List<ContigSupport> computeViralContigSupport(ViralReadAlignments viralReadAlignments)
    {
        LOGGER.info("Computing per-contig support");
        List<ContigSupport> contigSupports = new ContigSupportCalculator().compute(viralReadAlignments);
        LOGGER.info("Per-contig support complete");
        return contigSupports;
    }

    // For each oncology group (group of virus strains at interesting taxonomy granularity), select 1 viral genome which best represents
    // the virus present in the sample.
    private void selectRepresentativeViralContigs(ViralReadAlignments viralReadAlignments, List<ContigSupport> viralContigSupports)
    {
        LOGGER.info("Selecting representative contig per oncology group");
        PairwiseMargins pairwiseMargins = PairwiseMargins.from(viralReadAlignments);
        List<OncologyGroupRepresentativeSelection> selections = RepresentativeContigSelector.select(
                viralContigSupports, pairwiseMargins, viralReadAlignments.readCountsByOncologyGroup());
        logSelections(selections);

        OutputWriter.writeContigInfo(outputFile(CONTIG_INFO_TSV_SUFFIX), selections);
        if(mConfig.verboseOutput())
        {
            OutputWriter.writePairwiseMargins(outputFile(PAIRWISE_MARGINS_TSV_SUFFIX), pairwiseMargins, selections);
        }
        LOGGER.info("Selecting representative contigs complete");
    }

    // TODO: placeholder stages, each replaced by its implementation as it lands.
    private void callHostIntegrations()
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
        Map<CandidateIntegration, ViralSequenceAlignment> alignments = alignCandidateIntegrations(candidates);
        LOGGER.info("Aligned {} of {} inserted sequences to a viral contig", alignments.size(), candidates.size());

        LOGGER.info("Annotating {} integrations and writing TSV (stub)", candidates.size());
    }

    private Map<CandidateIntegration, ViralSequenceAlignment> alignCandidateIntegrations(List<CandidateIntegration> candidates)
    {
        List<String> insertSequences = candidates.stream().map(CandidateIntegration::insertSequence).toList();
        ViralSequenceAligner aligner = ViralSequenceAligner.create(mConfig, mViralReference);
        List<ViralSequenceAlignment> alignments = aligner.alignAll(insertSequences);

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
