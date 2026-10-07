package com.hartwig.hmftools.viridian.app;

import static java.lang.System.exit;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.summingInt;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.ALL_ALIGNMENTS_BAM_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.APP_NAME;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.CANDIDATES_FASTA_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.CONTIG_SUPPORT_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATIONS_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.PAIRWISE_MARGINS_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_MIN_SOFT_CLIP_BASES;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_DETECTION_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_REF_CONTIGS;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_REPORT_TSV_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.viridian.common.UserInputError;
import com.hartwig.hmftools.viridian.detection.DetectedVirus;
import com.hartwig.hmftools.viridian.detection.align.AllContigsReadAlignments;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAligner;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignmentStore;
import com.hartwig.hmftools.viridian.detection.assign.RepresentativeReadAssigner;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.ContigStatsCalculator;
import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.detection.extract.CandidateReadExtractor;
import com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilter;
import com.hartwig.hmftools.viridian.detection.extract.ViralKmerIndex;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.detection.select.PairwiseMargins;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeContigSelector;
import com.hartwig.hmftools.viridian.detection.support.ContigSupport;
import com.hartwig.hmftools.viridian.detection.support.ContigSupportCalculator;
import com.hartwig.hmftools.viridian.integration.Integration;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAligner;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.HostVariantCandidate;
import com.hartwig.hmftools.viridian.integration.extract.HostVariantExtractor;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reference.VirusReference;
import com.hartwig.hmftools.viridian.reporting.ClonalCoverage;
import com.hartwig.hmftools.viridian.reporting.VirusReport;
import com.hartwig.hmftools.viridian.reporting.VirusReporter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViridianApplication
{
    private final ViridianConfig mConfig;
    private final VirusReference mVirusReference;

    private static final Logger LOGGER = LogManager.getLogger(ViridianApplication.class);

    public ViridianApplication(ViridianConfig config)
    {
        mConfig = config;

        LOGGER.info("Loading virus reference data");
        mVirusReference = VirusReference.load(config.virusRefFile(), config.virusRefInfoFile(), config.oncologyGroupInfoFile());
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting {}", APP_NAME);
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        BwaMemAligner.initLibrary(mConfig.bwaLibPath());

        AllContigsReadAlignments allContigsReadAlignments = getAllContigsReadAlignments();

        List<ContigSupport> contigSupports = computeViralContigSupport(allContigsReadAlignments);

        List<OncologyGroupRepresentativeSelection> representativeSelections =
                selectRepresentativeViralContigs(allContigsReadAlignments, contigSupports);

        ViralReadAlignmentStore representativeAlignments =
                assignReadsToRepresentatives(allContigsReadAlignments.store(), representativeSelections);

        Map<ViralContig, ContigStats> representativeContigStats =
                calculateRepresentativeContigStats(representativeAlignments, allContigsReadAlignments.metrics().originClippedReads());

        List<DetectedVirus> detectedViruses = writeDetectedViruses(
                representativeSelections, representativeContigStats, allContigsReadAlignments.metrics().readCountsByOncologyGroup());

        Map<OncologyGroup, Integer> integrationCounts = callHostIntegrations();

        reportViruses(detectedViruses, integrationCounts);

        LOGGER.info("{} complete, mins({})", APP_NAME, runTimeMinsStr(startTimeMs));
    }

    private AllContigsReadAlignments getAllContigsReadAlignments()
    {
        String allAlignmentsBamFile = outputFile(ALL_ALIGNMENTS_BAM_SUFFIX);
        // Alignment is pretty slow, so allow reusing the cached BAM for a rerun.
        if(!canReuseExistingFile(mConfig.reuseReadAlignments(), allAlignmentsBamFile, "all-alignments BAM"))
        {
            alignCandidateReadsToViralContigs(allAlignmentsBamFile);
        }
        return AllContigsReadAlignments.load(allAlignmentsBamFile, mVirusReference);
    }

    private String getCandidateReads()
    {
        String candidateFastaFile = outputFile(CANDIDATES_FASTA_SUFFIX);
        // Read extraction is very slow for large samples, so allow reusing the cached FASTA for a rerun.
        if(!canReuseExistingFile(mConfig.reuseCandidateReads(), candidateFastaFile, "candidate read FASTA"))
        {
            extractCandidateReads(candidateFastaFile);
        }
        return candidateFastaFile;
    }

    // Extract reads which may be viral into a FASTA.
    private void extractCandidateReads(String candidateFastaFile)
    {
        ViralKmerIndex kmerIndex = null;
        if(mConfig.kmerFilterEnabled())
        {
            LOGGER.info("Building virus genome k-mer index");
            kmerIndex = ViralKmerIndex.build(mConfig.virusRefFile());
        }

        LOGGER.info("Extracting candidate viral reads from tumor BAM using read source {}", mConfig.candidateReadSource());
        CandidateReadFilter candidateFilter = new CandidateReadFilter(
                mConfig.candidateReadSource(), VIRAL_READ_MIN_SOFT_CLIP_BASES, VIRUS_REF_CONTIGS, kmerIndex);
        CandidateReadExtractor mCandidateExtractor = new CandidateReadExtractor(
                mConfig.refGenomeFile(), candidateFilter, mConfig.threads());
        mCandidateExtractor.extract(mConfig.tumorBam(), candidateFastaFile);
    }

    // Align potentially viral reads to all virus genomes, so we can decide which viruses are present.
    private void alignCandidateReadsToViralContigs(String allAlignmentsBamFile)
    {
        String candidateReadFasta = getCandidateReads();

        LOGGER.info("Aligning candidate reads to virus genomes");
        ViralReadAligner viralReadAligner = ViralReadAligner.create(
                mVirusReference, mConfig.virusBwaIndexImage(), mConfig.threads(), mConfig.alignmentBatchSize());
        viralReadAligner.align(candidateReadFasta, allAlignmentsBamFile);
    }

    // Compute support information for each virus genome and decide which genomes may be present.
    private List<ContigSupport> computeViralContigSupport(AllContigsReadAlignments alignments)
    {
        LOGGER.info("Computing per-contig support");
        Map<ViralContig, ContigStats> contigStats = ContigStatsCalculator.calculate(
                alignments.store().byContig(), alignments.metrics().originClippedReads());
        return new ContigSupportCalculator().calculate(alignments, contigStats);
    }

    // For each oncology group (group of virus strains at interesting taxonomy granularity), select 1 virus genome which best represents
    // the virus present in the sample.
    private List<OncologyGroupRepresentativeSelection> selectRepresentativeViralContigs(
            AllContigsReadAlignments alignments, List<ContigSupport> viralContigSupports)
    {
        LOGGER.info("Selecting representative contig per oncology group");
        PairwiseMargins pairwiseMargins = PairwiseMargins.from(alignments.store().byRead());
        List<OncologyGroupRepresentativeSelection> selections = RepresentativeContigSelector.select(
                viralContigSupports, pairwiseMargins, alignments.metrics().readCountsByOncologyGroup());

        LOGGER.info("Writing contig info output");
        OutputWriter.writeContigSupport(outputFile(CONTIG_SUPPORT_TSV_SUFFIX), selections);
        if(mConfig.verboseOutput())
        {
            LOGGER.info("Writing pairwise margins output");
            OutputWriter.writePairwiseMargins(outputFile(PAIRWISE_MARGINS_TSV_SUFFIX), pairwiseMargins, selections);
        }

        return selections;
    }

    // For each oncology group, assign the read alignments to only the selected representative virus genome.
    // This is simply filtering down to each read's alignment to that genome contig (or nothing, if it didn't align there at all).
    private ViralReadAlignmentStore assignReadsToRepresentatives(
            ViralReadAlignmentStore alignments, List<OncologyGroupRepresentativeSelection> representativeSelections)
    {
        LOGGER.info("Assigning reads to representative contigs");
        Set<ViralContig> representatives = representativeSelections.stream()
                .map(OncologyGroupRepresentativeSelection::representative)
                .filter(Objects::nonNull)
                .collect(toSet());
        Map<ReadId, ViralReadAlignment> assignments = RepresentativeReadAssigner.assign(
                alignments, representatives, outputFile(ALL_ALIGNMENTS_BAM_SUFFIX),
                outputFile(REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX));

        return new ViralReadAlignmentStore(assignments.values());
    }

    // After the representative virus genome is selected for each oncology group, measure final support statistics for the representative.
    private Map<ViralContig, ContigStats> calculateRepresentativeContigStats(
            ViralReadAlignmentStore representativeAlignments, Map<ViralContig, Integer> allOriginClippedReads)
    {
        LOGGER.info("Computing representative contig stats");
        // Origin clipped reads where dropped previously, but the information is useful to carry through for visualisation.
        Map<ViralContig, Integer> originClippedReads = allOriginClippedReads.entrySet().stream()
                .filter(entry -> representativeAlignments.byContig().containsKey(entry.getKey()))
                .collect(toMap(Map.Entry::getKey, Map.Entry::getValue));
        Map<ViralContig, ContigStats> contigStats =
                ContigStatsCalculator.calculate(representativeAlignments.byContig(), originClippedReads);
        return contigStats;
    }

    private List<DetectedVirus> writeDetectedViruses(
            List<OncologyGroupRepresentativeSelection> selections, Map<ViralContig, ContigStats> representativeContigStats,
            Map<OncologyGroup, Integer> oncologyGroupReadCounts)
    {
        LOGGER.info("Writing virus detection output");
        List<DetectedVirus> detections = DetectedVirus.from(selections, representativeContigStats, oncologyGroupReadCounts);
        OutputWriter.writeDetectedViruses(outputFile(VIRUS_DETECTION_TSV_SUFFIX), detections);
        return detections;
    }

    // Find where a virus inserted itself into the host genome, from the SVs ESVEE called.
    // Returns the count of integrations per oncology group.
    private Map<OncologyGroup, Integer> callHostIntegrations()
    {
        LOGGER.info("Extracting integration variant candidates from ESVEE VCF");
        List<HostVariantCandidate> candidates = new HostVariantExtractor(mConfig.sampleId()).extract(mConfig.esveeUnfilteredVcf());

        List<Integration> integrations = alignHostVariantCandidates(candidates);

        LOGGER.info("Writing integrations output");
        // Note that all integrations are written for informational purposes, but only the plausibly aligned
        // integrations are used.
        OutputWriter.writeIntegrations(outputFile(INTEGRATIONS_TSV_SUFFIX), integrations);

        // Take the integrations with plausible alignments. Only the counts are needed for reporting.
        return integrations.stream()
                .filter(Integration::isPlausible)
                .collect(groupingBy(i -> i.alignment().contig().oncologyGroup(), HashMap::new, summingInt(i -> 1)));
    }

    private List<Integration> alignHostVariantCandidates(List<HostVariantCandidate> candidates)
    {
        LOGGER.info("Aligning candidate viral integration sequences to virus genomes");

        List<String> insertSequences = candidates.stream().map(HostVariantCandidate::insertSequence).toList();
        ViralInsertAligner aligner = ViralInsertAligner.create(
                mVirusReference, mConfig.virusBwaIndexImage(), mConfig.threads());
        List<ViralInsertAlignment> alignments = aligner.align(insertSequences);

        List<Integration> integrations = new ArrayList<>(candidates.size());
        for(int i = 0; i < candidates.size(); ++i)
        {
            integrations.add(new Integration(candidates.get(i), alignments.get(i)));
        }
        return integrations;
    }

    // Decide which viruses are reported to the downstream pipeline.
    private List<VirusReport> reportViruses(List<DetectedVirus> detectedViruses, Map<OncologyGroup, Integer> integrationCounts)
            throws IOException
    {
        LOGGER.info("Loading reporting inputs");
        PurityContext purity = PurityContextFile.read(mConfig.purpleDir(), mConfig.sampleId());
        BamMetricSummary tumorMetrics = BamMetricSummary.read(
                BamMetricSummary.generateFilename(mConfig.bamMetricsTumorDir(), mConfig.sampleId()));

        Double expectedViralDepthPerCopy = ClonalCoverage.expectedViralDepthPerCopy(purity, tumorMetrics);
        if(expectedViralDepthPerCopy == null)
        {
            LOGGER.warn("Purple fit unusable; cannot report based on viral depth");
        }

        LOGGER.info("Calculating reporting statuses");
        List<VirusReport> reports = VirusReporter.report(
                detectedViruses, integrationCounts, mVirusReference::oncologyGroupInfo, expectedViralDepthPerCopy);

        LOGGER.info("Writing virus report output");
        OutputWriter.writeVirusReports(outputFile(VIRUS_REPORT_TSV_SUFFIX), reports);

        return reports;
    }

    private String outputFile(String suffix)
    {
        String outputId = mConfig.outputId();
        String f = mConfig.outputDir() + mConfig.sampleId();
        if(outputId != null)
        {
            f += "." + outputId;
        }
        f += "." + APP_NAME.toLowerCase() + suffix;
        return f;
    }

    // Dev/debug skip: if requested, if the file is already cached, use that rather than recomputing it.
    private static boolean canReuseExistingFile(boolean reuseRequested, String file, String description)
    {
        if(!reuseRequested)
        {
            return false;
        }
        // Check the length rather than only existence, as zero-length files may be left from a run that was stopped midway.
        else if(new File(file).length() > 0)
        {
            LOGGER.info("Reusing existing {}: {}", description, file);
            return true;
        }
        else
        {
            LOGGER.debug("Existing {} not present, regenerating: {}", description, file);
            return false;
        }
    }

    public static void main(@NotNull String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        ViridianConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        try
        {
            ViridianConfig config = ViridianConfig.fromConfigBuilder(configBuilder);
            ViridianApplication app = new ViridianApplication(config);
            app.run();
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
