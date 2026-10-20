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
import static com.hartwig.hmftools.viridian.common.ViridianConstants.DEFAULT_VIRUS_DECOY_CONTIGS;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.GENOME_SUPPORT_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATION_CANDIDATES_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.READ_VOTE_CORRECT_BASE_PROBABILITY;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.RUN_MANIFEST_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_MIN_SOFT_CLIP_BASES;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_PROFILE_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_REPORT_TSV_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.viridian.common.UserInputError;
import com.hartwig.hmftools.viridian.integration.IntegrationCandidate;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAligner;
import com.hartwig.hmftools.viridian.integration.align.ViralInsertAlignment;
import com.hartwig.hmftools.viridian.integration.extract.CandidateSv;
import com.hartwig.hmftools.viridian.integration.extract.CandidateSvExtractor;
import com.hartwig.hmftools.viridian.profiling.VirusProfile;
import com.hartwig.hmftools.viridian.profiling.align.AllGenomesReadAlignments;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAligner;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.profiling.align.ViralReadAlignmentStore;
import com.hartwig.hmftools.viridian.profiling.assign.RepresentativeReadAssigner;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.GenomeStatsCalculator;
import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.profiling.extract.CandidateReadExtractor;
import com.hartwig.hmftools.viridian.profiling.extract.CandidateReadFilter;
import com.hartwig.hmftools.viridian.profiling.extract.ViralKmerIndex;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeGenome;
import com.hartwig.hmftools.viridian.profiling.select.RepresentativeGenomeSelector;
import com.hartwig.hmftools.viridian.profiling.support.GenomeSupport;
import com.hartwig.hmftools.viridian.profiling.support.GenomeSupportCalculator;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusReference;
import com.hartwig.hmftools.viridian.reporting.Clonality;
import com.hartwig.hmftools.viridian.reporting.VirusReport;
import com.hartwig.hmftools.viridian.reporting.VirusReporter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

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

        List<VirusProfile> profiles = profileViruses();

        Map<OncologyGroup, Integer> integrationCounts = callHostIntegrations();

        reportViruses(profiles, integrationCounts);

        OutputWriter.writeRunManifest(outputFile(RUN_MANIFEST_TSV_SUFFIX), mConfig);

        LOGGER.info("{} complete, mins({})", APP_NAME, runTimeMinsStr(startTimeMs));
    }

    // Find which viruses are present from the reads.
    private List<VirusProfile> profileViruses()
    {
        Set<String> virusDecoyContigs = resolveVirusDecoyContigs();

        String candidateReadFastaFile = extractCandidateReads(virusDecoyContigs);

        AllGenomesReadAlignments allGenomesReadAlignments = alignCandidateReadsToVirusGenomes(candidateReadFastaFile);

        List<GenomeSupport> genomeSupports = computeGenomeSupport(allGenomesReadAlignments);

        List<RepresentativeGenome> representatives = selectRepresentativeGenomes(genomeSupports);

        ViralReadAlignmentStore representativeAlignments =
                assignReadsToRepresentatives(allGenomesReadAlignments.store(), representatives);

        Map<VirusGenome, GenomeStats> representativeGenomeStats =
                calculateRepresentativeGenomeStats(representativeAlignments, allGenomesReadAlignments.metrics().originClippedReads());

        List<VirusProfile> profiles = VirusProfile.from(representatives, representativeGenomeStats);

        LOGGER.info("Writing virus profile output");
        OutputWriter.writeVirusProfiles(outputFile(VIRUS_PROFILE_TSV_SUFFIX), profiles);

        return profiles;
    }

    // Extract reads which may be viral into a FASTA.
    private String extractCandidateReads(Set<String> virusDecoyContigs)
    {
        String candidateFastaFile = outputFile(CANDIDATES_FASTA_SUFFIX);

        ViralKmerIndex kmerIndex = null;
        if(mConfig.kmerFilterEnabled())
        {
            LOGGER.info("Building virus genome k-mer index");
            kmerIndex = ViralKmerIndex.build(mConfig.virusRefFile());
        }

        LOGGER.info("Extracting candidate viral reads from tumor BAM using read source {}", mConfig.candidateReadSource());
        CandidateReadFilter candidateFilter = new CandidateReadFilter(
                mConfig.candidateReadSource(), VIRAL_READ_MIN_SOFT_CLIP_BASES, virusDecoyContigs, kmerIndex);
        CandidateReadExtractor candidateExtractor = new CandidateReadExtractor(
                mConfig.refGenomeFile(), candidateFilter, mConfig.threads());
        candidateExtractor.extract(mConfig.tumorBam(), candidateFastaFile);

        return candidateFastaFile;
    }

    private Set<String> resolveVirusDecoyContigs()
    {
        checkVirusDecoyContigsInTumorBam(mConfig.tumorBam(), mConfig.refGenomeFile(), mConfig.extraVirusDecoyContigs());
        Set<String> virusDecoyContigs = new HashSet<>(DEFAULT_VIRUS_DECOY_CONTIGS);
        virusDecoyContigs.addAll(mConfig.extraVirusDecoyContigs());
        return virusDecoyContigs;
    }

    static void checkVirusDecoyContigsInTumorBam(String tumorBam, @Nullable String refGenomeFile, Set<String> contigs)
    {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        if(refGenomeFile != null)
        {
            readerFactory = readerFactory.referenceSequence(new File(refGenomeFile));
        }
        try(SamReader reader = readerFactory.open(new File(tumorBam)))
        {
            SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();
            List<String> missingContigs = contigs.stream().filter(contig -> dictionary.getSequence(contig) == null).sorted().toList();
            if(!missingContigs.isEmpty())
            {
                throw new UserInputError("Virus decoy contigs not in BAM: " + String.join(", ", missingContigs));
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to open tumor BAM", e);
        }
    }

    // Align potentially viral reads to all virus genomes, so we can decide which viruses are present.
    private AllGenomesReadAlignments alignCandidateReadsToVirusGenomes(String candidateReadFasta)
    {
        LOGGER.info("Aligning candidate reads to virus genomes");
        String allAlignmentsBamFile = outputFile(ALL_ALIGNMENTS_BAM_SUFFIX);
        ViralReadAligner viralReadAligner = ViralReadAligner.create(
                mVirusReference, mConfig.virusBwaIndexImage(), mConfig.threads(), mConfig.alignmentBatchSize());
        viralReadAligner.align(candidateReadFasta, allAlignmentsBamFile);
        return AllGenomesReadAlignments.load(allAlignmentsBamFile, mVirusReference);
    }

    // Compute support information for each virus genome and decide which genomes are candidates.
    private List<GenomeSupport> computeGenomeSupport(AllGenomesReadAlignments alignments)
    {
        LOGGER.info("Computing per-genome support");
        Map<VirusGenome, GenomeStats> genomeStats = GenomeStatsCalculator.calculate(
                alignments.store().byGenome(), alignments.metrics().originClippedReads());
        return new GenomeSupportCalculator(READ_VOTE_CORRECT_BASE_PROBABILITY).calculate(alignments, genomeStats);
    }

    // For each oncology group (group of virus strains at interesting taxonomy granularity), select 1 virus genome which best represents
    // the virus present in the sample.
    private List<RepresentativeGenome> selectRepresentativeGenomes(List<GenomeSupport> genomeSupports)
    {
        LOGGER.info("Selecting representative genome per oncology group");
        List<RepresentativeGenome> representatives = RepresentativeGenomeSelector.select(genomeSupports);

        LOGGER.info("Writing genome support output");
        OutputWriter.writeGenomeSupport(outputFile(GENOME_SUPPORT_TSV_SUFFIX), genomeSupports, representatives);

        return representatives;
    }

    // Assign each read to the representative virus genome it fits best.
    private ViralReadAlignmentStore assignReadsToRepresentatives(
            ViralReadAlignmentStore alignments, List<RepresentativeGenome> representatives)
    {
        LOGGER.info("Assigning reads to representative genomes");
        Set<VirusGenome> representativeGenomes =
                representatives.stream().map(RepresentativeGenome::genome).collect(toSet());
        Map<ReadId, ViralReadAlignment> assignments = RepresentativeReadAssigner.assign(
                alignments, representativeGenomes, outputFile(ALL_ALIGNMENTS_BAM_SUFFIX),
                outputFile(REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX));

        return new ViralReadAlignmentStore(assignments.values());
    }

    // After the representative virus genome is selected for each oncology group, measure final support statistics for the representative.
    private Map<VirusGenome, GenomeStats> calculateRepresentativeGenomeStats(
            ViralReadAlignmentStore representativeAlignments, Map<VirusGenome, Integer> allOriginClippedReads)
    {
        LOGGER.info("Computing representative genome stats");
        // Origin clipped reads were dropped previously, but the information is useful to carry through for visualisation.
        Map<VirusGenome, Integer> originClippedReads = allOriginClippedReads.entrySet().stream()
                .filter(entry -> representativeAlignments.byGenome().containsKey(entry.getKey()))
                .collect(toMap(Map.Entry::getKey, Map.Entry::getValue));
        Map<VirusGenome, GenomeStats> genomeStats =
                GenomeStatsCalculator.calculate(representativeAlignments.byGenome(), originClippedReads);
        return genomeStats;
    }

    // Find where a virus inserted itself into the host genome, from the SVs ESVEE called.
    // Returns the count of integrations per oncology group.
    private Map<OncologyGroup, Integer> callHostIntegrations()
    {
        LOGGER.info("Extracting integration candidate SVs from ESVEE VCF");
        List<CandidateSv> candidateSvs =
                new CandidateSvExtractor(mConfig.sampleId()).extract(mConfig.esveeUnfilteredVcf());

        List<IntegrationCandidate> integrationCandidates = alignCandidateSvs(candidateSvs);

        LOGGER.info("Writing integration candidates output");
        // Note that all integration candidates are written for informational purposes.
        // But only the aligned candidates are the called integrations which flow through to the rest of the code.
        OutputWriter.writeIntegrationCandidates(outputFile(INTEGRATION_CANDIDATES_TSV_SUFFIX), integrationCandidates);

        // Take the integration candidates with plausible alignments. Only the counts are needed for reporting.
        return integrationCandidates.stream()
                .filter(IntegrationCandidate::isIntegration)
                .collect(groupingBy(i -> i.alignment().genome().oncologyGroup(), HashMap::new, summingInt(i -> 1)));
    }

    private List<IntegrationCandidate> alignCandidateSvs(List<CandidateSv> candidateSvs)
    {
        LOGGER.info("Aligning integration candidate insert sequences to virus genomes");

        List<String> insertSequences = candidateSvs.stream().map(CandidateSv::insertSequence).toList();
        ViralInsertAligner aligner = ViralInsertAligner.create(
                mVirusReference, mConfig.virusBwaIndexImage(), mConfig.threads());
        List<ViralInsertAlignment> alignments = aligner.align(insertSequences);

        List<IntegrationCandidate> integrationCandidates = new ArrayList<>(candidateSvs.size());
        for(int i = 0; i < candidateSvs.size(); ++i)
        {
            integrationCandidates.add(new IntegrationCandidate(candidateSvs.get(i), alignments.get(i)));
        }
        return integrationCandidates;
    }

    // Decide which viruses are reported to the downstream pipeline.
    private List<VirusReport> reportViruses(List<VirusProfile> profiles, Map<OncologyGroup, Integer> integrationCounts)
            throws IOException
    {
        LOGGER.info("Loading reporting inputs");
        PurityContext purity = PurityContextFile.read(mConfig.purpleDir(), mConfig.sampleId());
        BamMetricSummary tumorMetrics = BamMetricSummary.read(
                BamMetricSummary.generateFilename(mConfig.bamMetricsTumorDir(), mConfig.sampleId()));

        Double expectedViralDepthPerCopy = Clonality.expectedViralDepthPerCopy(purity, tumorMetrics);
        if(expectedViralDepthPerCopy == null)
        {
            LOGGER.warn("Purple fit unusable; cannot report based on viral depth");
        }

        LOGGER.info("Calculating reporting statuses");
        List<VirusReport> reports = VirusReporter.report(
                profiles, integrationCounts, mVirusReference::oncologyGroupInfo, expectedViralDepthPerCopy);

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
