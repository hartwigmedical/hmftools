package com.hartwig.hmftools.viridian.app;

import static java.lang.System.exit;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.ALL_ALIGNMENTS_BAM_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.APP_NAME;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.CANDIDATES_FASTA_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.CONTIG_INFO_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATIONS_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.PAIRWISE_MARGINS_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_REF_CONTIGS;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.viridian.common.UserInputError;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAligner;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignments;
import com.hartwig.hmftools.viridian.detection.assign.RepresentativeReadAssigner;
import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.ContigStatsCalculator;
import com.hartwig.hmftools.viridian.detection.extract.CandidateReadExtractor;
import com.hartwig.hmftools.viridian.detection.extract.CandidateReadFilter;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupResolution;
import com.hartwig.hmftools.viridian.detection.select.PairwiseMargins;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeContigSelector;
import com.hartwig.hmftools.viridian.detection.support.ContigSupport;
import com.hartwig.hmftools.viridian.detection.support.ContigSupportCalculator;
import com.hartwig.hmftools.viridian.integration.align.ViralSequenceAligner;
import com.hartwig.hmftools.viridian.integration.align.ViralSequenceAlignment;
import com.hartwig.hmftools.viridian.integration.extract.CandidateIntegration;
import com.hartwig.hmftools.viridian.integration.extract.CandidateIntegrationExtractor;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reference.ViralReference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViridianApplication
{
    private final ViridianConfig mConfig;
    private final ViralReference mViralReference;

    private static final Logger LOGGER = LogManager.getLogger(ViridianApplication.class);

    public ViridianApplication(ViridianConfig config)
    {
        mConfig = config;

        LOGGER.info("Loading viral reference data");
        mViralReference = ViralReference.load(config.viralRefFile(), config.viralRefInfoFile());
    }

    public void run() throws IOException
    {
        LOGGER.info("Starting {}", APP_NAME);
        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        checkCreateOutputDir(mConfig.outputDir());

        BwaMemAligner.initLibrary(mConfig.bwaLibPath());

        ViralReadAlignments viralReadAlignments = getViralReadAlignments();

        List<ContigSupport> viralContigSupports = computeViralContigSupport(viralReadAlignments);

        List<OncologyGroupRepresentativeSelection> selections =
                selectRepresentativeViralContigs(viralReadAlignments, viralContigSupports);

        assignReadsToRepresentatives(viralReadAlignments, selections);

        callHostIntegrations();

        LOGGER.info("{} complete, mins({})", APP_NAME, runTimeMinsStr(startTimeMs));
    }

    private ViralReadAlignments getViralReadAlignments()
    {
        String allAlignmentsBamFile = outputFile(ALL_ALIGNMENTS_BAM_SUFFIX);
        // Alignment is pretty slow, so allow reusing the cached BAM for a rerun.
        if(!canReuseExistingFile(mConfig.reuseAlignments(), allAlignmentsBamFile, "all-alignments BAM"))
        {
            alignCandidateReadsToViralContigs(allAlignmentsBamFile);
        }
        return ViralReadAlignments.load(allAlignmentsBamFile, mViralReference);
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
        LOGGER.info("Extracting candidate viral reads from tumor BAM");
        CandidateReadFilter candidateFilter = new CandidateReadFilter(VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT, VIRAL_REF_CONTIGS);
        CandidateReadExtractor mCandidateExtractor = new CandidateReadExtractor(
                mConfig.refGenomeFile(), candidateFilter, mConfig.threads());
        mCandidateExtractor.extract(mConfig.tumorBam(), candidateFastaFile);
    }

    // Align potentially viral reads to all virus genomes, so we can decide which viruses are present.
    private void alignCandidateReadsToViralContigs(String allAlignmentsBamFile)
    {
        String candidateReadFasta = getCandidateReads();

        LOGGER.info("Aligning candidate reads to viral genomes");
        ViralReadAligner viralReadAligner = ViralReadAligner.create(
                mViralReference, mConfig.viralBwaIndexImage(), mConfig.threads(), mConfig.alignmentBatchSize());
        viralReadAligner.align(candidateReadFasta, allAlignmentsBamFile);
    }

    // Compute support information for each virus genome and decide which genomes may be present.
    private List<ContigSupport> computeViralContigSupport(ViralReadAlignments viralReadAlignments)
    {
        LOGGER.info("Computing per-contig support");
        Map<ViralContig, ContigStats> contigStats = ContigStatsCalculator.calculate(viralReadAlignments);
        return new ContigSupportCalculator().compute(viralReadAlignments, contigStats);
    }

    // For each oncology group (group of virus strains at interesting taxonomy granularity), select 1 viral genome which best represents
    // the virus present in the sample.
    private List<OncologyGroupRepresentativeSelection> selectRepresentativeViralContigs(
            ViralReadAlignments viralReadAlignments, List<ContigSupport> viralContigSupports)
    {
        LOGGER.info("Selecting representative contig per oncology group");
        PairwiseMargins pairwiseMargins = PairwiseMargins.from(viralReadAlignments);
        List<OncologyGroupRepresentativeSelection> selections = RepresentativeContigSelector.select(
                viralContigSupports, pairwiseMargins, viralReadAlignments.readCountsByOncologyGroup());
        logRepresentativeContigSelections(selections);

        LOGGER.info("Writing contig info output");
        OutputWriter.writeContigInfo(outputFile(CONTIG_INFO_TSV_SUFFIX), selections);
        if(mConfig.verboseOutput())
        {
            LOGGER.info("Writing pairwise margins output");
            OutputWriter.writePairwiseMargins(outputFile(PAIRWISE_MARGINS_TSV_SUFFIX), pairwiseMargins, selections);
        }

        return selections;
    }

    // For each oncology group, assign the read alignments to only the selected representative viral genome.
    // This is simply filtering down to each read's alignment to that genome contig (or nothing, if it didn't align there at all).
    private void assignReadsToRepresentatives(
            ViralReadAlignments viralReadAlignments, List<OncologyGroupRepresentativeSelection> selections)
    {
        LOGGER.info("Assigning reads to representative contigs");
        Set<ViralContig> representatives = selections.stream()
                .map(OncologyGroupRepresentativeSelection::representative)
                .filter(Objects::nonNull)
                .collect(toSet());
        RepresentativeReadAssigner.assign(
                viralReadAlignments, representatives, outputFile(ALL_ALIGNMENTS_BAM_SUFFIX),
                outputFile(REPRESENTATIVE_ALIGNMENTS_BAM_SUFFIX));
    }

    private static void logRepresentativeContigSelections(List<OncologyGroupRepresentativeSelection> selections)
    {
        for(OncologyGroupRepresentativeSelection selection : selections)
        {
            ViralContig representative = selection.representative();
            if(representative != null)
            {
                LOGGER.debug("oncologyGroup({}) representative({})", selection.oncologyGroup(), representative);
            }
            else if(selection.resolution() == OncologyGroupResolution.UNRESOLVED)
            {
                LOGGER.debug("oncologyGroup({}) unresolved({})", selection.oncologyGroup(), selection.outcome());
            }
        }
    }

    // Find where a virus inserted itself into the host genome, from the SVs ESVEE called.
    private void callHostIntegrations()
    {
        String esveeVcf = mConfig.esveeUnfilteredVcf();
        if(esveeVcf == null)
        {
            LOGGER.info("ESVEE VCF not specified; skipping integration site calling");
            return;
        }

        LOGGER.info("Extracting integration candidates from ESVEE VCF");
        List<CandidateIntegration> candidates = new CandidateIntegrationExtractor(mConfig.sampleId()).extract(esveeVcf);

        Map<CandidateIntegration, ViralSequenceAlignment> alignments = alignCandidateIntegrations(candidates);

        LOGGER.info("Writing integrations output");
        OutputWriter.writeIntegrations(outputFile(INTEGRATIONS_TSV_SUFFIX), candidates, alignments);
    }

    private Map<CandidateIntegration, ViralSequenceAlignment> alignCandidateIntegrations(List<CandidateIntegration> candidates)
    {
        LOGGER.info("Aligning candidate viral integration sequences to viral reference");

        List<String> insertSequences = candidates.stream().map(CandidateIntegration::insertSequence).toList();
        ViralSequenceAligner aligner = ViralSequenceAligner.create(
                mViralReference, mConfig.viralBwaIndexImage(), mConfig.threads());
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
