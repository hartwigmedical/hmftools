package com.hartwig.hmftools.viridian.app;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.ALIGNED_READ_BAM_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.APP_NAME;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.CANDIDATE_READ_FASTA_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.CONTIG_INFO_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.INTEGRATIONS_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.PAIRWISE_MARGINS_TSV_SUFFIX;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_MIN_SOFT_CLIP_BASES_DEFAULT;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_REF_CONTIGS;

import java.io.File;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.viridian.common.UserInputError;
import com.hartwig.hmftools.viridian.detection.contig_support.ContigSupport;
import com.hartwig.hmftools.viridian.detection.contig_support.ContigSupportCalculator;
import com.hartwig.hmftools.viridian.detection.read_align.ViralReadAligner;
import com.hartwig.hmftools.viridian.detection.read_align.ViralReadAlignments;
import com.hartwig.hmftools.viridian.detection.read_extract.CandidateReadExtractor;
import com.hartwig.hmftools.viridian.detection.read_extract.CandidateReadFilter;
import com.hartwig.hmftools.viridian.integration.seq_align.ViralSequenceAligner;
import com.hartwig.hmftools.viridian.integration.seq_align.ViralSequenceAlignment;
import com.hartwig.hmftools.viridian.integration.variant_extract.CandidateIntegration;
import com.hartwig.hmftools.viridian.integration.variant_extract.CandidateIntegrationExtractor;
import com.hartwig.hmftools.viridian.reference.ViralContig;
import com.hartwig.hmftools.viridian.reference.ViralReference;
import com.hartwig.hmftools.viridian.selection.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.selection.OncologyGroupResolution;
import com.hartwig.hmftools.viridian.selection.PairwiseMargins;
import com.hartwig.hmftools.viridian.selection.RepresentativeContigSelector;

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
        LOGGER.info("Viral reference data loaded");
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

        selectRepresentativeViralContigs(viralReadAlignments, viralContigSupports);

        callHostIntegrations();

        LOGGER.info("{} complete, mins({})", APP_NAME, runTimeMinsStr(startTimeMs));
    }

    private ViralReadAlignments getViralReadAlignments()
    {
        String viralReadBamFile = outputFile(ALIGNED_READ_BAM_SUFFIX);
        // Alignment is pretty slow, so allow reusing the cached BAM for a rerun.
        if(!canReuseExistingFile(mConfig.reuseReadsBam(), viralReadBamFile, "aligned read BAM"))
        {
            alignCandidateReadsToViralContigs(viralReadBamFile);
        }
        return ViralReadAlignments.load(viralReadBamFile, mViralReference);
    }

    private String getCandidateReads()
    {
        String candidateFastaFile = outputFile(CANDIDATE_READ_FASTA_SUFFIX);
        // Read extraction is very slow for large samples, so allow reusing the cached FASTA for a rerun.
        if(!canReuseExistingFile(mConfig.reuseReadsFasta(), candidateFastaFile, "candidate read FASTA"))
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
        mCandidateExtractor.extractToFasta(mConfig.tumorBam(), candidateFastaFile);
        LOGGER.info("Candidate read extraction complete");
    }

    // Align potentially viral reads to all virus genomes, so we can decide which viruses are present.
    private void alignCandidateReadsToViralContigs(String viralReadBamFile)
    {
        String candidateReadFasta = getCandidateReads();

        LOGGER.info("Aligning candidate reads to viral genomes");
        ViralReadAligner viralReadAligner = ViralReadAligner.create(
                mViralReference, mConfig.viralBwaIndexImage(), mConfig.threads(), mConfig.alignmentBatchSize());
        viralReadAligner.align(candidateReadFasta, viralReadBamFile);
        LOGGER.info("Candidate read alignment complete");
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
        logRepresentativeContigSelections(selections);

        OutputWriter.writeContigInfo(outputFile(CONTIG_INFO_TSV_SUFFIX), selections);
        if(mConfig.verboseOutput())
        {
            OutputWriter.writePairwiseMargins(outputFile(PAIRWISE_MARGINS_TSV_SUFFIX), pairwiseMargins, selections);
        }
        LOGGER.info("Selecting representative contigs complete");
    }

    private static void logRepresentativeContigSelections(List<OncologyGroupRepresentativeSelection> selections)
    {
        for(OncologyGroupRepresentativeSelection selection : selections)
        {
            ViralContig representative = selection.representative();
            if(representative != null)
            {
                LOGGER.info("oncologyGroup({}) representative({})", selection.oncologyGroup(), representative);
            }
            else if(selection.resolution() == OncologyGroupResolution.UNRESOLVED)
            {
                LOGGER.warn("oncologyGroup({}) unresolved({})", selection.oncologyGroup(), selection.outcome());
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

        LOGGER.info("Extracting integration candidates from ESVEE VCF: {}", esveeVcf);
        List<CandidateIntegration> candidates = new CandidateIntegrationExtractor(mConfig.sampleId()).extract(esveeVcf);

        LOGGER.info("Aligning {} candidate viral integration sequences to the viral reference", candidates.size());
        Map<CandidateIntegration, ViralSequenceAlignment> alignments = alignCandidateIntegrations(candidates);
        LOGGER.info("Aligned {} of {} inserted sequences to a viral contig", alignments.size(), candidates.size());

        OutputWriter.writeIntegrations(outputFile(INTEGRATIONS_TSV_SUFFIX), candidates, alignments);
        LOGGER.info("Integration site calling complete");
    }

    private Map<CandidateIntegration, ViralSequenceAlignment> alignCandidateIntegrations(List<CandidateIntegration> candidates)
    {
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
