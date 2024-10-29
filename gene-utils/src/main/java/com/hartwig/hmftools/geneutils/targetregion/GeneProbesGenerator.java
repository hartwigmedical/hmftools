package com.hartwig.hmftools.geneutils.targetregion;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.blastn.BlastnMatch;
import com.hartwig.hmftools.common.blastn.BlastnRunner;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class GeneProbesGenerator
{
    private static class GeneNameTranscriptId
    {
        public final String GeneName;
        public final String EnsemblTranscriptId;

        public GeneNameTranscriptId(final String geneName, final String ensemblTranscriptId)
        {
            GeneName = geneName;
            EnsemblTranscriptId = ensemblTranscriptId;
        }
    }

    private static final double PROBE_GC_MIN = 0.35;
    private static final double PROBE_GC_MAX = 0.55;
    private static final double MAX_PROBE_SUM_BLASTN_BITSCORE = 2500;

    private static final int PROBE_LENGTH = 120;

    private static final int MAX_CANDIDATE_PROBES = 8;
    private static final int MIN_INTRON_LENGTH = 3000;
    private static final int LONG_INTRON_LENGTH = 5000;

    private static final int MAX_EXONS_TO_ADD_INTRON = 19;

    private static final int FLANKING_DISTANCE = 1000;
    private static final int CANDIDATE_REGION_SIZE = PROBE_LENGTH * MAX_CANDIDATE_PROBES;

    private static final int BLASTN_WORD_SIZE = 15;

    private static final int MIN_BLAST_ALIGNMENT_LENGTH = 30;

    private static final String GENE_TRANSCRIPT_FILE = "gene_transcript_file";
    private static final String OUTPUT_PREFIX = "output_prefix";

    private static final String BLAST = "blast";

    private static final String BLAST_DB = "blast_db";

    private final String mOutputPrefix;

    private final String mOutputDir;

    private final int mThreads;

    private final String mBlast;

    private final String mBlastDb;

    private final RefGenomeInterface mRefGenome;

    private final EnsemblDataCache mEnsemblDataCache;

    private final List<GeneNameTranscriptId> mGeneNameTranscriptIds = new ArrayList<>();

    public GeneProbesGenerator(final ConfigBuilder configBuilder) throws FileNotFoundException
    {
        ConfigUtils.setLogLevel(configBuilder);
        mOutputPrefix = configBuilder.getValue(OUTPUT_PREFIX);
        mOutputDir = parseOutputDir(configBuilder);
        mThreads = parseThreads(configBuilder);
        mBlast = configBuilder.getValue(BLAST);
        mBlastDb = configBuilder.getValue(BLAST_DB);

        final String refGenomeFile = configBuilder.getValue(REF_GENOME);
        mRefGenome = new RefGenomeSource(new IndexedFastaSequenceFile(new File(refGenomeFile)));

        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        try (DelimFileReader reader = new DelimFileReader(configBuilder.getValue(GENE_TRANSCRIPT_FILE)))
        {
            reader.stream().forEach(row -> mGeneNameTranscriptIds.add(
                    new GeneNameTranscriptId(row.get(FLD_GENE_NAME), row.get(FLD_TRANS_ID))));
        }
    }

    public int run()
    {
        List<TargetedGene> targetedGenes = new ArrayList<>();

        // from the gene, find all the gene regions
        for(GeneNameTranscriptId geneNameTranscriptId : mGeneNameTranscriptIds)
        {
            String geneName = geneNameTranscriptId.GeneName;
            String transcriptId = geneNameTranscriptId.EnsemblTranscriptId;

            GU_LOGGER.info("processing gene: {}, ensembl transcript id: {}", geneName, transcriptId);

            // get all the gene regions from the ensembl
            GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);

            if(geneData == null)
            {
                GU_LOGGER.error("transcript data for gene: {} not found", geneName);
                throw new RuntimeException(String.format("transcript data for gene: %s not found", geneName));
            }

            TranscriptData transcriptData = mEnsemblDataCache.getTranscriptData(geneData.GeneId, transcriptId);

            if(transcriptData == null)
            {
                GU_LOGGER.error("transcript id({}) for gene({}) not found", transcriptId, geneName);
                throw new RuntimeException(String.format("transcript id(%s) for gene(%s) not found", transcriptId, geneName));
            }

            TargetedGene targetedGene = new TargetedGene(geneData, transcriptData);
            targetedGenes.add(targetedGene);
        }

        for(TargetedGene targetedGene : targetedGenes)
        {
            populateTargetedGeneRegions(targetedGene);
            populateCandidateProbes(targetedGene);
        }

        runBlastnOnProbeCandidates(targetedGenes);

        // now choose probes for each region
        selectProbeCandidates(targetedGenes);

        // write the outputs
        List<TargetedGeneRegion> regions = targetedGenes.stream().flatMap(o -> o.getRegions().stream()).collect(Collectors.toList());
        GeneProbeCandidateFileWriter.write(mOutputDir, mOutputPrefix, regions);
        GeneProbeRegionFileWriter.write(mOutputDir, mOutputPrefix, regions);

        return 0;
    }

    /*
    | Region Type | Gene                  | Design                                                                                         |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | Coding      | All                   | - Cover the full coding region of each exon                                                    |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | UTR         | All                   | - Add a 120nt probe centred on each non-coding exon                                            |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | Upstream/   | All                   | - Create 8 candidate 120nt probes from 1kb to 2kb bases upstream of gene.                      |
    | Downstream  |                       | - Filter for  0.35<GC<0.55 AND SUM_BLASTN<1000                                                 |
    |             |                       | - Choose 1 probe per region with lowest SUM_BLASTN score.                                      |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | Intronic    | Genes with < 20 exons | - For each intron of >5kb, create 8 candidate120nt probes from 1kb-2kb from EACH flanking exon |
    |             |                       | - For each intron of 3kb-5kb, create 8 candidate 120nt probes centred on intron                |
    |             |                       | - Filter for 0.35<GC<0.55 AND SUM_BLASTN<1000                                                  |
    |             |                       | - Choose 1 probe  per region with lowest SUM_BLASTN score.                                     |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    */
    static void populateTargetedGeneRegions(TargetedGene targetedGene)
    {
        // first we create a region 1-2kb upstream of the gene
        // TODO: should this be before first exon or before the gene start?
        targetedGene.addRegion(targetedGene.getGeneData().forwardStrand() ? TargetedGeneRegion.Type.UP_STREAM : TargetedGeneRegion.Type.DOWN_STREAM,
                targetedGene.getGeneData().GeneStart - FLANKING_DISTANCE - CANDIDATE_REGION_SIZE,
                targetedGene.getGeneData().GeneStart - FLANKING_DISTANCE - 1);

        int lastExonEnd = -1;

        TranscriptData transcript = targetedGene.getTranscriptData();

        for(ExonData exonData : transcript.exons())
        {
            if(lastExonEnd != -1 && targetedGene.getTranscriptData().exons().size() <= MAX_EXONS_TO_ADD_INTRON)
            {
                int intronLength = exonData.Start - lastExonEnd;

                if(intronLength > MIN_INTRON_LENGTH)
                {
                    if(intronLength > LONG_INTRON_LENGTH)
                    {
                        // for long intron, we want to split into two 1k regions each flank the exons
                        targetedGene.addRegion(TargetedGeneRegion.Type.INTRONIC_LONG,
                                lastExonEnd + 1 + FLANKING_DISTANCE,
                                lastExonEnd + FLANKING_DISTANCE + CANDIDATE_REGION_SIZE);

                        targetedGene.addRegion(TargetedGeneRegion.Type.INTRONIC_LONG,
                                exonData.Start - FLANKING_DISTANCE - CANDIDATE_REGION_SIZE,
                                exonData.Start - FLANKING_DISTANCE - 1);
                    }
                    else
                    {
                        // for smaller intron we add a region in the middle 1k
                        int intronMid = (lastExonEnd + 1 + exonData.Start) / 2;
                        targetedGene.addRegion(TargetedGeneRegion.Type.INTRONIC_SHORT,
                                intronMid - CANDIDATE_REGION_SIZE / 2,
                                intronMid + CANDIDATE_REGION_SIZE / 2 - 1);
                    }
                }
            }

            // not all exons are coding, but I am not sure how to get non coding ones

            // from Junran's script, an exon is coding of it overlaps with coding start / coding end
            boolean isCoding = transcript.CodingStart != null &&
                               transcript.CodingEnd != null &&
                               (exonData.Start < transcript.CodingEnd && exonData.End > transcript.CodingStart);

            if(isCoding)
            {
                // we take the coding part
                targetedGene.addRegion(TargetedGeneRegion.Type.CODING,
                        Math.max(exonData.Start, transcript.CodingStart),
                        Math.min(exonData.End, transcript.CodingEnd));
            }
            else
            {
                // just one probe at the middle, so just the 120 bases
                int exonMid = (exonData.Start + exonData.End + 1) / 2;
                targetedGene.addRegion(TargetedGeneRegion.Type.UTR, exonMid - PROBE_LENGTH / 2, exonMid + PROBE_LENGTH / 2 - 1);
            }

            lastExonEnd = exonData.End;
        }

        targetedGene.addRegion(targetedGene.getGeneData().forwardStrand() ? TargetedGeneRegion.Type.DOWN_STREAM : TargetedGeneRegion.Type.UP_STREAM,
                targetedGene.getGeneData().GeneEnd + 1 + FLANKING_DISTANCE,
                targetedGene.getGeneData().GeneEnd + FLANKING_DISTANCE + CANDIDATE_REGION_SIZE);
    }

    // from the targeted gene regions, we find the probes
    void populateCandidateProbes(TargetedGene targetedGene)
    {
        for(TargetedGeneRegion region : targetedGene.getRegions())
        {
            populateCandidateProbes(region, mRefGenome);
        }
    }

    static void populateCandidateProbes(TargetedGeneRegion region, RefGenomeInterface refGenomeInterface)
    {
        switch(region.getType())
        {
            case CODING:
            case UTR:
                // use whole region, so do not add probe candidate
                break;
            case UP_STREAM:
            case DOWN_STREAM:
            case INTRONIC_LONG:
            case INTRONIC_SHORT:
                // create 8 candidate probes
                for(int i = 0; i < MAX_CANDIDATE_PROBES; ++i)
                {
                    int start = region.getStart() + i * PROBE_LENGTH;
                    int end = start + PROBE_LENGTH - 1; // end is inclusive
                    ProbeCandidate probe = createProbeCandidate(new ChrBaseRegion(region.getGene().getGeneData().Chromosome, start, end),
                            refGenomeInterface);
                    region.getProbeCandidates().add(probe);
                }
                break;
        }
    }

    static ProbeCandidate createProbeCandidate(ChrBaseRegion chrBaseRegion, RefGenomeInterface refGenomeInterface)
    {
        // get the sequence from the ref genome
        String sequence = refGenomeInterface.getBaseString(chrBaseRegion.chromosome(), chrBaseRegion.start(), chrBaseRegion.end());
        double gcContent = GcCalcs.calcGcPercent(sequence);

        return new ProbeCandidate(chrBaseRegion, sequence, gcContent);
    }

    public void runBlastnOnProbeCandidates(Collection<TargetedGene> targetedGeneList)
    {
        int nextProbeKey = 0;

        Map<Integer, ProbeCandidate> probeCandidateMap = new HashMap<>();

        // get all the probes into a map with a keyed sequence
        for(TargetedGene targetedGene : targetedGeneList)
        {
            for(TargetedGeneRegion targetedGeneRegion : targetedGene.getRegions())
            {
                for(ProbeCandidate probeCandidate : targetedGeneRegion.getProbeCandidates())
                {
                    if(checkGcContent(probeCandidate))
                    {
                        int probeKey = nextProbeKey++;
                        probeCandidateMap.put(probeKey, probeCandidate);
                    }
                }
            }
        }

        GU_LOGGER.info("running blastn on {} probe candidates", probeCandidateMap.size());

        // run blastn on those
        Multimap<Integer, BlastnMatch> result = new BlastnRunner.Builder()
                .withTask("megablast")
                .withPrefix(mOutputPrefix)
                .withBlastDir(mBlast)
                .withBlastDb(mBlastDb)
                .withOutputDir(mOutputDir)
                .withKeepOutput(true)
                .withWordSize(BLASTN_WORD_SIZE)
                .withSubjectBestHit(true)
                .withNumThreads(mThreads)
                .build()
                .run(probeCandidateMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().getSequence())));

        // now go back to process the results
        for(int k : probeCandidateMap.keySet())
        {
            // get the probe
            ProbeCandidate probeCandidate = probeCandidateMap.get(k);
            Collection<BlastnMatch> matches = result.get(k);

            if(matches.isEmpty())
            {
                // this probe is not good
                GU_LOGGER.info("probe({}) no blastn match found, probably low complexity sequence", probeCandidate.toString());
                continue;
            }

            double sumBitScore = 0;

            // now process all the matches and sum up the bit score, but remove the one with best match
            Optional<BlastnMatch> bestMatch = matches.stream()
                    .filter(GeneProbesGenerator::isPrimaryBlastnMatch)
                    .max(Comparator.comparing(BlastnMatch::getBitScore));

            if(bestMatch.isPresent())
            {
                sumBitScore -= bestMatch.get().getBitScore();
            }

            for(BlastnMatch m : matches)
            {
                if(m.getAlignmentLength() >= MIN_BLAST_ALIGNMENT_LENGTH && isPrimaryBlastnMatch(m))
                {
                    sumBitScore += m.getBitScore();
                }
            }

            probeCandidate.setSumBlastnBitScore(sumBitScore);
        }

        GU_LOGGER.info("finished blastn on {} probes", probeCandidateMap.size());
    }

    // write out the candidate probe for each region we are interested
    public void selectProbeCandidates(Collection<TargetedGene> targetedGeneList)
    {
        for(TargetedGene targetedGene : targetedGeneList)
        {
            GU_LOGGER.debug("gene: {}", targetedGene.getGeneData().GeneName);

            for(TargetedGeneRegion targetedGeneRegion : targetedGene.getRegions())
            {
                // select the best probe
                GU_LOGGER.debug("gene: {}, region: type({}) {}:{}-{}",
                        targetedGene.getGeneData().GeneName, targetedGeneRegion.getType(),
                        targetedGeneRegion.getChromosome(), targetedGeneRegion.getStart(), targetedGeneRegion.getEnd());

                ProbeCandidate selectedProbe = null;

                for(ProbeCandidate probeCandidate : targetedGeneRegion.getProbeCandidates())
                {
                    checkSetProbeCandidateFilter(probeCandidate);

                    if(!probeCandidate.passFilter())
                    {
                        continue;
                    }

                    if(!Double.isNaN(probeCandidate.getSumBlastnBitScore()) &&
                       (selectedProbe == null || probeCandidate.getSumBlastnBitScore() < selectedProbe.getSumBlastnBitScore()))
                    {
                        selectedProbe = probeCandidate;
                    }
                }

                if(selectedProbe != null)
                {
                    targetedGeneRegion.setSelectedProbe(selectedProbe);

                    // now log the selected probe
                    GU_LOGGER.debug("region: type({}) {}:{}-{} selected probe: {}-{}", targetedGeneRegion.getType(),
                            targetedGeneRegion.getChromosome(), targetedGeneRegion.getStart(), targetedGeneRegion.getEnd(),
                            selectedProbe.getStart(), selectedProbe.getEnd());
                }
            }
        }
    }

    public static boolean checkGcContent(ProbeCandidate probeCandidate)
    {
        return probeCandidate.getGcContent() > PROBE_GC_MIN && probeCandidate.getGcContent() < PROBE_GC_MAX;
    }

    public static void checkSetProbeCandidateFilter(ProbeCandidate probeCandidate)
    {
        if(probeCandidate.getGcContent() <= PROBE_GC_MIN)
        {
            probeCandidate.setFilterReason(String.format("gc <= %.2f", PROBE_GC_MIN));
        }
        else if(probeCandidate.getGcContent() >= PROBE_GC_MAX)
        {
            probeCandidate.setFilterReason(String.format("gc >= %.2f", PROBE_GC_MAX));
        }
        else if(probeCandidate.getSumBlastnBitScore() >= MAX_PROBE_SUM_BLASTN_BITSCORE)
        {
            probeCandidate.setFilterReason(String.format("sum(bitscore) >= %g", MAX_PROBE_SUM_BLASTN_BITSCORE));
        }
    }

    public static boolean isPrimaryBlastnMatch(BlastnMatch match)
    {
        return match.getSubjectTitle().contains("Primary Assembly") ||
                match.getSubjectTitle().contains("unlocalized genomic scaffold") ||
                match.getSubjectTitle().contains("unplaced genomic scaffold");
    }

    public static void main(@NotNull final String[] args) throws FileNotFoundException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(GENE_TRANSCRIPT_FILE, true, "Gene name file");
        addRefGenomeVersion(configBuilder);
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC);
        configBuilder.addConfigItem(OUTPUT_PREFIX, true, "prefix of output BED filename");
        addOutputDir(configBuilder);
        EnsemblDataCache.addEnsemblDir(configBuilder, true);
        configBuilder.addPath(BLAST, true, "Location of blast installation");
        configBuilder.addPath(BLAST_DB, true, "Location of blast database");
        addThreadOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GeneProbesGenerator geneProbesGenerator = new GeneProbesGenerator(configBuilder);
        System.exit(geneProbesGenerator.run());
    }
}
