package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConstants.APP_NAME;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.ContigSidecar;

// builds one contig sequence per chromosome (ex: trans1 + xN + trans2 + ...) with sidecar file for liftback
public class SpliceFastaBuilder
{
    // N-padding between transcripts on the same alt contig
    private static final int SPACER_LENGTH = 500;
    private static final String TRANSCRIPT_CONTIGS_FILE_ID = ".fasta";
    private static final String CONTIG_MAPPINGS_FILE_ID = ".rna_contigs_mappings.tsv";

    private final SpliceFastaConfig mConfig;
    private final EnsemblDataCache mEnsemblDataCache;

    public SpliceFastaBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceFastaConfig(configBuilder);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        RefGenomeSource refGenome = RefGenomeSource.loadRefGenome(mConfig.RefGenomeFile);
        if(refGenome == null)
        {
            TARS_LOGGER.error("failed to load ref genome: {}", mConfig.RefGenomeFile);
            System.exit(1);
        }

        TARS_LOGGER.info("loading ensembl cache from {}", mConfig.EnsemblDataDir);
        mEnsemblDataCache.setRequiredData(true, false, false, false); // exons only
        mEnsemblDataCache.load(false);

        Map<String, List<GeneData>> chrGeneMap = mEnsemblDataCache.getChrGeneDataMap();

        String fastaFile = mConfig.formFilename(TRANSCRIPT_CONTIGS_FILE_ID);
        String mappingsFile = mConfig.formFilename(CONTIG_MAPPINGS_FILE_ID);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(refGenome);

        List<ContigEntry> contigEntries = new ArrayList<>();

        int genesProcessed = 0;
        int contigsWritten = 0;
        int skippedSingleExon = 0;

        try(BufferedWriter fastaWriter = createBufferedWriter(fastaFile))
        {
            for(Map.Entry<String, List<GeneData>> chrEntry : chrGeneMap.entrySet())
            {
                String chromosome = chrEntry.getKey();
                List<TranscriptContigBuilder.TranscriptContigResult> chromosomeTranscripts = new ArrayList<>();

                for(GeneData gene : chrEntry.getValue())
                {
                    ++genesProcessed;

                    List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(gene.GeneId);
                    if(transcripts == null || transcripts.isEmpty())
                    {
                        TARS_LOGGER.error("ensembl cache returned no transcripts for gene {}", gene.GeneId);
                        System.exit(1);
                    }

                    for(TranscriptData transcript : transcripts)
                    {
                        // single-exon transcripts have no junctions, so no contig; still record their exon as an
                        // annotation-only row so the liftback exon index matches the full ensembl cache.
                        if(transcript.exons().size() < 2)
                        {
                            ++skippedSingleExon;
                            contigEntries.add(ContigEntry.annotationOnly(
                                    gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand,
                                    TranscriptContigBuilder.sortedExonSpans(transcript)));
                            continue;
                        }

                        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene, transcript);
                        if(result == null)
                            continue;

                        // all-N contig is alignment-useless, so no contig; keep its exons as annotation only.
                        if(result.sequence().chars().allMatch(c -> c == 'N' || c == 'n'))
                        {
                            contigEntries.add(ContigEntry.annotationOnly(
                                    gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand,
                                    result.exonSpans()));
                            continue;
                        }

                        chromosomeTranscripts.add(result);
                    }

                    if(genesProcessed % 1000 == 0)
                    {
                        TARS_LOGGER.debug("processed {} genes", genesProcessed);
                    }
                }

                if(chromosomeTranscripts.isEmpty())
                    continue;

                String altContig = chromosome + ALT_CONTIG_SUFFIX;
                String sequence = packChromosomeContig(altContig, chromosomeTranscripts, contigEntries, SPACER_LENGTH);
                writeFastaContig(fastaWriter, altContig, sequence);
                contigsWritten += chromosomeTranscripts.size();
            }
        }
        catch(IOException e)
        {
            TARS_LOGGER.error("failed to write transcript contigs FASTA: {}", e.toString());
            System.exit(1);
        }

        TARS_LOGGER.info(
                "wrote {} transcript contigs to {} (skipped {} single-exon)",
                contigsWritten, fastaFile, skippedSingleExon);

        ContigSidecar.write(mappingsFile, contigEntries);

        TARS_LOGGER.info("SpliceFastaBuilder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // builds the chromosome's alt-contig sequence, recording where each transcript lands in it
    static String packChromosomeContig(
            final String altContig, final List<TranscriptContigBuilder.TranscriptContigResult> transcripts,
            final List<ContigEntry> entries, final int spacerLength)
    {
        String spacer = "N".repeat(spacerLength);
        StringBuilder sequence = new StringBuilder();

        for(TranscriptContigBuilder.TranscriptContigResult transcript : transcripts)
        {
            if(sequence.length() > 0)
            {
                sequence.append(spacer);
            }

            int contigStart = sequence.length() + 1;
            sequence.append(transcript.sequence());
            int contigEnd = sequence.length();

            entries.add(new ContigEntry(
                    altContig, contigStart, contigEnd,
                    transcript.geneId(), transcript.geneName(), transcript.transName(), transcript.chromosome(),
                    transcript.strand(), transcript.exonSpans()));
        }

        return sequence.toString();
    }

    private static void writeFastaContig(final BufferedWriter writer, final String contigName, final String sequence) throws IOException
    {
        writer.write(">" + contigName);
        writer.newLine();
        writer.write(sequence);
        writer.newLine();
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceFastaConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SpliceFastaBuilder builder = new SpliceFastaBuilder(configBuilder);
        builder.run();
    }
}
