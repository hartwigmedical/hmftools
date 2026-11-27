package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.ALT_CONTIG_SUFFIX;
import static com.hartwig.hmftools.tars.common.TarsConstants.APP_NAME;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.ContigSidecar;
import com.hartwig.hmftools.tars.fasta.TranscriptContigBuilder.TranscriptContigResult;

// builds one contig sequence per chromosome (trans1 + N-spacer + trans2 + ...) plus a sidecar file for liftback
public class SpliceFastaBuilder
{
    // N-padding between transcripts on the same alt contig
    private static final int SPACER_LENGTH = 500;
    private static final String TRANSCRIPT_CONTIGS_FILE_ID = ".fasta";
    private static final String CONTIG_MAPPINGS_FILE_ID = ".rna_contigs_mappings.tsv";

    private final SpliceFastaConfig mConfig;
    private final EnsemblDataCache mEnsemblDataCache;

    private int mGenesProcessed;
    private int mSkippedSingleExon;

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
            throw new IllegalStateException(String.format("failed to load ref genome: %s", mConfig.RefGenomeFile));
        }

        TARS_LOGGER.info("loading ensembl cache from {}", mConfig.EnsemblDataDir);
        mEnsemblDataCache.setRequiredData(true, false, false, false); // exons only
        mEnsemblDataCache.load(false);

        String fastaFile = mConfig.formFilename(TRANSCRIPT_CONTIGS_FILE_ID);
        String mappingsFile = mConfig.formFilename(CONTIG_MAPPINGS_FILE_ID);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(refGenome);

        List<ContigEntry> contigEntries = new ArrayList<>();

        int contigsWritten = 0;

        try(BufferedWriter fastaWriter = createBufferedWriter(fastaFile))
        {
            for(Map.Entry<String, List<GeneData>> chrEntry : mEnsemblDataCache.getChrGeneDataMap().entrySet())
            {
                List<TranscriptContigResult> chromosomeTranscripts = buildChromosomeTranscripts(
                        builder, chrEntry.getValue(), contigEntries);

                if(!chromosomeTranscripts.isEmpty())
                {
                    String altContig = chrEntry.getKey() + ALT_CONTIG_SUFFIX;
                    String sequence = packChromosomeContig(altContig, chromosomeTranscripts, contigEntries, SPACER_LENGTH);
                    writeFastaContig(fastaWriter, altContig, sequence);
                    contigsWritten += chromosomeTranscripts.size();
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException("failed to write transcript contigs FASTA", e);
        }

        TARS_LOGGER.info(
                "wrote {} transcript contigs to {} (skipped {} single-exon)",
                contigsWritten, fastaFile, mSkippedSingleExon);

        ContigSidecar.write(mappingsFile, contigEntries);

        TARS_LOGGER.info("SpliceFastaBuilder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<TranscriptContigResult> buildChromosomeTranscripts(
            final TranscriptContigBuilder builder, final List<GeneData> genes, final List<ContigEntry> contigEntries)
    {
        List<TranscriptContigResult> chromosomeTranscripts = new ArrayList<>();

        for(GeneData gene : genes)
        {
            ++mGenesProcessed;

            List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(gene.GeneId);
            if(transcripts == null || transcripts.isEmpty())
            {
                throw new IllegalStateException(String.format(
                        "ensembl cache returned no transcripts for gene %s", gene.GeneId));
            }

            for(TranscriptData transcript : transcripts)
            {
                // single-exon transcripts have no junctions, so no contig; record the exon as an annotation-only
                // row so the liftback exon index covers the full ensembl cache
                if(transcript.exons().size() < 2)
                {
                    ++mSkippedSingleExon;
                    contigEntries.add(ContigEntry.annotationOnly(
                            gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand,
                            TranscriptContigBuilder.sortedExonSpans(transcript)));
                }
                else
                {
                    TranscriptContigResult result = builder.build(gene, transcript);

                    // an all-N sequence cannot be aligned to, so no contig; keep its exons as annotation only
                    if(result.sequence().chars().allMatch(c -> c == 'N' || c == 'n'))
                    {
                        contigEntries.add(ContigEntry.annotationOnly(
                                gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand,
                                result.exonSpans()));
                    }
                    else
                    {
                        chromosomeTranscripts.add(result);
                    }
                }
            }

            if(mGenesProcessed % 1000 == 0)
            {
                TARS_LOGGER.debug("processed {} genes", mGenesProcessed);
            }
        }

        return chromosomeTranscripts;
    }

    static String packChromosomeContig(
            final String altContig, final List<TranscriptContigResult> transcripts,
            final List<ContigEntry> entries, final int spacerLength)
    {
        String spacer = "N".repeat(spacerLength);
        StringBuilder sequence = new StringBuilder();

        for(TranscriptContigResult transcript : transcripts)
        {
            if(sequence.length() > 0)
            {
                sequence.append(spacer);
            }

            int contigStart = sequence.length() + 1;
            sequence.append(transcript.sequence());
            int contigEnd = sequence.length();

            entries.add(new ContigEntry(
                    new ChrBaseRegion(altContig, contigStart, contigEnd),
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

        try
        {
            SpliceFastaBuilder builder = new SpliceFastaBuilder(configBuilder);
            builder.run();
        }
        catch(RuntimeException e)
        {
            TARS_LOGGER.error("Runtime error", e);
            System.exit(1);
        }
    }
}
