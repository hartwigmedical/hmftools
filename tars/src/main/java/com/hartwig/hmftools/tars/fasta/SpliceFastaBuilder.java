package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.CONTIG_MAPPINGS_FILE_ID;
import static com.hartwig.hmftools.tars.common.TarsConstants.TRANSCRIPT_CONTIGS_FILE_ID;
import static com.hartwig.hmftools.tars.common.TarsConstants.APP_NAME;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.tars.common.ContigEntry;
import com.hartwig.hmftools.tars.common.ContigSidecar;

public class SpliceFastaBuilder
{
    // N-padding between transcripts on the same alt contig
    private static final int SPACER_LENGTH = 500;

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
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        Map<String, List<GeneData>> chrGeneMap = mEnsemblDataCache.getChrGeneDataMap();

        String fastaFile = mConfig.formFilename(TRANSCRIPT_CONTIGS_FILE_ID);
        String mappingsFile = mConfig.formFilename(CONTIG_MAPPINGS_FILE_ID);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(refGenome);
        AltContigPacker packer = new AltContigPacker(SPACER_LENGTH);

        List<ContigEntry> contigEntries = new ArrayList<>();

        int genesProcessed = 0;
        int contigsWritten = 0;
        int skippedSingleExon = 0;
        int skippedAllN = 0;

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
                        throw new IllegalStateException("ensembl cache returned no transcripts for gene " + gene.GeneId);
                    }

                    for(TranscriptData transcript : transcripts)
                    {
                        // single-exon transcripts have no junctions, so no contig; still record their exon as an
                        // annotation-only row so the liftback exon index matches the full ensembl cache.
                        if(transcript.exons().size() < 2)
                        {
                            ++skippedSingleExon;
                            contigEntries.add(annotationEntry(gene, transcript, exonSpansOf(transcript)));
                            continue;
                        }

                        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene, transcript);
                        if(result == null)
                            continue;

                        // all-N contig is alignment-useless, so no contig; keep its exons as annotation only.
                        if(isAllN(result.sequence()))
                        {
                            ++skippedAllN;
                            contigEntries.add(annotationEntry(gene, transcript, result.exonSpans()));
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

                AltContigPacker.PackResult packed = packer.pack(chromosome, chromosomeTranscripts);
                writeFastaContig(fastaWriter, packed.altContig(), packed.sequence());
                contigEntries.addAll(packed.entries());
                contigsWritten += packed.entries().size();
            }
        }
        catch(IOException e)
        {
            TARS_LOGGER.error("failed to write transcript contigs FASTA: {}", e.toString());
            System.exit(1);
        }

        TARS_LOGGER.info("wrote {} transcript contigs to {} (skipped {} single-exon, {} all-N)",
                contigsWritten, fastaFile, skippedSingleExon, skippedAllN);

        ContigSidecar.write(mappingsFile, contigEntries);

        TARS_LOGGER.info(
                "next: cat <ref.fasta> {} > ref_genome_<ver>_rna_contigs.fasta"
                        + " && bwa-mem2 index ref_genome_<ver>_rna_contigs.fasta",
                fastaFile);

        TARS_LOGGER.info("SpliceFastaBuilder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static ContigEntry annotationEntry(final GeneData gene, final TranscriptData transcript, final List<BaseRegion> exonSpans)
    {
        return ContigEntry.annotationOnly(
                gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand, exonSpans);
    }

    private static List<BaseRegion> exonSpansOf(final TranscriptData transcript)
    {
        List<ExonData> ordered = new ArrayList<>(transcript.exons());
        ordered.sort(Comparator.comparingInt(exon -> exon.Start));

        List<BaseRegion> spans = new ArrayList<>(ordered.size());
        for(ExonData exon : ordered)
        {
            spans.add(new BaseRegion(exon.Start, exon.End));
        }
        return spans;
    }

    private static boolean isAllN(final String sequence)
    {
        for(int i = 0; i < sequence.length(); ++i)
        {
            char c = sequence.charAt(i);
            if(c != 'N' && c != 'n')
            {
                return false;
            }
        }
        return true;
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
        SpliceFastaConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SpliceFastaBuilder builder = new SpliceFastaBuilder(configBuilder);
        builder.run();
    }
}
