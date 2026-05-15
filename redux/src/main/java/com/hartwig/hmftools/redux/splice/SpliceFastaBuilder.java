package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.CONTIG_SIDECAR_FILE_ID;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.TRANSCRIPT_CONTIGS_FILE_ID;

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

public class SpliceFastaBuilder
{
    private final SpliceFastaConfig mConfig;
    private final EnsemblDataCache mEnsemblDataCache;

    public SpliceFastaBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceFastaConfig(configBuilder);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        long startTimeMs = System.currentTimeMillis();

        RefGenomeSource refGenome = RefGenomeSource.loadRefGenome(mConfig.RefGenomeFile);
        if(refGenome == null)
        {
            RD_LOGGER.error("failed to load ref genome: {}", mConfig.RefGenomeFile);
            System.exit(1);
        }

        RD_LOGGER.info("loading ensembl cache from {}", mConfig.EnsemblDataDir);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        Map<String, List<GeneData>> chrGeneMap = mEnsemblDataCache.getChrGeneDataMap();

        String fastaFile = mConfig.formFilename(TRANSCRIPT_CONTIGS_FILE_ID);
        String sidecarFile = mConfig.formFilename(CONTIG_SIDECAR_FILE_ID);

        TranscriptContigBuilder builder = new TranscriptContigBuilder(refGenome);

        List<ContigEntry> contigEntries = new ArrayList<>();

        int genesProcessed = 0;
        int contigsWritten = 0;
        int skippedSingleExon = 0;
        int skippedAllN = 0;

        try(BufferedWriter fastaWriter = createBufferedWriter(fastaFile))
        {
            for(List<GeneData> genes : chrGeneMap.values())
            {
                for(GeneData gene : genes)
                {
                    ++genesProcessed;

                    List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(gene.GeneId);
                    // TODO: do i really need this here? CHECK if ensembledatacache already covers this
                    if(transcripts == null || transcripts.isEmpty())
                        continue;

                    for(TranscriptData transcript : transcripts)
                    {
                        // single-exon transcripts have no junctions — bwa-mem2 already aligns them fine against genomic ref
                        if(transcript.exons().size() < 2)
                        {
                            ++skippedSingleExon;
                            continue;
                        }

                        TranscriptContigBuilder.TranscriptContigResult result = builder.build(gene, transcript);
                        if(result == null)
                            continue;

                        // exonic positions all masked to N in the ref — contig is alignment-useless
                        if(isAllN(result.Sequence))
                        {
                            ++skippedAllN;
                            continue;
                        }

                        writeFastaContig(fastaWriter, result.Entry.ContigName, result.Sequence);
                        contigEntries.add(result.Entry);
                        ++contigsWritten;
                    }

                    if(genesProcessed % 1000 == 0)
                        RD_LOGGER.debug("processed {} genes, {} contigs", genesProcessed, contigsWritten);
                }
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to write transcript contigs FASTA: {}", e.toString());
            System.exit(1);
        }

        RD_LOGGER.info("wrote {} transcript contigs to {} (skipped {} single-exon, {} all-N)",
                contigsWritten, fastaFile, skippedSingleExon, skippedAllN);

        ContigSidecar.write(sidecarFile, contigEntries);

        RD_LOGGER.info("next: cat <ref.fasta> {} > ref_genome_<ver>_rna_contigs.fasta && bwa-mem2 index ref_genome_<ver>_rna_contigs.fasta",
                fastaFile);

        RD_LOGGER.info("SpliceFastaBuilder complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static boolean isAllN(final String sequence)
    {
        for(int i = 0; i < sequence.length(); ++i)
        {
            char c = sequence.charAt(i);
            if(c != 'N' && c != 'n')
                return false;
        }
        return true;
    }

    private static void writeFastaContig(final BufferedWriter writer, final String contigName, final String sequence) throws IOException
    {
        // TODO: confirm this with Charles
        // unwrapped: matches BlastnRunner.writeBlastFasta. bwa-mem2 / samtools accept any line length.
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
