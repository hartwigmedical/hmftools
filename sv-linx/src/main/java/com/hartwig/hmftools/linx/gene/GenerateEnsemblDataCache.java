package com.hartwig.hmftools.linx.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.gene.EnsemblDAO.ENSEMBL_TRANS_SPLICE_DATA_FILE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class GenerateEnsemblDataCache
{
    private static final String LOG_DEBUG = "log_debug";
    private static final String OUTPUT_DIR = "output_dir";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        writeEnsemblDataFiles(cmd);
    }

    public static void writeEnsemblDataFiles(final CommandLine cmd)
    {
        final String outputDir = cmd.getOptionValue(OUTPUT_DIR);

        final RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, String.valueOf(V37)));

        LNX_LOGGER.info("writing Ensembl data files to {}", outputDir);

        if(EnsemblDAO.hasDatabaseConfig(cmd))
        {
            EnsemblDAO ensemblDAO = new EnsemblDAO(cmd);

            if(!ensemblDAO.isValid())
            {
                LNX_LOGGER.info("invalid Ensembl DAO");
                return;
            }

            ensemblDAO.writeDataCacheFiles(outputDir);
        }

        LNX_LOGGER.debug("reloading transcript data to generate splice acceptor positions");

        // create the transcript splice acceptor position data
        EnsemblDataCache geneTransCache = new EnsemblDataCache(outputDir, refGenomeVersion);
        geneTransCache.load(false);

        createTranscriptPreGenePositionData(
                geneTransCache.getChrGeneDataMap(), geneTransCache.getTranscriptDataMap(), PRE_GENE_PROMOTOR_DISTANCE, outputDir);

        LNX_LOGGER.info("Ensembl data cache complete");
    }

    private static void createTranscriptPreGenePositionData(
            final Map<String, List<EnsemblGeneData>> chrGeneDataMap, final Map<String, List<TranscriptData>> transcriptDataMap,
            int preGenePromotorDistance, final String outputDir)
    {
        // generate a cache file of the nearest upstream splice acceptor from another gene for each transcript
        try
        {
            String outputFile = outputDir;
            if (!outputFile.endsWith(File.separator))
                outputFile += File.separator;

            outputFile += ENSEMBL_TRANS_SPLICE_DATA_FILE;

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,TransId,TransName,TransStartPos,PreSpliceAcceptorPosition,Distance");
            writer.newLine();

            // for each gene, collect up any gene which overlaps it or is within the specified distance upstream from it
            for (Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
            {
                final String chromosome = entry.getKey();

                LNX_LOGGER.debug("calculating pre-gene positions for chromosome({})", chromosome);

                final List<EnsemblGeneData> geneList = entry.getValue();

                for (final EnsemblGeneData gene : geneList)
                {
                    List<String> proximateGenes = Lists.newArrayList();

                    int geneRangeStart;
                    int geneRangeEnd;

                    if (gene.Strand == POS_STRAND)
                    {
                        geneRangeStart = gene.GeneStart - preGenePromotorDistance;
                        geneRangeEnd = gene.GeneEnd;
                    }
                    else
                    {
                        geneRangeStart = gene.GeneStart;
                        geneRangeEnd = gene.GeneEnd + preGenePromotorDistance;
                    }

                    for (final EnsemblGeneData otherGene : geneList)
                    {
                        if (otherGene.Strand != gene.Strand)
                            continue;

                        //if (otherGene.GeneId.equals(gene.GeneId)) // skip same gene
                        //    continue;

                        if(positionsOverlap(geneRangeStart, geneRangeEnd, otherGene.GeneStart, otherGene.GeneEnd))
                        {
                            proximateGenes.add(otherGene.GeneId);
                        }
                    }

                    if(proximateGenes.isEmpty())
                        continue;

                    // now set the preceding splice acceptor position for each transcript in this gene
                    final List<TranscriptData> transDataList = transcriptDataMap.get(gene.GeneId);

                    if (transDataList == null || transDataList.isEmpty())
                        continue;

                    for(final TranscriptData transData : transDataList)
                    {
                        long transStartPos = gene.Strand == POS_STRAND ? transData.TransStart : transData.TransEnd;

                        long firstSpliceAcceptorPos = findFirstSpliceAcceptor(transData.TransId, transStartPos, gene.Strand, proximateGenes, transcriptDataMap);

                        if(firstSpliceAcceptorPos < 0)
                            continue;

                        long distance = gene.Strand == 1 ? transStartPos - firstSpliceAcceptorPos : firstSpliceAcceptorPos - transStartPos;

                        // cache value and continue
                        writer.write(String.format("%s,%d,%s,%d,%d,%d",
                                gene.GeneId, transData.TransId, transData.TransName, transStartPos, firstSpliceAcceptorPos, distance));

                        writer.newLine();
                    }
                }
            }

            LNX_LOGGER.info("pre-gene positions written to file: {}", outputFile);
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("error writing Ensembl trans splice data file: {}", e.toString());
        }
    }

    private static long findFirstSpliceAcceptor(
            int transId, long transStartPos, int strand,
            final List<String> proximateGenes, final Map<String, List<TranscriptData>> transDataMap)
    {
        long closestPosition = -1;

        for(final String geneId : proximateGenes)
        {
            final List<TranscriptData> transDataList = transDataMap.get(geneId);

            if(transDataList == null)
                continue;

            for(final TranscriptData transData : transDataList)
            {
                if(transId == transData.TransId)
                    continue;

                if(transData.exons().size() <= 1) // must have a splice acceptor
                    continue;

                for(ExonData exon : transData.exons())
                {
                    if(exon.Rank == 1)
                        continue;

                    if(strand == POS_STRAND && exon.Start < transStartPos)
                    {
                        if (closestPosition == -1 || exon.Start > closestPosition)
                            closestPosition = exon.Start;
                    }
                    else if(strand == NEG_STRAND && exon.End > transStartPos)
                    {
                        if (closestPosition == -1 || exon.End < closestPosition)
                            closestPosition = exon.End;
                    }
                }
            }
        }

        return closestPosition;
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(OUTPUT_DIR, true, "Directory to write Ensembl data files");
        options.addOption(LOG_DEBUG, false, "Log in verbose mode");
        EnsemblDAO.addCmdLineArgs(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }


}