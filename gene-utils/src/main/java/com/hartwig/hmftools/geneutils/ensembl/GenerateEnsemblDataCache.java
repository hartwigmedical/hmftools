package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.ENSEMBL_TRANS_SPLICE_DATA_FILE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.DEFAULT_PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

import org.jetbrains.annotations.NotNull;

public class GenerateEnsemblDataCache
{
    public static final String HGNC_GENE_DATA_FILE = "hgnc_gene_data_file";
    public static final String REF_ENSEMBL_DIR = "ref_ensembl_dir";

    public static void writeEnsemblDataFiles(final ConfigBuilder configBuilder)
    {
        String outputDir = parseOutputDir(configBuilder);

        RefGenomeVersion refGenomeVersion = RefGenomeVersion.from(configBuilder);

        GU_LOGGER.info("writing Ensembl version({}) data files to {}", refGenomeVersion, outputDir);

        if(EnsemblDAO.hasDatabaseConfig(configBuilder))
        {
            EnsemblDAO ensemblDAO = new EnsemblDAO(configBuilder);

            if(!ensemblDAO.isValid())
            {
                GU_LOGGER.info("invalid Ensembl DAO");
                System.exit(1);
            }

            ensemblDAO.writeDataCacheFiles(outputDir);
        }

        GU_LOGGER.debug("reloading transcript data to generate splice acceptor positions");

        // create the transcript splice acceptor position data
        EnsemblDataCache geneTransCache = new EnsemblDataCache(outputDir, refGenomeVersion);
        geneTransCache.load(false);

        createTranscriptPreGenePositionData(
                geneTransCache.getChrGeneDataMap(), geneTransCache.getTranscriptDataMap(), DEFAULT_PRE_GENE_PROMOTOR_DISTANCE, outputDir);

        GU_LOGGER.info("Ensembl data cache complete");
    }

    private static void createTranscriptPreGenePositionData(
            final Map<String, List<GeneData>> chrGeneDataMap, final Map<String, List<TranscriptData>> transcriptDataMap,
            int preGenePromotorDistance, final String outputDir)
    {
        // generate a cache file of the nearest upstream splice acceptor from another gene for each transcript
        try
        {
            String outputFile = outputDir + ENSEMBL_TRANS_SPLICE_DATA_FILE;

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,TransId,TransName,TransStartPos,PreSpliceAcceptorPosition,Distance");
            writer.newLine();

            // for each gene, collect up any gene which overlaps it or is within the specified distance upstream from it
            for(Map.Entry<String, List<GeneData>> entry : chrGeneDataMap.entrySet())
            {
                final String chromosome = entry.getKey();

                GU_LOGGER.debug("calculating pre-gene positions for chromosome({})", chromosome);

                final List<GeneData> geneList = entry.getValue();

                for(final GeneData gene : geneList)
                {
                    List<String> proximateGenes = Lists.newArrayList();

                    int geneRangeStart;
                    int geneRangeEnd;

                    if(gene.Strand == POS_STRAND)
                    {
                        geneRangeStart = gene.GeneStart - preGenePromotorDistance;
                        geneRangeEnd = gene.GeneEnd;
                    }
                    else
                    {
                        geneRangeStart = gene.GeneStart;
                        geneRangeEnd = gene.GeneEnd + preGenePromotorDistance;
                    }

                    for(final GeneData otherGene : geneList)
                    {
                        if(otherGene.Strand != gene.Strand)
                            continue;

                        //if(otherGene.GeneId.equals(gene.GeneId)) // skip same gene
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

                    if(transDataList == null || transDataList.isEmpty())
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

            GU_LOGGER.info("pre-gene positions written to file: {}", outputFile);
            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            GU_LOGGER.error("error writing Ensembl trans splice data file: {}", e.toString());
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
                        if(closestPosition == -1 || exon.Start > closestPosition)
                            closestPosition = exon.Start;
                    }
                    else if(strand == NEG_STRAND && exon.End > transStartPos)
                    {
                        if(closestPosition == -1 || exon.End < closestPosition)
                            closestPosition = exon.End;
                    }
                }
            }
        }

        return closestPosition;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        configBuilder.addPath(HGNC_GENE_DATA_FILE, false, "HGNC gene data file");
        configBuilder.addPath(REF_ENSEMBL_DIR, false, "Ensembl data dir for canonical trans to use");
        EnsemblDAO.addCmdLineArgs(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
        addOutputDir(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        writeEnsemblDataFiles(configBuilder);
    }
}