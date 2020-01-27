package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.formOutputPath;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RnaExpression
{
    private static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    private static final String GENE_ID_FILE = "gene_id_file";
    private static final String ALL_TRANSCRIPTS = "all_transcripts";
    private static final String SAMPLE = "sample";

    private static final Logger LOGGER = LogManager.getLogger(RnaExpression.class);

    private final String mSampledId;
    private final RnaBamReader mRnaBamReader;
    private final SvGeneTranscriptCollection mGeneTransCache;
    private final List<String> mRestrictedGeneIds;
    private final String mOutputDir;
    private final boolean mAllTranscripts;
    private final List<GeneReadData> mGeneReadDatalist;
    private BufferedWriter mWriter;

    public RnaExpression(final CommandLine cmd)
    {
        mRnaBamReader = new RnaBamReader(cmd);

        mSampledId = cmd.getOptionValue(SAMPLE);

        mGeneTransCache = new SvGeneTranscriptCollection();
        mGeneTransCache.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));

        mRestrictedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE));
            mGeneTransCache.setRestrictedGeneIdList(mRestrictedGeneIds);
        }

        mAllTranscripts = cmd.hasOption(ALL_TRANSCRIPTS);

        mGeneTransCache.setRequiredData(true, false, false, !cmd.hasOption(ALL_TRANSCRIPTS));
        mGeneTransCache.loadEnsemblData(false);

        mOutputDir = formOutputPath(cmd.getOptionValue(DATA_OUTPUT_DIR));

        mGeneReadDatalist = Lists.newArrayList();
        mWriter = null;
    }

    public void runAnalysis()
    {
        // measure read counts of exonic regions for all specific genes
        for(Map.Entry<String,List<EnsemblGeneData>> entry : mGeneTransCache.getChrGeneDataMap().entrySet())
        {
            final List<EnsemblGeneData> genesDataList = entry.getValue();
            genesDataList.forEach(x -> processGene(x));
        }

        closeBufferedWriter(mWriter);
    }

    private void processGene(final EnsemblGeneData geneData)
    {
        GeneReadData geneReadData = new GeneReadData(geneData);

        final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);
        if(transDataList.isEmpty())
        {
            LOGGER.warn("no transcripts found for gene({}:{})", geneData.GeneId, geneData.GeneName);
            return;
        }

        // form a genomic region for each unique exon amongst the transcripts
        long minTransPos = -1;
        long maxTransPos = 0;

        for(final TranscriptData transData : transDataList)
        {
            for(ExonData exon : transData.exons())
            {
                if(geneReadData.hasRegionData(exon.ExonStart, exon.ExonEnd))
                    continue;

                GenomeRegion region = GenomeRegions.create(geneData.Chromosome, exon.ExonStart, exon.ExonEnd);
                RegionReadData regionReadData = new RegionReadData(
                        region, String.format("trans(%s) exon(%d)", transData.TransName, exon.ExonRank));
                geneReadData.addRegionReadData(regionReadData);
            }

            maxTransPos = max(transData.TransEnd, maxTransPos);

            if(minTransPos < 0 || transData.TransStart < minTransPos)
                minTransPos = transData.TransStart;
        }

        GenomeRegion geneRegion = GenomeRegions.create(geneData.Chromosome, minTransPos, maxTransPos);
        mRnaBamReader.readBamCounts(geneRegion);
        mRnaBamReader.analyseReads(geneReadData);

        // report evidence for each gene transcript
        for(final TranscriptData transData : transDataList)
        {
            writeResults(geneReadData, transData);
        }


        mGeneReadDatalist.add(geneReadData);
    }

    private void writeResults(final GeneReadData geneReadData, final TranscriptData transData)
    {
        if(mOutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mOutputDir + "RNA_GENE_EXPRESSION.csv";

                mWriter = createBufferedWriter(outputFileName, false);
                mWriter.write("SampleId,GeneId,GeneName,TransId,ExonCount,ExonsMatched,LinksMatched,AvgDepth");
                mWriter.newLine();
            }

            int exonsFound = 0;
            int linksFound = 0;
            double readDepthTotal = 0;

            final List<ExonData> exons = transData.exons();

            for(int i = 0; i < exons.size(); ++i)
            {
                ExonData exon = exons.get(i);

                final RegionReadData exonReadData = geneReadData.findRegionData(exon.ExonStart, exon.ExonEnd);
                if(exonReadData == null)
                    continue;

                ++exonsFound;

                ExonData prevExon = i > 0 ? exons.get(i - 1) : null;
                ExonData nextExon = i < exons.size() - 1 ? exons.get(i + 1) : null;

                RegionReadData prevRegion = prevExon != null ? geneReadData.findRegionData(prevExon.ExonStart, prevExon.ExonEnd) : null;
                RegionReadData nextRegion = nextExon != null ? geneReadData.findRegionData(nextExon.ExonStart, nextExon.ExonEnd) : null;

                boolean linked = (prevExon == null || exonReadData.getLinkedRegions().containsKey(prevRegion))
                        && (nextExon == null || exonReadData.getLinkedRegions().containsKey(nextRegion));

                if(linked)
                    ++linksFound;

                readDepthTotal += exonReadData.averageDepth();
            }

            double avgReadDepth = exonsFound > 0 ? readDepthTotal / exonsFound : 0;

            mWriter.write(String.format("%s,%s,%s,%s",
                    mSampledId, geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, transData.TransName));

            mWriter.write(String.format(",%d,%d,%d,%.0f",
                    exons.size(), exonsFound, linksFound, avgReadDepth));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write gene expression file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if(!validConfig(cmd))
        {
            LOGGER.error("missing config options, exiting");
            return;
        }

        RnaExpression rnaExpression = new RnaExpression(cmd);
        rnaExpression.runAnalysis();

        LOGGER.info("RNA expression analysis complete");
    }

    public static boolean validConfig(final CommandLine cmd)
    {
        if(!cmd.hasOption(SAMPLE) || !cmd.hasOption(GENE_TRANSCRIPTS_DIR))
            return false;

        return RnaBamReader.validConfig(cmd);
    }

    private void loadGeneIdsFile(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("GeneId"))
            {
                // check for header row
                fileContents.remove(0);
            }

            mRestrictedGeneIds.addAll(fileContents.stream().map(x -> x.split(",")[0]).collect(Collectors.toList()));

            LOGGER.info("file({}) loaded {} genes", filename, mRestrictedGeneIds.size());
        }
        catch (IOException e)
        {
            LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache");

        options.addOption(ALL_TRANSCRIPTS, false, "Check all transcripts, not just canonical");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        RnaBamReader.addCommandLineOptions(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
