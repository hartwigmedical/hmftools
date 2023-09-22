package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

// routine for reading in the full ref-genome wildtype McfFlurry binding predictions (ran by Francisco)
// takes the per-transcript results and consolidates them into a single file
public class WildtypePeptidePredictions
{
    private final EnsemblDataCache mEnsemblDataCache;

    private final String mPredictionsDir;
    private final String mOutputFile;
    private BufferedWriter mWriter;

    private static final String PREDICTIONS_DIR = "predictions_dir";
    private static final String OUTPUT_FILE = "output_file";

    public WildtypePeptidePredictions(final ConfigBuilder configBuilder)
    {
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(false, false, false, true);

        mEnsemblDataCache.load(false);

        mPredictionsDir = configBuilder.getValue(PREDICTIONS_DIR);
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mWriter = initialiseWriter(parseOutputDir(configBuilder));
    }

    public void run()
    {
        NE_LOGGER.info("loading peptide predictions");

        for(Map.Entry<String,List<GeneData>> entry : mEnsemblDataCache.getChrGeneDataMap().entrySet())
        {
            for(GeneData geneData : entry.getValue())
            {
                TranscriptData transData = mEnsemblDataCache.getTranscriptData(geneData.GeneId, "");

                if(transData == null || transData.CodingStart == null)
                    continue;

                loadTranscriptFile(transData);
            }
        }

        closeBufferedWriter(mWriter);

        // NE_LOGGER.info("found {} random peptides from {} coding transcripts", totalPeptideCount, transCodingCount);
    }

    private void loadTranscriptFile(final TranscriptData transData)
    {
        /* ENST00000513876.predictions.tsv:

            peptide n_flank c_flank allele  mhcflurry_affinity      mhcflurry_affinity_percentile   pos_protein
            AVLYPLYFAR      LCATA   RECSP   HLA-A68:03      309.30710734621863      1.1514999999999993      103
         */

        final String filename = String.format("%s/%s.predictions.tsv", mPredictionsDir, transData.TransName);

        if(!Files.exists(Paths.get(filename)))
            return;

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(fileContents.get(0), "\t");
            fileContents.remove(0);
            int alleleIndex = fieldsIndexMap.get("allele");
            int peptideIndex = fieldsIndexMap.get("peptide");
            int affinityIndex = fieldsIndexMap.get("mhcflurry_affinity");
            int nfIndex = fieldsIndexMap.get("n_flank");
            int cfIndex = fieldsIndexMap.get("c_flank");

            for(String data : fileContents)
            {
                String[] items = data.split("\t");

                String peptide = items[peptideIndex];
                int peptideLength = peptide.length();

                if(peptideLength < 8 || peptideLength > 12)
                    continue;

                String allele = items[alleleIndex].replaceAll("HLA-","").replaceAll(":", "");
                double affinity = Double.parseDouble(items[affinityIndex]);

                writeData(transData.GeneId, transData.TransName, allele, peptide, affinity, items[nfIndex], items[cfIndex]);

            }

            NE_LOGGER.debug("gene({}) trans({}) loaded {} predictions", transData.GeneId, transData.TransName, fileContents.size());
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load transcript predictions file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            writer.write("GeneId,TransId,Allele,Peptide,Affinity,NFlank,CFlank");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    private void writeData(final String geneId, final String transId, final String allele, final String peptide,
            double affinity, final String nFlank, final String cFlank)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%.2f,%s,%s", geneId, transId, allele, peptide, affinity, nFlank, cFlank));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addEnsemblDir(configBuilder);
        configBuilder.addPath(OUTPUT_FILE, true, "Output filename");
        configBuilder.addPath(PREDICTIONS_DIR, true, "McfFlurry predictions directory");

        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        WildtypePeptidePredictions wtPredictions = new WildtypePeptidePredictions(configBuilder);
        wtPredictions.run();
    }
}
