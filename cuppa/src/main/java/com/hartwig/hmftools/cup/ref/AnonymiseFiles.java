package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_ALT_SJ_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_GENE_EXP_SAMPLE;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SAMPLE_POS_FREQ_COUNTS;
import static com.hartwig.hmftools.cup.CuppaRefFiles.REF_FILE_SNV_COUNTS;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_END;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_START;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class AnonymiseFiles
{
    private final AnonymiseConfig mConfig;

    private final SampleDataCache mSampleDataCache;
    private final Map<String,String> mSampleMapping;

    public AnonymiseFiles(final CommandLine cmd)
    {
        mConfig = new AnonymiseConfig(cmd);

        mSampleDataCache = new SampleDataCache();
        mSampleMapping = Maps.newHashMap();

        loadSampleData(cmd);
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadReferenceSampleData(mConfig.RefSampleDataFile);

        CUP_LOGGER.info("loaded {} reference samples", mSampleDataCache.RefSampleDataList.size());
    }

    public void run()
    {
        if(!mSampleDataCache.isValid())
        {
            CUP_LOGGER.error("invalid config");
            return;
        }

        CUP_LOGGER.info("writing anonymised ref data to {}", mConfig.OutputDir);

        for(int sampleIndex = 0; sampleIndex < mSampleDataCache.RefSampleDataList.size(); ++sampleIndex)
        {
            final SampleData sample = mSampleDataCache.RefSampleDataList.get(sampleIndex);
            getOrAddAnonSampleId(sample.Id);
        }

        rewriteCategoryFile(mConfig.RefSampleDataFile, REF_FILE_SAMPLE_DATA);

        rewriteMatrixFile(mConfig.RefSnvCountsFile, REF_FILE_SNV_COUNTS, Lists.newArrayList("BucketName"));

        rewriteMatrixFile(mConfig.RefSnvSamplePosFreqFile, REF_FILE_SAMPLE_POS_FREQ_COUNTS, Lists.newArrayList("BucketName"));

        rewriteMatrixFile(mConfig.RefGeneExpSampleFile, REF_FILE_GENE_EXP_SAMPLE, Lists.newArrayList(FLD_GENE_ID, FLD_GENE_NAME));

        rewriteMatrixFile(mConfig.RefAltSjSampleFile, REF_FILE_ALT_SJ_SAMPLE,
                Lists.newArrayList(FLD_GENE_ID, FLD_CHROMOSOME, FLD_POS_START, FLD_POS_END));

        CUP_LOGGER.info("CUP sampleIds anonymised");
    }

    private String getOrAddAnonSampleId(final String sampleId)
    {
        String anonSampleId = mSampleMapping.get(sampleId);

        if(anonSampleId == null)
        {
            anonSampleId = String.format("Sample_%04d", mSampleMapping.size());
            mSampleMapping.put(sampleId, anonSampleId);
        }

        return anonSampleId;
    }

    private void rewriteCategoryFile(final String refFilename, final String outputFilename)
    {
        if(!Files.exists(Paths.get(refFilename)))
            return;

        CUP_LOGGER.info("rewriting category file({})", refFilename);

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.OutputDir + outputFilename, false);

            final List<String> fileData = Files.readAllLines(new File(refFilename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            writer.write(header);
            writer.newLine();

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);
            final int sampleIdIndex = fieldsIndexMap.get("SampleId");

            for(final String line : fileData)
            {
                final String[] itemData = line.split(DATA_DELIM, -1);
                StringJoiner sj = new StringJoiner(DATA_DELIM);

                for(int i = 0; i < itemData.length; ++i)
                {
                    if(i == sampleIdIndex)
                    {
                        final String sampleId = itemData[i];
                        final String anonSampleId = getOrAddAnonSampleId(sampleId);
                        sj.add(anonSampleId);
                    }
                    else
                    {
                        sj.add(itemData[i]);
                    }
                }

                writer.write(sj.toString());
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read file({}) / write file({}): {}", refFilename, outputFilename, e.toString());
        }
    }

    private void rewriteMatrixFile(final String refFilename, final String outputFilename, final List<String> ignoreFields)
    {
        if(!Files.exists(Paths.get(refFilename)))
            return;

        CUP_LOGGER.info("rewriting matrix file({})", refFilename);

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.OutputDir + outputFilename, false);

            BufferedReader fileReader = createBufferedReader(refFilename);

            String header = fileReader.readLine();

            // convert any sampleIds in the column names
            final String[] itemData = header.split(DATA_DELIM, -1);
            StringJoiner sj = new StringJoiner(DATA_DELIM);

            for(int i = 0; i < itemData.length; ++i)
            {
                final String fieldName = itemData[i];

                if(ignoreFields.contains(fieldName))
                {
                    sj.add(fieldName);
                }
                else
                {
                    final String sampleId = fieldName;
                    final String anonSampleId = getOrAddAnonSampleId(sampleId);
                    sj.add(anonSampleId);
                }
            }

            writer.write(sj.toString());
            writer.newLine();

            String line = null;

            while((line = fileReader.readLine()) != null)
            {
                writer.write(line);
                writer.newLine();
            }

            closeBufferedWriter(writer);
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read file({}) / write file({}): {}", refFilename, outputFilename, e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        AnonymiseConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        AnonymiseFiles anonymiser = new AnonymiseFiles(cmd);
        anonymiser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
