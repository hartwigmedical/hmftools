package com.hartwig.hmftools.dnds;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.dnds.DndsCommon.DN_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;

public class SampleMutationalLoad
{
    public final double Purity; // currently unused
    public final int SnvBiallelic;
    public final int SnvNonBiallelic;
    public final int IndelBiallelic;
    public final int IndelNonBiallelic;

    public SampleMutationalLoad(
            final double purity, final int snvBiallelic, final int snvNonBiallelic, final int indelBiallelic, final int indelNonBiallelic)
    {
        Purity = purity;
        SnvBiallelic = snvBiallelic;
        SnvNonBiallelic = snvNonBiallelic;
        IndelBiallelic = indelBiallelic;
        IndelNonBiallelic = indelNonBiallelic;
    }

    public int snvTotal() { return SnvBiallelic + SnvNonBiallelic; }
    public int indelTotal() { return IndelBiallelic + IndelNonBiallelic; }

    private static final String SAMPLE_ID = "SampleId";
    private static final String PURITY = "Purity";
    private static final String SNV_BI = "SnvBiallelic";
    private static final String SNV_NON_BI = "SnvNonBiallelic";
    private static final String INDEL_BI = "IndelBiallelic";
    private static final String INDEL_NON_BI = "IndelNonBiallelic";

    public static final String COHORT_SAMPLE_MUT_LOAD = "dnds_cohort_mut_load.tsv";

    public static String cohortSampleMutationalLoadFilename(final String sourceDir)
    {
        return sourceDir + COHORT_SAMPLE_MUT_LOAD;
    }

    public static BufferedWriter initialiseWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(SAMPLE_ID).add(PURITY).add(SNV_BI).add(SNV_NON_BI).add(INDEL_BI).add(INDEL_NON_BI);
            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            DN_LOGGER.error("failed to create cohort sample mut-load writer: {}", e.toString());
            return null;
        }
    }

    public synchronized  static void writeSampleMutationalLoad(
            final BufferedWriter writer, final String sampleId, final SampleMutationalLoad mutLoad)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(sampleId);
            sj.add(format("%.4f", mutLoad.Purity));
            sj.add(String.valueOf(mutLoad.SnvBiallelic));
            sj.add(String.valueOf(mutLoad.SnvNonBiallelic));
            sj.add(String.valueOf(mutLoad.IndelBiallelic));
            sj.add(String.valueOf(mutLoad.IndelNonBiallelic));
            writer.write(sj.toString());
            writer.newLine();
        }
        catch (IOException e)
        {
            DN_LOGGER.error("failed to write cohort sample mut-load writer: {}", e.toString());
        }
    }

    public static Map<String,SampleMutationalLoad> loadCohortSampleMutationalLoads(final String filename, boolean expectExists)
    {
        Map<String,SampleMutationalLoad> sampleMutationalLoadMap = Maps.newHashMap();

        if(!Files.exists(Paths.get(filename)))
        {
            return expectExists ? null : sampleMutationalLoadMap;
        }

        try
        {
            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), TSV_DELIM);
            lines.remove(0);

            int sampleIndex = fieldsIndexMap.get(SAMPLE_ID);
            int purityIndex = fieldsIndexMap.get(PURITY);
            int snvIndex = fieldsIndexMap.get(SNV_BI);
            int snvNbIndex = fieldsIndexMap.get(SNV_NON_BI);
            int indelIndex = fieldsIndexMap.get(INDEL_BI);
            int indelNbIndex = fieldsIndexMap.get(INDEL_NON_BI);

            for(String line : lines)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                String sampleId = values[sampleIndex];

                sampleMutationalLoadMap.put(sampleId, new SampleMutationalLoad(
                        Double.parseDouble(values[purityIndex]), Integer.parseInt(values[snvIndex]), Integer.parseInt(values[snvNbIndex]),
                        Integer.parseInt(values[indelIndex]), Integer.parseInt(values[indelNbIndex])));
            }
        }
        catch(IOException e)
        {
            DN_LOGGER.error("failed to read cohort sample mut-load writer: {}", e.toString());
            return null;
        }

        return sampleMutationalLoadMap;
    }
}
