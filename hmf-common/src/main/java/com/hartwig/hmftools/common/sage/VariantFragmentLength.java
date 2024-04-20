package com.hartwig.hmftools.common.sage;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sage.FragmentLengthCounts.ALT_COUNT;
import static com.hartwig.hmftools.common.sage.FragmentLengthCounts.REF_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class VariantFragmentLength
{
    public final String SampleId;
    public final String VariantInfo;
    public final int Length;
    public final int RefCount;
    public final int AltCount;

    public static final String VARIANT_FRAG_LENGTHS_FILE_ID = ".frag_lengths.tsv.gz";

    public VariantFragmentLength(final String sampleId, final String variantInfo, final int length, final int refCount, final int altCount)
    {
        SampleId = sampleId;
        VariantInfo = variantInfo;
        Length = length;
        RefCount = refCount;
        AltCount = altCount;
    }

    private static final String FLD_VARIANT = "Variant";
    public static final String FLD_REF_COUNT = "RefCount";
    public static final String FLD_ALT_COUNT = "AltCount";
    public static final String FLD_LENGTH = "Length";

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(FLD_SAMPLE_ID);
        sj.add(FLD_VARIANT);
        sj.add(FLD_LENGTH);
        sj.add(FLD_REF_COUNT);
        sj.add(FLD_ALT_COUNT);
        return sj.toString();
    }

    public static void writeVariantFragmentLength(
            final BufferedWriter writer, final String sampleId, final String variantInfo, final FragmentLengthCounts fragmentLengthData)
            throws IOException
    {
        if(writer == null)
            return;

        for(Map.Entry<Integer,int[]> entry : fragmentLengthData.lengthCounts().entrySet())
        {
            writer.write(format("%s\t%s\t%d\t%d\t%d",
                    sampleId, variantInfo, entry.getKey(), entry.getValue()[REF_COUNT], entry.getValue()[ALT_COUNT]));

            writer.newLine();
        }
    }

    public static List<VariantFragmentLength> read(final String filename)
    {
        List<VariantFragmentLength> fragLengthData = Lists.newArrayList();

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);
            String header = fileReader.readLine();

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int sampleIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
            int varIndex = fieldsIndexMap.get(FLD_VARIANT);
            int refIndex = fieldsIndexMap.get(FLD_REF_COUNT);
            int altIndex = fieldsIndexMap.get(FLD_ALT_COUNT);
            int lengthIndex = fieldsIndexMap.get(FLD_LENGTH);

            String line;
            while((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(TSV_DELIM, -1);

                fragLengthData.add(new VariantFragmentLength(
                        values[sampleIndex], values[varIndex],Integer.parseInt(values[lengthIndex]),
                        Integer.parseInt(values[refIndex]), Integer.parseInt(values[altIndex])));
            }

            return fragLengthData;
        }
        catch(Exception e)
        {
            return null;
        }
    }
}
