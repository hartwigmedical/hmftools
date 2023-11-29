package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.evidence.FragmentLengthData.ALT_COUNT;
import static com.hartwig.hmftools.sage.evidence.FragmentLengthData.REF_COUNT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.hartwig.hmftools.sage.SageConfig;

public class FragmentLengths
{
    private final SageConfig mConfig;
    private final BufferedWriter mWriter;

    public FragmentLengths(final SageConfig config)
    {
        mConfig = config;
        mWriter = initialiseWriter();
    }

    public void close() { closeBufferedWriter(mWriter); }

    private BufferedWriter initialiseWriter()
    {
        if(!mConfig.WriteFragmentLengths)
            return null;

        String outputVcf = mConfig.OutputFile;
        String fileName = outputVcf.replace(".vcf.gz", ".frag_lengths.tsv.gz");

        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId\tVariant\tLength\tRefCount\tAltCount");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise fragment length writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeVariantFragmentLength(
            final String variantInfo, final String sampleId, final FragmentLengthData fragmentLengthData)
    {
        if(mWriter == null)
            return;

        try
        {
            for(Map.Entry<Integer,int[]> entry : fragmentLengthData.lengthCounts().entrySet())
            {
                mWriter.write(format("%s\t%s\t%d\t%d\t%d",
                        sampleId, variantInfo, entry.getKey(), entry.getValue()[REF_COUNT], entry.getValue()[ALT_COUNT]));

                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise fragment length writer: {}", e.toString());
        }
    }
}
