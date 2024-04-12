package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.sage.VariantFragmentLength.VARIANT_FRAG_LENGTHS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.common.sage.VariantFragmentLength;
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
        String fileName = outputVcf.replace(".vcf.gz", VARIANT_FRAG_LENGTHS_FILE_ID);

        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(VariantFragmentLength.header());
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
            final String variantInfo, final String sampleId, final FragmentLengthCounts fragmentLengthData)
    {
        if(mWriter == null)
            return;

        try
        {
            VariantFragmentLength.writeVariantFragmentLength(mWriter, sampleId, variantInfo, fragmentLengthData);
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write fragment length writer: {}", e.toString());
        }
    }
}
