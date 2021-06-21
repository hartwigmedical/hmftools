package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.sage.variant.SageVariant;

public class VariantFile
{
    public final BufferedWriter mWriter;

    public VariantFile(final String sampleId, final String outputDir)
    {
        mWriter = initialise(sampleId, outputDir);
    }

    private BufferedWriter initialise(final String sampleId, final String outputDir)
    {
        String fileName = outputDir + sampleId + ".sage.variants.csv";
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("Chromosome,Position,Ref,Alt");
            writer.write(",Pass,MaxReadDepth,Tier,TotalQuality");
            writer.write(",LocalPhaseSet,IsRealigned,LocalRealignSet");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return null;
        }
    }

    public void writeToFile(final SageVariant variant)
    {
        if(mWriter == null)
            return;

        try
        {
            mWriter.write(String.format("%s,%d,%s,%s", variant.chromosome(), variant.position(), variant.ref(), variant.alt()));

            mWriter.write(String.format(",%s,%d,%s,%d",
                    variant.isPassing(), variant.candidate().maxReadDepth(), variant.candidate().tier(), variant.totalQuality()));

            mWriter.write(String.format(",%d,%s,%d",
                    variant.localPhaseSet(), variant.isRealigned(), variant.localRealignSet()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write variant data: {}", e.toString());
        }
    }

    public void close() { closeBufferedWriter(mWriter); }

}
