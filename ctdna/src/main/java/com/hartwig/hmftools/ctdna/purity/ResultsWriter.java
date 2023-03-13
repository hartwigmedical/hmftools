package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.variant.VariantContextDecorator;

public class ResultsWriter
{
    private final InterpretConfig mConfig;
    private final BufferedWriter mVariantWriter;

    public ResultsWriter(final InterpretConfig config)
    {
        mConfig = config;
        mVariantWriter = initialiseWriter();
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mConfig.OutputDir + "patient_interpret_results";

            if(mConfig.OutputId != null)
                fileName += "." + mConfig.OutputId;

            fileName += ".csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("PatientId,SampleId,Chromosome,Position,Ref,Alt,Tier,Type,RepeatCount,Mappability,SubclonalPerc,AD,DP,QualPerAD");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeSampleVariant(
            final String patientId, final String sampleId, final VariantContextDecorator variant,
            double subclonalLikelihood, int alleleCount, int depth, double qualPerAlleleCount)
    {
        try
        {
            mVariantWriter.write(String.format("%s,%s,%s,%d,%s,%s",
                    patientId, sampleId, variant.chromosome(), variant.position(), variant.ref(), variant.alt()));

            mVariantWriter.write(String.format(",%s,%s,%d,%.3f,%.3f",
                    variant.tier(), variant.type(), variant.repeatCount(), variant.mappability(), subclonalLikelihood));

            mVariantWriter.write(String.format(",%d,%d,%.1f", alleleCount, depth, qualPerAlleleCount));

            mVariantWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }
    public void close()
    {
        closeBufferedWriter(mVariantWriter);
    }

}
