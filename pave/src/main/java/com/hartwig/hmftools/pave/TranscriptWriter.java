package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.pave.impact.CodingContext;
import com.hartwig.hmftools.pave.impact.ProteinContext;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

public class TranscriptWriter
{
    private final BufferedWriter mWriter;
    private final String mSampleId;

    public TranscriptWriter(final PaveConfig config)
    {
        mWriter = config.WriteTranscriptFile ? initialiseWriter(config) : null;
        mSampleId = config.SampleId;
    }

    private BufferedWriter initialiseWriter(final PaveConfig config)
    {
        try
        {
            String fileSuffix = ".pave.transcript" + TSV_EXTENSION;
            String transFileName = config.OutputDir + config.SampleId + fileSuffix;
            BufferedWriter writer = createBufferedWriter(transFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(VariantData.tsvHeader());
            sj.add(VariantData.extraDataHeader());

            sj.add("GeneId\tGeneName");
            sj.add(VariantTransImpact.tsvHeader());
            sj.add(CodingContext.tsvHeader());
            sj.add(ProteinContext.tsvHeader());

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeVariantData(final VariantData variant, final String geneName)
    {
        if(mWriter == null)
            return;

        List<VariantTransImpact> geneImpacts = variant.getGeneImpacts(geneName);

        if(geneImpacts == null)
            return;

        try
        {
            for(VariantTransImpact impact : geneImpacts)
            {
                if(impact.TransData == null)
                    continue;

                mWriter.write(String.format("%s\t%s", variant.tsvData(), variant.extraDataTsv(mSampleId)));

                mWriter.write(String.format("\t%s\t%s\t%s\t%s\t%s",
                        impact.TransData.GeneId, geneName, impact.toTsv(), impact.codingContext().toTsv(),
                        impact.proteinContext() != null ? impact.proteinContext().toTsv() : ProteinContext.empty()));

                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to write variant CSV file: {}", e.toString());
            return;
        }
    }

    public void close() { closeBufferedWriter(mWriter); }

}
