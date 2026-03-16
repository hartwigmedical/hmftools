package com.hartwig.hmftools.pavereverse.batch;

import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.RPV_LOGGER;
import static com.hartwig.hmftools.pavereverse.batch.VariantsEncoder.columns;
import static com.hartwig.hmftools.pavereverse.parse.ProteinVariantParser.trimInitialPdot;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.ReversePave;

import org.apache.logging.log4j.util.BiConsumer;

public class BatchProcessor
{
    private final ReversePave mReversePave;

    public BatchProcessor(ReversePave reversePave)
    {
        mReversePave = reversePave;
    }

    public void process(String inputTsv, String outputTsv) throws IOException
    {
        BiConsumer<VariantRow, DelimFileWriter.Row> encoder = new VariantsEncoder();
        new File(outputTsv).createNewFile();
        try(DelimFileReader reader = new DelimFileReader(inputTsv);
                DelimFileWriter<VariantRow> writer = new DelimFileWriter<>(outputTsv, columns(), encoder))
        {
            for(DelimFileReader.Row row : reader)
            {
                String gene = row.get(0);
                String transcript = row.get(1);
                String proteinVariant = row.get(2);
                try
                {
                    String trimmedVariant = trimInitialPdot(proteinVariant);
                    BaseSequenceVariants variants =
                            mReversePave.calculateProteinVariantAllowingMultipleNonCanonicalTranscriptMatches(gene, trimmedVariant);
                    writer.writeRow(new VariantRow(gene, proteinVariant, variants));
                }
                catch(IllegalArgumentException e)
                {
                    RPV_LOGGER.error("Failed to compute variants for " + gene + ", " + proteinVariant + " and " + transcript, e);
                }
            }
        }
    }
}
