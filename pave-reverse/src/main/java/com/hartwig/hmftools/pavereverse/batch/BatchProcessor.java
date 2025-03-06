package com.hartwig.hmftools.pavereverse.batch;

import static com.hartwig.hmftools.pavereverse.batch.BaseSequenceVariantsEncoder.columns;

import org.apache.logging.log4j.util.BiConsumer;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.ReversePave;

import org.jetbrains.annotations.NotNull;

public class BatchProcessor
{

    @NotNull
    private final ReversePave reversePave;

    public BatchProcessor(@NotNull final ReversePave reversePave)
    {
        this.reversePave = reversePave;
    }

    public void process(String inputTsv, String outputTsv)
    {
        BiConsumer<BaseSequenceVariants, DelimFileWriter.Row> encoder = new BaseSequenceVariantsEncoder();
        try(DelimFileReader reader = new DelimFileReader(inputTsv);
                DelimFileWriter<BaseSequenceVariants> writer = new DelimFileWriter<>(outputTsv, columns(), encoder))
        {
            for(DelimFileReader.Row row : reader)
            {
                String gene = row.get(0);
                String transcript = row.get(1);
                String proteinVariant = row.get(2);
                BaseSequenceVariants variants = reversePave.calculateVariant(gene, transcript, proteinVariant);
                writer.writeRow(variants);
            }
        }
    }
}
