package com.hartwig.hmftools.pavereverse.batch;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.logging.log4j.util.BiConsumer;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;

public class BaseSequenceVariantsEncoder implements BiConsumer<BaseSequenceVariants, DelimFileWriter.Row>
{
    static final String CHROMOSOME = "CHROMOSOME";
    static final String TRANSCRIPT = "TRANSCRIPT";
    static final String VARIANTS = "VARIANTS";

    static Iterable<String> columns()
    {
        return List.of(CHROMOSOME, TRANSCRIPT, VARIANTS);
    }

    @Override
    public void accept(final BaseSequenceVariants variants, final DelimFileWriter.Row row)
    {
        row.set(CHROMOSOME, variants.mChromosome);
        row.set(TRANSCRIPT, variants.transcriptName());
        row.set(VARIANTS, encodeChanges(variants.changes()));
    }

    static String encodeChanges(Set<BaseSequenceChange> changes)
    {
        return changes.stream().map(BaseSequenceChange::toCsv).collect(Collectors.joining(";"));
    }
}
