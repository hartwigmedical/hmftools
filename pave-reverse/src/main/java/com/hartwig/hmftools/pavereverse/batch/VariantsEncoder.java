package com.hartwig.hmftools.pavereverse.batch;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.apache.logging.log4j.util.BiConsumer;

public class VariantsEncoder implements BiConsumer<VariantRow, DelimFileWriter.Row>
{
    static final String GENE = "GENE";
    static final String HGVS_PROTEIN = "HGVS_PROTEIN";
    static final String CHROMOSOME = "CHROMOSOME";
    static final String TRANSCRIPT = "TRANSCRIPT";
    static final String VARIANTS = "VARIANTS";

    public static Iterable<String> columns()
    {
        return List.of(GENE, HGVS_PROTEIN, CHROMOSOME, TRANSCRIPT, VARIANTS);
    }

    @Override
    public void accept(final VariantRow variantRow, final DelimFileWriter.Row row)
    {
        row.set(GENE, variantRow.mGene);
        row.set(HGVS_PROTEIN, variantRow.mProteinHGSV);
        row.set(CHROMOSOME, variantRow.mVariants.Chromosome);
        row.set(TRANSCRIPT, variantRow.mVariants.transcriptName());
        row.set(VARIANTS, encodeChanges(variantRow.mVariants.changes()));
    }

    static String encodeChanges(Set<BaseSequenceChange> changes)
    {
        return changes.stream().map(BaseSequenceChange::toCsv).collect(Collectors.joining(";"));
    }
}
