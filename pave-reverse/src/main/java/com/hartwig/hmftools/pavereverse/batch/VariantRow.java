package com.hartwig.hmftools.pavereverse.batch;

import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;

public class VariantRow
{
    final String mGene;
    final String mProteinHGSV;
    final BaseSequenceVariants mVariants;

    public VariantRow(final String mGene, final String mProteinHGSV, final BaseSequenceVariants variants)
    {
        this.mGene = mGene;
        this.mProteinHGSV = mProteinHGSV;
        this.mVariants = variants;
    }
}
