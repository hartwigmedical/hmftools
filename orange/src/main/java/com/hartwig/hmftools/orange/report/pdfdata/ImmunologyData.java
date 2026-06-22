package com.hartwig.hmftools.orange.report.pdfdata;

import com.hartwig.hmftools.datamodel.hla.LilacRecord;

import org.jetbrains.annotations.Nullable;

public class ImmunologyData
{
    public final boolean hasPurpleFail;
    @Nullable
    public final LilacRecord lilac;
    public final boolean hasRna;
    public final boolean isWholeGenome;

    public ImmunologyData(
            final boolean hasPurpleFail,
            @Nullable final LilacRecord lilac,
            final boolean hasRna,
            final boolean isWholeGenome)
    {
        this.hasPurpleFail = hasPurpleFail;
        this.lilac = lilac;
        this.hasRna = hasRna;
        this.isWholeGenome = isWholeGenome;
    }
}
