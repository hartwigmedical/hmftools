package com.hartwig.hmftools.finding;

import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import jakarta.validation.constraints.NotNull;

public interface EventFactory
{
    @NotNull
    String variantEvent(@NotNull PurpleVariant variant);

    @NotNull
    String gainDeletionEvent(@NotNull PurpleGeneCopyNumber geneCopyNumber);

    @NotNull
    String gainDeletionEvent(@NotNull PurpleGainDeletion purpleGainDeletion);

    @NotNull
    String gainDeletionEvent(@NotNull PurpleLossOfHeterozygosity loh);

    @NotNull
    String disruptionEvent(@NotNull LinxBreakend breakend);

    @NotNull
    String fusionEvent(@NotNull LinxFusion fusion);

    @NotNull
    String virusEvent(@NotNull VirusInterpreterEntry virus);

    @NotNull
    String immunologyEvent(@NotNull LilacAllele lilacAllele);
}
