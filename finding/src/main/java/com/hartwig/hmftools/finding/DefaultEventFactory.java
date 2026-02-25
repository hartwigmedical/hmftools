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

// TODO: Still need to fill this in. What conventions to follow for these events could be location specific
//  (ACTIN and OncoAct already use different conventions), so hence it's "pluggable"
class DefaultEventFactory implements EventFactory
{
    @NotNull
    public String variantEvent(@NotNull PurpleVariant variant)
    {
        return variant.gene() + " " + impact(variant);
    }

    @NotNull
    private String impact(@NotNull PurpleVariant variant)
    {
        return "";
    }

    @NotNull
    public String gainDeletionEvent(@NotNull PurpleGeneCopyNumber geneCopyNumber)
    {
        return geneCopyNumber.gene() + " disruption ";
    }

    @NotNull
    public String gainDeletionEvent(@NotNull PurpleGainDeletion purpleGainDeletion)
    {
        return purpleGainDeletion.gene() + " disruption ";
    }

    @NotNull
    public String gainDeletionEvent(@NotNull PurpleLossOfHeterozygosity loh)
    {
        return loh.gene() + " disruption ";
    }

    @NotNull
    public String disruptionEvent(@NotNull LinxBreakend breakend)
    {
        return breakend.gene() + " disruption ";
    }

    @NotNull
    public String fusionEvent(@NotNull LinxFusion fusion)
    {
        return "";
    }

    @NotNull
    public String virusEvent(@NotNull VirusInterpreterEntry virus)
    {
        return "";
    }

    @NotNull
    public String immunologyEvent(@NotNull LilacAllele lilacAllele)
    {
        return "HLA-" + lilacAllele.allele();
    }
}
