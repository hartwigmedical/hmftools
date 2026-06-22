package com.hartwig.hmftools.orange.report.pdfdata;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.jetbrains.annotations.Nullable;

public class GermlineFindingsData
{
    public final boolean hasPurpleFail;
    public final boolean hasReferenceId;
    public final boolean hasRna;
    @Nullable
    public final List<PurpleDriver> germlineDrivers;
    @Nullable
    public final List<PurpleVariant> germlineVariants;
    @Nullable
    public final List<PurpleGainDeletion> germlineGainsDels;
    @Nullable
    public final List<LinxBreakend> germlineBreakends;
    public final Set<PurpleGermlineAberration> germlineAberrations;
    @Nullable
    public final Set<PeachGenotype> peach;

    public GermlineFindingsData(
            final boolean hasPurpleFail,
            final boolean hasReferenceId,
            final boolean hasRna,
            @Nullable final List<PurpleDriver> germlineDrivers,
            @Nullable final List<PurpleVariant> germlineVariants,
            @Nullable final List<PurpleGainDeletion> germlineGainsDels,
            @Nullable final List<LinxBreakend> germlineBreakends,
            final Set<PurpleGermlineAberration> germlineAberrations,
            @Nullable final Set<PeachGenotype> peach)
    {
        this.hasPurpleFail = hasPurpleFail;
        this.hasReferenceId = hasReferenceId;
        this.hasRna = hasRna;
        this.germlineDrivers = germlineDrivers;
        this.germlineVariants = germlineVariants;
        this.germlineGainsDels = germlineGainsDels;
        this.germlineBreakends = germlineBreakends;
        this.germlineAberrations = germlineAberrations;
        this.peach = peach;
    }
}
