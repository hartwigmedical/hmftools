package com.hartwig.hmftools.protect.consent;

import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

public class ReportableVariantConsent {

    private final LimsGermlineReportingLevel reportingLevel;
    private final GermlineReportingModel germlineReportingModel;

    public ReportableVariantConsent(@NotNull LimsGermlineReportingLevel reportingLevel,
            @NotNull GermlineReportingModel germlineReportingModel) {
        this.reportingLevel = reportingLevel;
        this.germlineReportingModel = germlineReportingModel;
    }

    public boolean consentToReport(@NotNull ReportableVariant reportable) {
        return reportable.source().equals(ReportableVariantSource.PURPLE) || germlineReportingModel.reportableGermlineGenes()
                .contains(reportable.gene());
    }

    public boolean notifyClinicalGeneticist(@NotNull ReportableVariant reportable) {
        return reportable.source().equals(ReportableVariantSource.PURPLE) || germlineReportingModel.notifyAboutGene(reportingLevel,
                reportable.gene());
    }
}
