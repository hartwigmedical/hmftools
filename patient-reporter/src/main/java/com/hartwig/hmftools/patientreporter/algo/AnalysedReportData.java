package com.hartwig.hmftools.patientreporter.algo;

import com.hartwig.hmftools.patientreporter.ReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryModel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class AnalysedReportData implements ReportData {

    @NotNull
    public abstract GermlineReportingModel germlineReportingModel();

    @NotNull
    public abstract SummaryModel summaryModel();

    @NotNull
    public abstract VirusDbModel virusDbModel();

    @NotNull
    public abstract VirusSummaryModel virusSummaryModel();

    @NotNull
    public abstract VirusBlacklistModel virusBlackListModel();
}
