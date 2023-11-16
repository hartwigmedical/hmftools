package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BiopsyData // TODO FormRepeatKey, what is that ?
{
    @NotNull
    public abstract String combinedKey(); // basically is subjectKey and registrationDate appended, so each entry has a unique PK.

    @NotNull
    public abstract String subjectKey();

    @NotNull
    public abstract LocalDate registrationDate();

    public abstract LocalDate sampleDate();

    public abstract String sampleSite();

    public abstract String sampleSiteDetails();

    public abstract String sampleCollectMethod(); //TODO it says hmf-type = sample type? What?

    public abstract String studyCode();

    public abstract String otherTrial();

    public abstract String otherTrialCode();
    public abstract LocalDate otherTrialDate();

    public abstract String diagnosis();

    public abstract String BDMWDPNR(); // TODO, find out this field

    public abstract String tNumber();

    public abstract boolean wasWgsSuccessful();

    public abstract String reasonWgsWasNotSuccessful();

    public abstract String sampleType();

    public abstract String wgsReportPipelineVersion();

    public abstract LocalDate hmfReportDate();

    public abstract String wgsFindingsSummary();

    public abstract String tumorTypeOncoTree();

    public abstract String tumorTypeHmf();

}
