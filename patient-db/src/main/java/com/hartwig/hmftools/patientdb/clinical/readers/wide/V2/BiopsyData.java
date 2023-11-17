package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;
import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class })
public abstract class BiopsyData // TODO FormRepeatKey, what is that ?
{
    @NotNull
    public abstract String combinedKey(); // basically is subjectKey and registrationDate appended, so each entry has a unique PK.

    @NotNull
    public abstract String subjectKey();

    @NotNull
    public abstract String formRepeatKey();

    public abstract Optional<LocalDate> sampleDate();

    public abstract Optional<String> sampleSite();

    public abstract Optional<String> sampleSiteDetails();

    public abstract Optional<String> sampleCollectMethod(); //TODO it says hmf-type = sample type? What?

    public abstract Optional<String> studyCode();

    public abstract Optional<String> otherTrial();

    public abstract Optional<String> otherTrialCode();

    public abstract Optional<LocalDate> otherTrialDate();

    public abstract Optional<String> diagnosis();

    public abstract Optional<String> BDMWDPNR(); // TODO, find out this field

    public abstract Optional<String> tNumber();

    public abstract Optional<Boolean> wasWgsSuccessful();

    public abstract Optional<String> reasonWgsWasNotSuccessful();

    public abstract Optional<String> sampleType();

    public abstract Optional<String> wgsReportPipelineVersion();

    public abstract Optional<LocalDate> hmfReportDate();

    public abstract Optional<String> wgsFindingsSummary();

    public abstract Optional<String> tumorTypeOncoTree();

    public abstract Optional<String> tumorTypeHmf();



    public static ImmutableBiopsyData.Builder builder() {
        return ImmutableBiopsyData.builder();
    }

}
