package com.hartwig.hmftools.patientdb.clinical.readers.drup;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.clinical.matchers.BiopsyMatcher;
import com.hartwig.hmftools.patientdb.clinical.matchers.MatchResult;
import com.hartwig.hmftools.patientdb.clinical.readers.EcrfPatientReader;

import org.jetbrains.annotations.NotNull;

public class DrupPatientReader implements EcrfPatientReader {

    @NotNull
    private final BaselineReader baselineReader;
    @NotNull
    private final BiopsyReader biopsyReader;

    public DrupPatientReader(@NotNull PrimaryTumorCurator primaryTumorCurator, @NotNull BiopsySiteCurator biopsySiteCurator) {
        this.baselineReader = new BaselineReader(primaryTumorCurator);
        this.biopsyReader = new BiopsyReader(biopsySiteCurator);
    }

    @NotNull
    @Override
    public Patient read(@NotNull EcrfPatient ecrfPatient, @NotNull List<SampleData> sequencedSamples, @NotNull String consentConfigTsv) {
        BaselineData baselineData = baselineReader.read(ecrfPatient);
        PreTreatmentData noPreTreatmentData = ImmutablePreTreatmentData.builder().formStatus(FormStatus.undefined()).build();
        List<BiopsyData> clinicalBiopsies = biopsyReader.read(ecrfPatient, baselineData.curatedPrimaryTumor());

        MatchResult<BiopsyData> matchedBiopsies =
                BiopsyMatcher.matchBiopsiesToTumorSamples(ecrfPatient.patientId(), sequencedSamples, clinicalBiopsies);

        return new Patient(ecrfPatient.patientId(),
                baselineData,
                noPreTreatmentData,
                sequencedSamples,
                matchedBiopsies.values(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                matchedBiopsies.findings());
    }
}
