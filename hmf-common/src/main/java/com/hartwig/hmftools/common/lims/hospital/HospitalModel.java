package com.hartwig.hmftools.common.lims.hospital;

import java.util.Map;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalModel {

    @NotNull
    abstract Map<String, HospitalAddress> hospitalAddressMap();

    @NotNull
    abstract Map<String, HospitalPersons> hospitalPersonsCPCT();

    @NotNull
    abstract Map<String, HospitalPersons> hospitalPersonsDRUP();

    @NotNull
    abstract Map<String, HospitalPersons> hospitalPersonsWIDE();

    @NotNull
    abstract Map<String, HospitalPersons> hospitalPersonsCOREDB();

    @NotNull
    abstract Map<String, String> sampleToHospitalMapping();

    @Nullable
    public HospitalContactData queryHospitalData(@NotNull String sampleId, @NotNull String coreRequesterName,
            @NotNull String coreRequesterEmail, @NotNull LimsCohortConfigData cohort) {
        String hospitalId = getHospitalIdForSample(sampleId, cohort);
        if (hospitalId == null) {
            return null;
        }

        if (cohort.cohortId().equals(Strings.EMPTY)) {
            return null;
        }

        HospitalAddress address = hospitalAddressMap().get(hospitalId);
        HospitalPersons persons = findPersonsForStudy(hospitalId, cohort);
        if (address == null || (persons == null && (cohort.hospitalPersonsStudy()))) {
            return null;
        }

        String requesterName = null;
        String requesterEmail = null;
        if (cohort.hospitalPersonsRequester()) {
            requesterName = coreRequesterName;
            requesterEmail = coreRequesterEmail;
        } else if (persons != null) {
            requesterName = persons.requesterName();
            requesterEmail = persons.requesterEmail();
        }

        return ImmutableHospitalContactData.builder()
                .hospitalPI(persons != null ? persons.hospitalPI() : Lims.NOT_AVAILABLE_STRING)
                .requesterName(requesterName != null ? requesterName : Lims.NOT_AVAILABLE_STRING)
                .requesterEmail(requesterEmail != null ? requesterEmail : Lims.NOT_AVAILABLE_STRING)
                .hospitalName(address.hospitalName())
                .hospitalAddress(address.hospitalZip() + " " + address.hospitalCity())
                .build();
    }

    @Nullable
    private HospitalPersons findPersonsForStudy(@NotNull String hospitalId, @NotNull LimsCohortConfigData cohort) {
        switch (cohort.cohortId()) {
            case "CPCT":
            case "CPCTpancreas":
                return hospitalPersonsCPCT().get(hospitalId);
            case "DRUP":
            case "DRUPstage3":
                return hospitalPersonsDRUP().get(hospitalId);
            case "WIDE":
                return hospitalPersonsWIDE().get(hospitalId);
            case "COREDB":
                return hospitalPersonsCOREDB().get(hospitalId);
            default:
                return null;
        }
    }

    @Nullable
    private String getHospitalIdForSample(@NotNull String sampleId, @NotNull LimsCohortConfigData cohort) {
        if (sampleToHospitalMapping().containsKey(sampleId)) {
            return sampleToHospitalMapping().get(sampleId);
        } else {

            if (cohort.hospitalId()) {
                // We assume all these projects follow a structure like CPCT##<hospitalId><identifier>
                return sampleId.substring(6, 8);
            } else {
                return null;
            }
        }
    }
}
