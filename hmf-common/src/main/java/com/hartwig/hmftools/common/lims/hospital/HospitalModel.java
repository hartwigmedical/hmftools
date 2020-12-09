package com.hartwig.hmftools.common.lims.hospital;

import java.util.Map;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsCohort;

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
    abstract Map<String, String> sampleToHospitalMapping();

    @Nullable
    public HospitalContactData queryHospitalData(@NotNull String sampleId, @NotNull String coreRequesterName,
            @NotNull String coreRequesterEmail, @NotNull LimsCohort cohort) {
        String hospitalId = getHospitalIdForSample(sampleId, cohort);
        if (hospitalId == null) {
            return null;
        }

        HospitalAddress address = hospitalAddressMap().get(hospitalId);
        HospitalPersons persons = findPersonsForStudy(hospitalId, cohort);
        if (address == null || (persons == null && (cohort == LimsCohort.CPCT || cohort == LimsCohort.WIDE || cohort == LimsCohort.DRUP))) {
            return null;
        }

        String requesterName = null;
        String requesterEmail = null;
        if (cohort == LimsCohort.CORE) {
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
    private HospitalPersons findPersonsForStudy(@NotNull String hospitalId, @NotNull LimsCohort cohort) {
        switch (cohort) {
            case CPCT:
                return hospitalPersonsCPCT().get(hospitalId);
            case DRUP:
                return hospitalPersonsDRUP().get(hospitalId);
            case WIDE:
                return hospitalPersonsWIDE().get(hospitalId);
            default:
                return null;
        }
    }

    @Nullable
    private String getHospitalIdForSample(@NotNull String sampleId, @NotNull LimsCohort cohort) {
        if (sampleToHospitalMapping().containsKey(sampleId)) {
            return sampleToHospitalMapping().get(sampleId);
        } else {

            if ((cohort == LimsCohort.DRUP || cohort == LimsCohort.CPCT || cohort == LimsCohort.WIDE || cohort == LimsCohort.CORE)
                    && sampleId.length() >= 12) {
                // We assume all these projects follow a structure like CPCT##<hospitalId><identifier>
                return sampleId.substring(6, 8);
            } else {
                return null;
            }
        }
    }
}
