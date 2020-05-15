package com.hartwig.hmftools.common.lims.hospital;

import java.util.Map;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsStudy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalModel {

    private static final Logger LOGGER = LogManager.getLogger(HospitalModel.class);

    @NotNull
    abstract Map<String, HospitalAddress> hospitalAddress();

    @NotNull
    abstract Map<String, HospitalContact> hospitalContactCPCT();

    @NotNull
    abstract Map<String, HospitalContact> hospitalContactDRUP();

    @NotNull
    abstract Map<String, HospitalContact> hospitalContactWIDE();

    @NotNull
    abstract Map<String, String> sampleToHospitalMapping();

    @Nullable
    public HospitalData queryHospitalData(@NotNull String sampleId, @NotNull String coreRequesterName,
            @NotNull String coreRequesterEmail) {
        String hospitalId = getHospitalIdForSample(sampleId);
        if (hospitalId == null) {
            return null;
        }

        HospitalAddress address = hospitalAddress().get(hospitalId);
        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        HospitalContact contact = findContactDataForStudy(hospitalId, study);
        if (address == null || (contact == null && (study == LimsStudy.CPCT || study == LimsStudy.WIDE || study == LimsStudy.DRUP))) {
            return null;
        }

        String requesterName = null;
        String requesterEmail = null;
        if (study == LimsStudy.CORE) {
            requesterName = coreRequesterName;
            requesterEmail = coreRequesterEmail;
        } else if (contact != null) {
            requesterName = contact.requesterName();
            requesterEmail = contact.requesterEmail();
        }

        return ImmutableHospitalData.builder()
                .hospitalPI(contact != null ? contact.hospitalPI() : Lims.NOT_AVAILABLE_STRING)
                .requesterName(requesterName != null ? requesterName : Lims.NOT_AVAILABLE_STRING)
                .requesterEmail(requesterEmail != null ? requesterEmail : Lims.NOT_AVAILABLE_STRING)
                .hospitalName(address.hospitalName())
                .hospitalAddress(address.hospitalZip() + " " + address.hospitalCity())
                .build();
    }

    @Nullable
    private HospitalContact findContactDataForStudy(@NotNull String hospitalId, @NotNull LimsStudy study) {
        switch (study) {
            case CPCT:
                return hospitalContactCPCT().get(hospitalId);
            case DRUP:
                return hospitalContactDRUP().get(hospitalId);
            case WIDE:
                return hospitalContactWIDE().get(hospitalId);
            default:
                return null;
        }
    }

    @Nullable
    private String getHospitalIdForSample(@NotNull String sampleId) {
        if (sampleToHospitalMapping().containsKey(sampleId)) {
            return sampleToHospitalMapping().get(sampleId);
        } else {
            LimsStudy type = LimsStudy.fromSampleId(sampleId);

            if ((type == LimsStudy.DRUP || type == LimsStudy.CPCT || type == LimsStudy.WIDE || type == LimsStudy.CORE)
                    && sampleId.length() >= 12) {
                // We assume all these projects follow a structure like CPCT##<hospitalId><identifier>
                return sampleId.substring(6, 8);
            } else {
                return null;
            }
        }
    }
}
