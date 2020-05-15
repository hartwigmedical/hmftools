package com.hartwig.hmftools.common.lims.hospital;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsStudy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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

    @NotNull
    public HospitalData queryHospitalData(@NotNull String sampleId, @NotNull String requesterNameCore,
            @NotNull String requesterEmailCore) {
        return ImmutableHospitalData.builder()
                .hospitalPI(extractHospitalPI(sampleId))
                .requesterName(extractRequestName(sampleId, requesterNameCore))
                .requesterEmail(extractRequestEmail(sampleId, requesterEmailCore))
                .hospitalName(extractHospitalName(sampleId))
                .hospitalAddress(extractHospitalAddress(sampleId))
                .build();
    }

    private boolean checkIdExistInInputFile(HospitalContact hospitalContact, @Nullable String hospitalID) {
        if (hospitalContact == null) {
            return false;
        } else {
            return true;
        }
    }

    @Nullable
    @VisibleForTesting
    String extractHospitalPI(@NotNull String sampleId) {
        HospitalContact hospitalContact;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            hospitalContact = hospitalContactCPCT().get(hospitalID);
            boolean existsInFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existsInFile ? hospitalContact.hospitalPI() : null;
        } else if (study == LimsStudy.DRUP) {
            hospitalContact = hospitalContactDRUP().get(hospitalID);
            boolean existsInFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existsInFile ? hospitalContact.hospitalPI() : null;
        } else if (study == LimsStudy.WIDE) {
            hospitalContact = hospitalContactWIDE().get(hospitalID);
            boolean existsInFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existsInFile ? hospitalContact.hospitalPI() : null;
        } else if (study == LimsStudy.CORE) {
            return null;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    String extractRequestName(@NotNull String sampleId, @NotNull String requesterNameCore) {
        HospitalContact hospitalContact;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            return null;
        } else if (study == LimsStudy.DRUP) {
            return null;
        } else if (study == LimsStudy.WIDE) {
            hospitalContact = hospitalContactWIDE().get(hospitalID);
            boolean existsInFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existsInFile ? hospitalContact.requesterName() : null;
        } else if (study == LimsStudy.CORE) {
            return requesterNameCore;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    String extractRequestEmail(@NotNull String sampleId, @NotNull String requesterEmailCore) {
        HospitalContact hospitalContact;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            return null;
        } else if (study == LimsStudy.DRUP) {
            return null;
        } else if (study == LimsStudy.WIDE) {
            hospitalContact = hospitalContactWIDE().get(hospitalID);
            boolean existsInFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existsInFile ? hospitalContact.requesterEmail() : null;
        } else if (study == LimsStudy.CORE) {
            return requesterEmailCore;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    private HospitalAddress findByHospitalCore(@NotNull String hospitalId) {
        for (HospitalAddress hospitalAddress : hospitalAddress().values()) {
            if (hospitalAddress.hospitalId().equals(hospitalId)) {
                return hospitalAddress;
            }
        }
        return null;
    }

    @Nullable
    private HospitalAddress findHospitalAddressModel(@NotNull String sampleId) {
        String hospitalSampleMapping = sampleToHospitalMapping().get(sampleId);
        HospitalAddress hospitalAddress;

        if (hospitalSampleMapping != null) {
            hospitalAddress = findByHospitalCore(hospitalSampleMapping);
            if (hospitalAddress == null) {
                return null;
            }
        } else {
            String hospitalID = extractHospitalIdFromSample(sampleId);
            if (hospitalID == null) {
                LOGGER.warn("Could not find hospital for sample '{}'", sampleId);
                return null;
            }
            hospitalAddress = hospitalAddress().get(hospitalID);
            if (hospitalAddress == null) {
                return null;
            }

        }
        checkAddresseeFields(hospitalAddress);
        return hospitalAddress;
    }

    @Nullable
    @VisibleForTesting
    String extractHospitalAddress(@NotNull String sampleId) {
        HospitalAddress hospitalAddress = findHospitalAddressModel(sampleId);

        if (hospitalAddress == null) {
            return null;
        } else {
            return hospitalAddress.hospitalName() + ", " + hospitalAddress.hospitalZip() + " " + hospitalAddress.hospitalCity();
        }
    }

    @Nullable
    @VisibleForTesting
    String extractHospitalName(@NotNull String sampleId) {
        HospitalAddress hospitalAddress = findHospitalAddressModel(sampleId);

        if (hospitalAddress == null) {
            return null;
        } else {
            return hospitalAddress.hospitalName();
        }

    }

    @Nullable
    private static String extractHospitalIdFromSample(@NotNull String sample) {
        LimsStudy type = LimsStudy.fromSampleId(sample);

        if (type == LimsStudy.DRUP || type == LimsStudy.CPCT || type == LimsStudy.WIDE || type == LimsStudy.CORE && sample.length() >= 12) {
            // We assume all these projects follow a structure like CPCT##<hospital><identifier>
            return sample.substring(6, 8);
        }

        return null;
    }

    private static void checkAddresseeFields(@NotNull HospitalAddress hospital) {
        List<String> missingFields = Lists.newArrayList();

        if (hospital.hospitalName().isEmpty()) {
            missingFields.add("name");
        }
        if (hospital.hospitalZip().isEmpty()) {
            missingFields.add("zip");
        }
        if (hospital.hospitalCity().isEmpty()) {
            missingFields.add("city");
        }

        if (!missingFields.isEmpty()) {
            LOGGER.warn("Some address fields ({}) are missing.", Strings.join(missingFields, ','));
        }
    }
}
