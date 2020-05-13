package com.hartwig.hmftools.common.hospital;

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
    abstract Map<String, HospitalSampleMapping> sampleHospitalMapping();

    @NotNull
    abstract Map<String, HospitalContact> hospitalContactCPCT();

    @NotNull
    abstract Map<String, HospitalContact> hospitalContactDRUP();

    @NotNull
    abstract Map<String, HospitalContact> hospitalContactWIDE();

    @NotNull
    abstract Map<String, HospitalAddress> hospitalAddress();

    @NotNull
    public HospitalQuery generateHospitalQuery(@NotNull String sampleId, @NotNull String requestNameCore,
            @NotNull String requestEmailCore) {
        return ImmutableHospitalQuery.builder()
                .hospitalPI(extractHospitalPI(sampleId))
                .requestName(extractRequestName(sampleId, requestNameCore))
                .requestEmail(extractRequestEmail(sampleId, requestEmailCore))
                .hospitalAddress(extractHospitalAddress(sampleId))
                .hospitalName(extractHospitalName(sampleId))
                .build();
    }

    private boolean checkIdExistInInputFile(HospitalContact hospitalContact, @Nullable String hospitalID) {
        if (hospitalContact == null) {
            LOGGER.warn("Hospital model does not contain id '{}'", hospitalID);
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
            boolean existinFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existinFile ? hospitalContact.hospitalPI() : null;
        } else if (study == LimsStudy.DRUP) {
            hospitalContact = hospitalContactDRUP().get(hospitalID);
            boolean existinFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existinFile ? hospitalContact.hospitalPI() : null;
        } else if (study == LimsStudy.WIDE) {
            hospitalContact = hospitalContactWIDE().get(hospitalID);
            boolean existinFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existinFile ? hospitalContact.hospitalPI() : null;
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
    String extractRequestName(@NotNull String sampleId, @NotNull String requestNameCore) {
        HospitalContact hospitalContact;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            return null;
        } else if (study == LimsStudy.DRUP) {
            return null;
        } else if (study == LimsStudy.WIDE) {
            hospitalContact = hospitalContactWIDE().get(hospitalID);
            boolean existinFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existinFile ? hospitalContact.requestName() : null;
        } else if (study == LimsStudy.CORE) {
            return requestNameCore;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    @VisibleForTesting
    String extractRequestEmail(@NotNull String sampleId, @NotNull String requestEmailCore) {
        HospitalContact hospitalContact;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            return null;
        } else if (study == LimsStudy.DRUP) {
            return null;
        } else if (study == LimsStudy.WIDE) {
            hospitalContact = hospitalContactWIDE().get(hospitalID);
            boolean existinFile = checkIdExistInInputFile(hospitalContact, hospitalID);
            return existinFile ? hospitalContact.requestEmail() : null;
        } else if (study == LimsStudy.CORE) {
            return requestEmailCore;
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
    @VisibleForTesting
    String extractHospitalAddress(@NotNull String sampleId) {
        HospitalAddress hospitalAddress;
        String hospitalID = extractHospitalIdFromSample(sampleId);
        LimsStudy study = LimsStudy.fromSampleId(sampleId);

        if (sampleId.startsWith("CORE19") || sampleId.contains("CORE18")) {
            // These are the old core names, we need to manually map them.
            HospitalSampleMapping hospitalSampleMapping = sampleHospitalMapping().get(sampleId);
            if (hospitalSampleMapping == null) {
                LOGGER.error("Cannot find sample hospital mapping for sample '{}'.", sampleId);
                return null;
            } else {
                hospitalAddress = findByHospitalCore(hospitalSampleMapping.hospitalId());
                if (hospitalAddress == null) {
                    LOGGER.error("Cannot find hospital details for sample '{}' using '{}'.", sampleId, hospitalSampleMapping.hospitalId());
                    return null;
                }
            }

            checkAddresseeFields(hospitalAddress);
            return hospitalAddress.hospitalName() + ", " + hospitalAddress.hospitalZip() + " " + hospitalAddress.hospitalCity();

        } else {
            if (hospitalID == null) {
                LOGGER.warn("Could not find hospital for sample '{}'", sampleId);
                return null;
            }

            hospitalAddress = hospitalAddress().get(hospitalID);
            if (hospitalAddress == null) {
                LOGGER.warn("Hospital model does not contain id '{}'", hospitalID);
                return null;
            }
            checkAddresseeFields(hospitalAddress);
            return hospitalAddress.hospitalName() + ", " + hospitalAddress.hospitalZip() + " " + hospitalAddress.hospitalCity();
        }
    }

    @Nullable
    @VisibleForTesting
    public String extractHospitalName(@NotNull String sampleId) {
        HospitalAddress hospitalAddress;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        if (sampleId.startsWith("CORE19") || sampleId.contains("CORE18")) {
            // These are the old core names, we need to manually map them.
            HospitalSampleMapping hospitalSampleMapping = sampleHospitalMapping().get(sampleId);
            if (hospitalSampleMapping == null) {
                LOGGER.error("Cannot find sample hospital mapping for sample '{}'.", sampleId);
                return null;
            } else {
                hospitalAddress = findByHospitalCore(hospitalSampleMapping.hospitalId());
                if (hospitalAddress == null) {
                    LOGGER.error("Cannot find hospital details for sample '{}' using '{}'.", sampleId, hospitalSampleMapping.hospitalId());
                    return null;
                }
            }

            checkAddresseeFields(hospitalAddress);
            return hospitalAddress.hospitalName();

        } else {
            if (hospitalID == null) {
                LOGGER.warn("Could not find hospital for sample '{}'", sampleId);
                return null;
            }

            hospitalAddress = hospitalAddress().get(hospitalID);
            if (hospitalAddress == null) {
                LOGGER.warn("Hospital model does not contain id '{}'", hospitalID);
                return null;
            }

            checkAddresseeFields(hospitalAddress);

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
