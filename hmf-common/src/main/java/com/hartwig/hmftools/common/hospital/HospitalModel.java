package com.hartwig.hmftools.common.hospital;

import java.util.List;
import java.util.Map;

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
    abstract Map<String, HospitalData> hospitalDataCPCT();

    @NotNull
    abstract Map<String, HospitalData> hospitalDataDRUP();

    @NotNull
    abstract Map<String, HospitalData> hospitalDataWIDE();

    @NotNull
    abstract Map<String, HospitalAdress> hospitalAdress();

    @Nullable
    public String extractHospitalPI(@NotNull String sampleId) {
        HospitalData hospitalData;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            hospitalData = hospitalDataCPCT().get(hospitalID);
            return hospitalData.hospitalPI();
        } else if (study == LimsStudy.DRUP) {
            hospitalData = hospitalDataDRUP().get(hospitalID);
            return hospitalData.hospitalPI();
        } else if (study == LimsStudy.WIDE) {
            hospitalData = hospitalDataWIDE().get(hospitalID);
            return hospitalData.hospitalPI();
        } else if (study == LimsStudy.CORE) {
            return null;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    public String extractRequestName(@NotNull String sampleId, @NotNull String requestNameCore) {
        HospitalData hospitalData;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            return null;
        } else if (study == LimsStudy.DRUP) {
            return null;
        } else if (study == LimsStudy.WIDE) {
            hospitalData = hospitalDataWIDE().get(hospitalID);
            return hospitalData.requestName();
        } else if (study == LimsStudy.CORE) {
            return requestNameCore;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    public String extractRequestEmail(@NotNull String sampleId, @NotNull String requestEmailCore) {
        HospitalData hospitalData;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        LimsStudy study = LimsStudy.fromSampleId(sampleId);
        if (study == LimsStudy.CPCT) {
            return null;
        } else if (study == LimsStudy.DRUP) {
            return null;
        } else if (study == LimsStudy.WIDE) {
            hospitalData = hospitalDataWIDE().get(hospitalID);
            return hospitalData.requestEmail();
        } else if (study == LimsStudy.CORE) {
            return requestEmailCore;
        } else if (study == LimsStudy.NON_CANCER_STUDY) {
            return null;
        } else {
            return null;
        }
    }

    @Nullable
    private HospitalAdress findByHospitalCore(@NotNull String hospital) {
        for (HospitalAdress hospitalAdress : hospitalAdress().values()) {
            if (hospitalAdress.hospitalName().equals(hospital)) {
                return hospitalAdress;
            }
        }
        return null;
    }

    @NotNull
    private HospitalAdress extractHospitalAdressDataOldCoreNames(@NotNull String sampleId) {
        HospitalAdress hospitalAdress;
        // These are the old core names, we need to manually map them.
        HospitalSampleMapping hospitalSampleMapping = sampleHospitalMapping().get(sampleId);
        if (hospitalSampleMapping == null) {
            LOGGER.error("Cannot find sample hospital mapping for sample '{}'.", sampleId);
            return null;
        } else {
            hospitalAdress = findByHospitalCore(hospitalSampleMapping.internalHospitalName());
            if (hospitalAdress == null) {
                LOGGER.error("Cannot find hospital details for sample '{}' using '{}'.",
                        sampleId,
                        hospitalSampleMapping.internalHospitalName());
                return null;
            }
        }
        return hospitalAdress;
    }

    @Nullable
    private HospitalAdress checkInputHospital(@Nullable String hospitalID, @NotNull String sampleId) {
        HospitalAdress hospitalAdress;
        if (hospitalID == null) {
            LOGGER.warn("Could not find hospital for sample '{}'", sampleId);
            return null;
        }

        hospitalAdress = hospitalAdress().get(hospitalID);
        if (hospitalAdress == null) {
            LOGGER.warn("Hospital model does not contain id '{}'", hospitalID);
            return null;
        }
        return hospitalAdress;
    }

    @Nullable
    public String extractHospitalAdress(@NotNull String sampleId) {
        HospitalAdress hospitalAdress;
        String hospitalID = extractHospitalIdFromSample(sampleId);
        LimsStudy study = LimsStudy.fromSampleId(sampleId);

        if (sampleId.startsWith("CORE19") || sampleId.contains("CORE18")) {
            hospitalAdress = extractHospitalAdressDataOldCoreNames(sampleId);

            checkAddresseeFields(hospitalAdress);
            return hospitalAdress.hospitalName() + ", " + hospitalAdress.hospitalZip() + " "
                    + hospitalAdress.hospitalCity();

        } else {
            hospitalAdress = checkInputHospital(hospitalID, sampleId);
            checkAddresseeFields(hospitalAdress);
            if (study == LimsStudy.CORE) {
                return hospitalAdress.hospitalName() + ", " + hospitalAdress.hospitalZip() + " "
                        + hospitalAdress.hospitalCity();
            } else {
                return extractHospitalPI(sampleId) + ", " + hospitalAdress.hospitalName() + ", " + hospitalAdress.hospitalZip() + " "
                        + hospitalAdress.hospitalCity();
            }
        }
    }

    @Nullable
    public String extractHospitalName(@NotNull String sampleId) {
        HospitalAdress hospitalAdress;
        String hospitalID = extractHospitalIdFromSample(sampleId);

        if (sampleId.startsWith("CORE19") || sampleId.contains("CORE18")) {
            hospitalAdress = extractHospitalAdressDataOldCoreNames(sampleId);

            checkAddresseeFields(hospitalAdress);
            return hospitalAdress.hospitalName();

        } else {
            hospitalAdress = checkInputHospital(hospitalID, sampleId);

            checkAddresseeFields(hospitalAdress);

            return hospitalAdress.hospitalName();
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

    private static void checkAddresseeFields(@NotNull HospitalAdress hospital) {
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

    @NotNull
    public HospitalQuery generateHospitalQuery(@NotNull String sampleId, @NotNull String requestNameCore, @NotNull String requestEmailCore) {
        return ImmutableHospitalQuery.builder()
                .hospitalPI(extractHospitalPI(sampleId))
                .analyseRequestName(extractRequestName(sampleId, requestNameCore))
                .analyseRequestEmail(extractRequestEmail(sampleId, requestEmailCore))
                .fullHospitalString(extractHospitalAdress(sampleId))
                .hospital(extractHospitalName(sampleId))
                .build();
    }

}
