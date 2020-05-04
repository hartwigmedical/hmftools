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
    private static final String NA_STRING = "N/A";

    @NotNull
    abstract Map<String, HospitalData> hospitalPerId();

    @NotNull
    abstract Map<String, HospitalSampleMapping> sampleHospitalMapping();

    @NotNull
    abstract Map<String, HospitalCore> hospitalCoreMap();

    public int hospitalCount() {
        return hospitalPerId().values().size();
    }

    @NotNull
    public HospitalQuery queryHospitalDataForSample(@NotNull String sampleId) {
        LimsStudy study = LimsStudy.fromSampleId(sampleId);

        if (study == LimsStudy.CORE) {
            HospitalCore hospitalCore = findHospitalForSampleCore(sampleId);
            return ImmutableHospitalQuery.builder()
                    .hospitalPA(NA_STRING)
                    .analyseRequestName(NA_STRING)
                    .analyseRequestEmail(NA_STRING)
                    .hospitalId(NA_STRING)
                    .hospitalName(NA_STRING)
                    .hospitalAdres(NA_STRING)
                    .build();
        } else {
            HospitalData hospital = findHospitalForSample(sampleId);
            return ImmutableHospitalQuery.builder()
                    .hospitalPA(NA_STRING)
                    .analyseRequestName(NA_STRING)
                    .analyseRequestEmail(NA_STRING)
                    .hospitalId(NA_STRING)
                    .hospitalName(NA_STRING)
                    .hospitalAdres(NA_STRING)
                    .build();
        }
    }

    @NotNull
    private static String fullAddresseeStringCore(@NotNull HospitalCore hospital) {
        return hospital.externalHospitalName() + ", " + hospital.addressZip() + " " + hospital.addressCity();
    }

    @NotNull
    private static String fullAddresseeString(@NotNull String sample, @NotNull HospitalData hospital) {
        checkAddresseeFields(sample, hospital);

        String hospitalPI = determinePIName(sample, hospital);
        String hospitalAddress = hospital.externalHospitalName() + ", " + hospital.addressZip() + " " + hospital.addressCity();
        return hospitalPI.isEmpty() ? hospitalAddress : hospitalPI + ", " + hospitalAddress;
    }

    @Nullable
    @VisibleForTesting
    HospitalData hospitalPerId(@Nullable String hospitalId) {
        return hospitalPerId().get(hospitalId);
    }

    @Nullable
    private HospitalCore findHospitalForSampleCore(@NotNull String sample) {
        HospitalCore hospital;
        if (sample.startsWith("CORE19") || sample.contains("CORE18")) {
            // These are the old core names, we need to manually map them.
            HospitalSampleMapping hospitalSampleMapping = sampleHospitalMapping().get(sample);
            if (hospitalSampleMapping == null) {
                LOGGER.error("Cannot find sample hospital mapping for sample '{}'.", sample);
                return null;
            } else {
                hospital = findByHospitalCore(hospitalSampleMapping.internalHospitalName());
                if (hospital == null) {
                    LOGGER.error("Cannot find hospital details for sample '{}' using '{}'.",
                            sample,
                            hospitalSampleMapping.internalHospitalName());
                    return null;
                }
            }
        } else {
            String hospitalId = extractHospitalIdFromSample(sample);
            if (hospitalId == null) {
                LOGGER.warn("Could not find hospital for sample '{}'", sample);
                return null;
            }

            hospital = hospitalCoreMap().get(hospitalId);
            if (hospital == null) {
                LOGGER.warn("Hospital model does not contain id '{}'", hospitalId);
                return null;
            }
        }

        return hospital;
    }

    @Nullable
    private HospitalCore findByHospitalCore(@NotNull String hospital) {
        for (HospitalCore hospitalCore : hospitalCoreMap().values()) {
            if (hospitalCore.internalHospitalName().equals(hospital)) {
                return hospitalCore;
            }
        }
        return null;
    }

    @Nullable
    private HospitalData findHospitalForSample(@NotNull String sample) {
        HospitalData hospital;
        String hospitalId = extractHospitalIdFromSample(sample);
        if (hospitalId == null) {
            LOGGER.warn("Could not find hospital for sample '{}'", sample);
            return null;
        }

        hospital = hospitalPerId().get(hospitalId);
        if (hospital == null) {
            LOGGER.warn("Hospital model does not contain id '{}'", hospitalId);
            return null;
        }
        return hospital;
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

    private static void checkAddresseeFields(@NotNull String sample, @NotNull HospitalData hospital) {
        List<String> missingFields = Lists.newArrayList();

        LimsStudy type = LimsStudy.fromSampleId(sample);
        if (type != LimsStudy.CORE && determinePIName(sample, hospital).isEmpty()) {
            missingFields.add("requester");
        }
        if (hospital.externalHospitalName().isEmpty()) {
            missingFields.add("name");
        }
        if (hospital.addressZip().isEmpty()) {
            missingFields.add("zip");
        }
        if (hospital.addressCity().isEmpty()) {
            missingFields.add("city");
        }

        if (!missingFields.isEmpty()) {
            LOGGER.warn("Some address fields ({}) are missing.", Strings.join(missingFields, ','));
        }
    }

    @NotNull
    @VisibleForTesting
    static String determinePIName(@NotNull String sample, @NotNull HospitalData hospital) {
        LimsStudy type = LimsStudy.fromSampleId(sample);

        if (type == LimsStudy.CPCT) {
            return hospital.cpctPI();
        } else if (type == LimsStudy.DRUP) {
            String drupPi = hospital.drupPI();
            if (drupPi.trim().equals("*")) {
                return hospital.cpctPI();
            }
            return hospital.drupPI();
        } else if (type == LimsStudy.WIDE) {
            return hospital.widePI();
        } else if (type == LimsStudy.CORE || type == LimsStudy.NON_CANCER_STUDY) {
            return Strings.EMPTY;
        }

        return NA_STRING;
    }

    @NotNull
    @VisibleForTesting
    static String determinePIEmail(@NotNull String sample, @NotNull HospitalData hospital) {
        LimsStudy type = LimsStudy.fromSampleId(sample);

        if (type == LimsStudy.CPCT) {
            return extractPIEmailFromRecipientList(hospital.cpctRecipients());
        } else if (type == LimsStudy.DRUP) {
            if (hospital.drupPI().trim().equals("*")) {
                return extractPIEmailFromRecipientList(hospital.cpctRecipients());
            }
            return extractPIEmailFromRecipientList(hospital.drupRecipients());
        } else if (type == LimsStudy.WIDE) {
            return extractPIEmailFromRecipientList(hospital.wideRecipients());
        }

        return Strings.EMPTY;
    }

    @NotNull
    private static String extractPIEmailFromRecipientList(@NotNull String recipients) {
        // We assume the first email in a list separated with ";" is the PI email.
        return recipients.split(";")[0];
    }
}
