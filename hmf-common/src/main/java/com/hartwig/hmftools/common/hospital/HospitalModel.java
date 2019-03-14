package com.hartwig.hmftools.common.hospital;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsSampleType;

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
    protected abstract Map<String, HospitalData> hospitalPerId();

    @NotNull
    protected abstract Map<String, HospitalSampleMapping> sampleHospitalMapping();

    @Nullable
    public String addresseeStringForSample(@NotNull final String sample, @NotNull final String requesterName) {
        HospitalData hospital = findHospitalForSample(sample);

        if (hospital == null) {
            return null;
        }

        checkAddresseeFields(sample, hospital, requesterName);

        String requester = determineRequester(sample, hospital, requesterName);
        String hospitalAddress = hospital.addressName() + ", " + hospital.addressZip() + " " + hospital.addressCity();
        return requester.isEmpty() ? hospitalAddress : requester + ", " + hospitalAddress;
    }

    public int hospitalCount() {
        return hospitalPerId().values().size();
    }

    @Nullable
    @VisibleForTesting
    HospitalData hospitalPerId(@Nullable final String hospitalId) {
        return hospitalPerId().get(hospitalId);
    }

    @Nullable
    private HospitalData findHospitalForSample(@NotNull String sample) {
        final HospitalData hospital;
        if (sample.startsWith("CORE19") || sample.contains("CORE18")) {
            // These are the old core names, we need to manually map them.
            final HospitalSampleMapping hospitalSampleMapping = sampleHospitalMapping().get(sample);
            if (hospitalSampleMapping == null) {
                LOGGER.error("Cannot find sample hospital mapping for sample " + sample);
                return null;
            } else {
                hospital = findByHospital(hospitalSampleMapping.hospital());
                if (hospital == null) {
                    LOGGER.error("Cannot find hospital details for sample " + sample + " using " + hospitalSampleMapping.hospital());
                    return null;
                }
            }
        } else {
            final String hospitalId = extractHospitalIdFromSample(sample);
            if (hospitalId == null) {
                LOGGER.error("Could not extract hospital ID for sample " + sample);
                return null;
            }

            hospital = hospitalPerId().get(hospitalId);
            if (hospital == null) {
                LOGGER.error("Hospital model does not contain id " + hospitalId);
                return null;
            }
        }

        return hospital;
    }

    @Nullable
    private HospitalData findByHospital(@NotNull String hospital) {
        for (HospitalData hospitalData : hospitalPerId().values()) {
            if (hospitalData.hospital().equals(hospital)) {
                return hospitalData;
            }
        }
        return null;
    }

    @Nullable
    private static String extractHospitalIdFromSample(@NotNull final String sample) {
        LimsSampleType type = LimsSampleType.fromSampleId(sample);

        if (type == LimsSampleType.DRUP || type == LimsSampleType.CPCT || type == LimsSampleType.WIDE
                || type == LimsSampleType.CORE && sample.length() >= 12) {
            // We assume all these projects follow a structure like CPCT##<hospital><identifier>
            return sample.substring(6, 8);
        }

        LOGGER.error("Sample parameter: " + sample + " is not in CPCT/DRUP/WIDE/CORE format");
        return null;
    }

    private static void checkAddresseeFields(@NotNull final String sample, @NotNull final HospitalData hospital,
            @NotNull final String requesterName) {
        final List<String> missingFields = Lists.newArrayList();
        if (determineRequester(sample, hospital, requesterName).isEmpty()) {
            missingFields.add("requester");
        }
        if (hospital.addressName().isEmpty()) {
            missingFields.add("name");
        }
        if (hospital.addressZip().isEmpty()) {
            missingFields.add("zip");
        }
        if (hospital.addressCity().isEmpty()) {
            missingFields.add("city");
        }
        if (!missingFields.isEmpty()) {
            LOGGER.warn("Some address fields (" + Strings.join(missingFields, ',') + ") are missing.");
        }
    }

    @NotNull
    @VisibleForTesting
    static String determineRequester(@NotNull final String sample, @NotNull final HospitalData hospital,
            @NotNull final String requesterName) {
        LimsSampleType type = LimsSampleType.fromSampleId(sample);

        if (type == LimsSampleType.CPCT) {
            return hospital.cpctPI();
        } else if (type == LimsSampleType.DRUP) {
            final String drupPi = hospital.drupPI();
            if (drupPi.trim().equals("*")) {
                return hospital.cpctPI();
            }
            return hospital.drupPI();
        } else if (type == LimsSampleType.WIDE) {
            return requesterName;
        }
        return Strings.EMPTY;
    }
}
