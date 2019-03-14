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
    protected abstract Map<String, HospitalData> hospitalPerHospital();

    @NotNull
    protected abstract Map<String, HospitalSampleMapping> hospitalPerIdManual();

    @Nullable
    public String addresseeStringForSample(@NotNull final String contactNames, @NotNull final String sample) {
        String address;
        if (sample.startsWith("CORE19") || sample.contains("CORE18")) { // These are the old core names
            final HospitalSampleMapping hospitalSampleMapping = hospitalPerIdManual().get(sample);
            final HospitalData hospital = hospitalPerHospital().get(hospitalSampleMapping.addressName());
            if (hospital == null) {
                LOGGER.error("HospitalModelFactory model cannot find hospital details for project " + hospitalSampleMapping);
                return null;
            }
            address = hospital.addressName() + ", " + hospital.addressZip() + " " + hospital.addressCity();

        } else {
            final String HospitalId = getHospitalIdFromSample(sample);
            if (HospitalId == null) {
                return null;
            }
            final HospitalData hospital = hospitalPerId(HospitalId);
            if (hospital == null) {
                LOGGER.error("Hospital model does not contain id " + HospitalId);
                return null;
            }
            checkAddresseeFields(sample, hospital, contactNames);
            address = getPI(sample, hospital, contactNames) + ", " + hospital.addressName() + ", " + hospital.addressZip() + " "
                    + hospital.addressCity();
        }
        return address;
    }

    @Nullable
    @VisibleForTesting
    HospitalData hospitalPerId(@Nullable final String hospitalId) {
        return hospitalPerId().get(hospitalId);
    }

    @Nullable
    private static String getHospitalIdFromSample(@NotNull final String sample) {
        final String ucSample = sample.toUpperCase();
        LimsSampleType type = LimsSampleType.fromSampleId(ucSample);

        if (type == LimsSampleType.DRUP || type == LimsSampleType.CPCT || type == LimsSampleType.WIDE
                || type == LimsSampleType.CORE && sample.length() >= 12) {
            return sample.substring(6, 8);
        }

        LOGGER.warn("Sample parameter: " + sample + " is not in CPCT/DRUP/WIDE/CORE format");
        return null;
    }

    private static void checkAddresseeFields(@NotNull final String sample, @NotNull final HospitalData hospital,
            @NotNull final String contactNames) {
        final List<String> missingFields = Lists.newArrayList();
        if (getPI(sample, hospital, contactNames).isEmpty()) {
            missingFields.add("PI");
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
    static String getPI(@NotNull final String sample, @NotNull final HospitalData hospital, @NotNull final String contactNames) {
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
            return contactNames;
        }
        return Strings.EMPTY;
    }
}
