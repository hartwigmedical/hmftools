package com.hartwig.hmftools.common.center;

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
public abstract class CenterModel {
    private static final Logger LOGGER = LogManager.getLogger(CenterModel.class);

    @NotNull
    protected abstract Map<String, CenterData> centerPerId();

    @NotNull
    protected abstract Map<String, CenterData> centerPerHospital();

    @NotNull
    protected abstract Map<String, CenterDataManualMapping> centerPerIdManual();

    @Nullable
    public String addresseeStringForSample(@NotNull final String contactNames, @NotNull final String sample) {
        String address;
        if (sample.startsWith("CORE19") || sample.contains("CORE18")) { // This are the old core names
            final CenterDataManualMapping centerDataManualMapping = centerPerIdManual().get(sample);
            final CenterData center = centerPerHospital().get(centerDataManualMapping.addressName());
            if (center == null) {
                LOGGER.error("CenterModelFactory model cannot find center details for project " + centerDataManualMapping);
                return null;
            }
            address = center.addressName() + ", " + center.addressZip() + " " + center.addressCity();

        } else {
            final String centerId = getCenterIdFromSample(sample);
            if (centerId == null) {
                return null;
            }
            final CenterData center = centerPerId(centerId);
            if (center == null) {
                LOGGER.error("CenterModelFactory model does not contain id " + centerId);
                return null;
            }
            checkAddresseeFields(sample, center, contactNames);
            address = getPI(sample, center, contactNames) + ", " + center.addressName() + ", " + center.addressZip() + " "
                    + center.addressCity();
        }
        return address;
    }

    @Nullable
    @VisibleForTesting
    CenterData centerPerId(@Nullable final String centerId) {
        return centerPerId().get(centerId);
    }

    @Nullable
    private static String getCenterIdFromSample(@NotNull final String sample) {
        final String ucSample = sample.toUpperCase();
        LimsSampleType type = LimsSampleType.fromSampleId(ucSample);

        if (type == LimsSampleType.DRUP || type == LimsSampleType.CPCT || type == LimsSampleType.WIDE
                || type == LimsSampleType.CORE && sample.length() >= 12) {
            return sample.substring(6, 8);
        }

        LOGGER.warn("Sample parameter: " + sample + " is not in CPCT/DRUP/WIDE/CORE format");
        return null;
    }

    private static void checkAddresseeFields(@NotNull final String sample, @NotNull final CenterData center,
            @NotNull final String contactNames) {
        final List<String> missingFields = Lists.newArrayList();
        if (getPI(sample, center, contactNames).isEmpty()) {
            missingFields.add("PI");
        }
        if (center.addressName().isEmpty()) {
            missingFields.add("name");
        }
        if (center.addressZip().isEmpty()) {
            missingFields.add("zip");
        }
        if (center.addressCity().isEmpty()) {
            missingFields.add("city");
        }
        if (!missingFields.isEmpty()) {
            LOGGER.warn("Some address fields (" + Strings.join(missingFields, ',') + ") are missing.");
        }
    }

    @NotNull
    @VisibleForTesting
    static String getPI(@NotNull final String sample, @NotNull final CenterData center, @NotNull final String contactNames) {
        LimsSampleType type = LimsSampleType.fromSampleId(sample);

        if (type == LimsSampleType.CPCT) {
            return center.cpctPI();
        } else if (type == LimsSampleType.DRUP) {
            final String drupPi = center.drupPI();
            if (drupPi.trim().equals("*")) {
                return center.cpctPI();
            }
            return center.drupPI();
        } else if (type == LimsSampleType.WIDE) {
            return contactNames;
        }
        return Strings.EMPTY;
    }
}
