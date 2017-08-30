package com.hartwig.hmftools.common.center;

import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

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

    @Nullable
    String getCpctRecipients(@NotNull final String centerId) {
        final CenterData center = centerPerId().get(centerId);
        if (center == null) {
            LOGGER.error("Center model does not contain id " + centerId);
            return null;
        }
        return center.cpctRecipients();
    }

    @Nullable
    public String getCpctRecipientsFromSample(@NotNull final String sample) {
        return getCpctRecipients(getCenterIdFromSample(sample));
    }

    @Nullable
    String getDrupRecipients(@NotNull final String centerId) {
        final CenterData center = centerPerId().get(centerId);
        if (center == null) {
            LOGGER.error("Center model does not contain id " + centerId);
            return null;
        }
        final String drupRecipients = center.drupRecipients();
        if (drupRecipients.trim().equals("*")) {
            return center.cpctRecipients();
        }
        return drupRecipients;
    }

    @Nullable
    public String getDrupRecipientsFromSample(@NotNull final String sample) {
        return getDrupRecipients(getCenterIdFromSample(sample));
    }

    @Nullable
    CenterData centerPerId(@NotNull final String centerId) {
        return centerPerId().get(centerId);
    }

    // MIVO: expects sample in DRUP/CPCT format
    @NotNull
    private static String getCenterIdFromSample(@NotNull final String sample) {
        final String ucSample = sample.toUpperCase();
        if ((ucSample.startsWith("DRUP") || ucSample.startsWith("CPCT")) && sample.length() >= 12) {
            return sample.substring(6, 8);
        }
        throw new IllegalArgumentException("Sample parameter: " + sample + "was not in CPCT/DRUP format");
    }

    @Nullable
    public String getAddresseeStringForSample(@NotNull final String sample) {
        final String centerId = getCenterIdFromSample(sample);
        final CenterData center = centerPerId(centerId);
        if (center == null) {
            LOGGER.error("Center model does not contain id " + centerId);
            return null;
        }
        checkAddresseeFields(sample, center);
        return getPI(sample, center) + ", " + center.addressName() + ", " + center.addressZip() + " " + center.addressCity();
    }

    private static void checkAddresseeFields(@NotNull final String sample, @NotNull final CenterData center) {
        final List<String> missingFields = Lists.newArrayList();
        if (getPI(sample, center).isEmpty()) {
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
    static String getPI(@NotNull final String sample, @NotNull final CenterData center) {
        if (sample.toUpperCase().startsWith("CPCT")) {
            return center.cpctPI();
        } else if (sample.toUpperCase().startsWith("DRUP")) {
            final String drupPi = center.drupPI();
            if (drupPi.trim().equals("*")) {
                return center.cpctPI();
            }
            return center.drupPI();
        }
        return "";
    }
}
