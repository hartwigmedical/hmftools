package com.hartwig.hmftools.common.centra;

import java.util.List;
import java.util.Map;

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
public abstract class CentraModel {
    private static final Logger LOGGER = LogManager.getLogger(CentraModel.class);

    @NotNull
    protected abstract Map<String, CentraData> centraPerId();

    @Nullable
    String getCpctRecipients(@NotNull final String centraId) {
        final CentraData centra = centraPerId().get(centraId);
        if (centra == null) {
            LOGGER.error("Centra model does not contain id " + centraId);
            return null;
        }
        return centra.cpctRecipients();
    }

    @Nullable
    public String getCpctRecipientsFromSample(@NotNull final String sample) {
        return getCpctRecipients(getCentraIdFromSample(sample));
    }

    @Nullable
    String getDrupRecipients(@NotNull final String centraId) {
        final CentraData centra = centraPerId().get(centraId);
        if (centra == null) {
            LOGGER.error("Centra model does not contain id " + centraId);
            return null;
        }
        return centra.drupRecipients();
    }

    @Nullable
    public String getDrupRecipientsFromSample(@NotNull final String sample) {
        return getDrupRecipients(getCentraIdFromSample(sample));
    }

    @Nullable
    CentraData centraPerId(@NotNull final String centraId) {
        return centraPerId().get(centraId);
    }

    // MIVO: expects sample in DRUP/CPCT format
    @NotNull
    private static String getCentraIdFromSample(@NotNull final String sample) {
        final String ucSample = sample.toUpperCase();
        if ((ucSample.startsWith("DRUP") || ucSample.startsWith("CPCT")) && sample.length() >= 12) {
            return sample.substring(6, 8);
        }
        throw new IllegalArgumentException("Sample parameter: " + sample + "was not in CPCT/DRUP format");
    }

    @Nullable
    public String getAddressStringForSample(@NotNull final String sample) {
        final String centraId = getCentraIdFromSample(sample);
        final CentraData centra = centraPerId(centraId);
        if (centra == null) {
            LOGGER.error("Centra model does not contain id " + centraId);
            return null;
        }
        checkAddressFields(centra);
        return centra.addressName() + ", " + centra.addressZip() + " " + centra.addressCity();
    }

    private void checkAddressFields(@NotNull final CentraData centra) {
        final List<String> missingFields = Lists.newArrayList();
        if (centra.addressName().isEmpty()) {
            missingFields.add("name");
        }

        if (centra.addressZip().isEmpty()) {
            missingFields.add("zip");
        }

        if (centra.addressCity().isEmpty()) {
            missingFields.add("city");
        }
        if (!missingFields.isEmpty()) {
            LOGGER.warn("Some address fields (" + Strings.join(missingFields, ',') + ") are missing.");
        }
    }
}
