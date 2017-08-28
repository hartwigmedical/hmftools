package com.hartwig.hmftools.common.centra;

import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
        return getCpctRecipients(getCentraFromSample(sample));
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
        return getDrupRecipients(getCentraFromSample(sample));
    }

    @Nullable
    CentraData centraPerId(@NotNull final String centraId) {
        return centraPerId().get(centraId);
    }

    // MIVO: expects sample in DRUP/CPCT format
    @NotNull
    private static String getCentraFromSample(@NotNull final String sample) {
        final String ucSample = sample.toUpperCase();
        if ((ucSample.startsWith("DRUP") || ucSample.startsWith("CPCT")) && sample.length() >= 12) {
            return sample.substring(6, 8);
        }
        throw new IllegalArgumentException("Sample parameter: " + sample + "was not in CPCT/DRUP format");
    }
}
