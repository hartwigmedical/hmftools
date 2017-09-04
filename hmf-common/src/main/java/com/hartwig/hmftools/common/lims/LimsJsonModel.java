package com.hartwig.hmftools.common.lims;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsJsonModel {
    private static final Logger LOGGER = LogManager.getLogger(LimsJsonModel.class);
    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    public abstract Multimap<String, LimsJsonData> dataPerPatient();

    public static LimsJsonModel readModelFromFile(@NotNull final String limsJson) throws FileNotFoundException {
        final Gson gson = LimsGsonAdapter.buildGson();
        final JsonObject jsonObject = new JsonParser().parse(new FileReader(limsJson)).getAsJsonObject();
        final Set<Map.Entry<String, JsonElement>> jsonSamples = jsonObject.getAsJsonObject("samples").entrySet();
        final Multimap<String, LimsJsonData> limsDataPerPatient = ArrayListMultimap.create();
        jsonSamples.forEach(jsonSample -> {
            final JsonObject jsonSampleObject = jsonSample.getValue().getAsJsonObject();
            final String sampleLabel = jsonSampleObject.get("label").getAsString();
            if (sampleLabel.equals("CPCT") || sampleLabel.equals("DRUP")) {
                final LimsJsonData limsJsonData = gson.fromJson(jsonSample.getValue(), LimsJsonData.class);
                limsDataPerPatient.put(limsJsonData.patient(), limsJsonData);
            }
        });
        return ImmutableLimsJsonModel.of(limsDataPerPatient);
    }

    @Nullable
    public String barcodeForSample(@NotNull final String sample) {
        final Collection<LimsJsonData> dataPerPatient = patientDataForSample(sample);
        for (final LimsJsonData data : dataPerPatient) {
            if (data.sampleName().equals(sample)) {
                return data.sampleBarcode().toUpperCase();
            }
        }
        LOGGER.warn("Could not find barcode for sample: " + sample + " in Lims");
        return null;
    }

    @Nullable
    public String bloodBarcodeForSample(@NotNull final String sample) {
        final Collection<LimsJsonData> dataPerPatient = patientDataForSample(sample);
        for (final LimsJsonData data : dataPerPatient) {
            if (data.sampleSource().toLowerCase().equals("blood")) {
                return data.sampleBarcode().toUpperCase();
            }
        }
        LOGGER.warn("Could not find blood barcode for sample: " + sample + " in Lims");
        return null;
    }

    @Nullable
    public LocalDate arrivalDateForSample(@NotNull final String sample) {
        final Collection<LimsJsonData> dataPerPatient = patientDataForSample(sample);
        for (final LimsJsonData data : dataPerPatient) {
            if (data.sampleName().equals(sample)) {
                final LocalDate arrivalDate = getNullableDate(data.arrivalDateString(), DATE_FORMATTER);
                if (arrivalDate == null) {
                    LOGGER.warn("Lims arrival date: " + data.arrivalDateString() + " is not a valid date.");
                }
                return arrivalDate;
            }
        }
        LOGGER.warn("Could not find arrival date for sample: " + sample + " in Lims");
        return null;
    }

    @Nullable
    public LocalDate bloodArrivalDateForSample(@NotNull final String sample) {
        final Collection<LimsJsonData> dataPerPatient = patientDataForSample(sample);
        for (final LimsJsonData data : dataPerPatient) {
            if (data.sampleSource().toLowerCase().equals("blood")) {
                final LocalDate arrivalDate = getNullableDate(data.arrivalDateString(), DATE_FORMATTER);
                if (arrivalDate == null) {
                    LOGGER.warn("Lims arrival date: " + data.arrivalDateString() + " is not a valid date.");
                }
                return arrivalDate;
            }
        }
        LOGGER.warn("Could not find blood arrival date for sample: " + sample + " in Lims");
        return null;
    }

    @Nullable
    public LocalDate samplingDateForSample(@NotNull final String sample) {
        final Collection<LimsJsonData> dataPerPatient = patientDataForSample(sample);
        for (final LimsJsonData data : dataPerPatient) {
            if (data.sampleName().equals(sample)) {
                final LocalDate samplingDate = getNullableDate(data.samplingDateString(), DATE_FORMATTER);
                if (samplingDate == null) {
                    LOGGER.warn("Lims sampling date: " + data.samplingDateString() + " is not a valid date.");
                }
                return samplingDate;
            }
        }
        LOGGER.warn("Could not find sampling date for sample: " + sample + " in Lims");
        return null;
    }

    @Nullable
    public Double tumorPercentageForSample(@NotNull final String sample) {
        final Collection<LimsJsonData> dataPerPatient = patientDataForSample(sample);
        for (final LimsJsonData data : dataPerPatient) {
            if (data.sampleName().equals(sample)) {
                if (data.sampleSource().equals("Blood")) {
                    return null;
                } else {
                    return Double.parseDouble(data.tumorPercentage()) / 100D;
                }
            }
        }
        LOGGER.warn("Could not find tumor percentage for sample: " + sample + " in Lims");
        return null;
    }

    @NotNull
    private Collection<LimsJsonData> patientDataForSample(@NotNull final String sample) {
        final String patientId = getPatientIdFromSample(sample);
        if (patientId == null) {
            LOGGER.error("could not retrieve patientId from sample: " + sample);
        }
        return dataPerPatient().get(patientId);
    }

    @Nullable
    private static String getPatientIdFromSample(@NotNull final String sample) {
        return sample.length() >= 12 ? sample.substring(0, 12) : null;
    }

    @NotNull
    public static LimsJsonModel buildEmptyModel() {
        return ImmutableLimsJsonModel.of(ArrayListMultimap.create());
    }

    @Nullable
    private static LocalDate getNullableDate(@Nullable final String dateString, @NotNull final DateTimeFormatter dateFormatter) {
        if (dateString == null) {
            return null;
        }
        try {
            return LocalDate.parse(dateString, dateFormatter);
        } catch (DateTimeParseException e) {
            return null;
        }
    }
}
