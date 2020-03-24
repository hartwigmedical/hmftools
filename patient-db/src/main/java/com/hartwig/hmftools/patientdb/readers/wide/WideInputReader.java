package com.hartwig.hmftools.patientdb.readers.wide;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class WideInputReader {

    private static final Logger LOGGER = LogManager.getLogger(WideInputReader.class);
    private static final String FIELD_SEPARATOR = ";";

    private static final int TREATMENT_DATA_COUNT = 5;
    private static final int TREATMENT_DATA_ID = 0;
    private static final int TREATMENT_DATA_CODE = 1;
    private static final int TREATMENT_DATA_DRUG = 2;
    private static final int TREATMENT_DATA_START_DATE = 3;
    private static final int TREATMENT_DATA_END_DATE = 4;

    private static final int PRE_TREATMENT_DATA_COUNT = 7;
    private static final int PRE_TREATMENT_DATA_PATIENT_ID = 0;
    private static final int PRE_TREATMENT_DATA_PREVIOUS_THERAPY = 1;
    private static final int PRE_TREATMENT_DATA_DRUG1 = 2;
    private static final int PRE_TREATMENT_DATA_DRUG2 = 3;
    private static final int PRE_TREATMENT_DATA_DRUG3 = 4;
    private static final int PRE_TREATMENT_DATA_DRUG4 = 5;
    private static final int PRE_TREATMENT_DATA_DATE_LAST_SYSTEMIC_THERAPY = 6;

    private static final int BIOPSY_DATA_COUNT = 6;
    private static final int BIOPSY_DATA_PATIENT_ID = 0;
    private static final int BIOPSY_DATA_WIDE_ID = 1;
    private static final int BIOPSY_DATA_DATA_AVAILABLE = 2;
    private static final int BIOPSY_DATA_TISSUE_ID = 3;
    private static final int BIOPSY_DATA_DATE = 4;
    private static final int BIOPSY_DATA_WGS_SUCCESFUL = 4;

    private WideInputReader() {
    }

    @NotNull
    public static Map<Integer, WideTreatmentData> buildTreatmentData(@NotNull final String pathToCsv) throws IOException {

        final Map<Integer, WideTreatmentData> treatmentDataPerEntry = Maps.newHashMap();
        int count = 1;
        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines.subList(1, lines.size())) {
            final String[] parts = line.split(FIELD_SEPARATOR, TREATMENT_DATA_COUNT);
            if (parts.length == TREATMENT_DATA_COUNT) {
                WideTreatmentData wideTreatmentData = ImmutableWideTreatmentData.of(parts[TREATMENT_DATA_ID].replace("-", ""),
                        parts[TREATMENT_DATA_CODE],
                        parts[TREATMENT_DATA_DRUG],
                        parts[TREATMENT_DATA_START_DATE],
                        parts[TREATMENT_DATA_END_DATE]);

                treatmentDataPerEntry.put(count, wideTreatmentData);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE treatment csv: " + line);
            }
            count +=1;
        }
        LOGGER.info(treatmentDataPerEntry);
        return treatmentDataPerEntry;
    }

    @NotNull
    public static Map<Integer, WidePreTreatmentData> buildPreTreatmentData(@NotNull final String pathToCsv) throws IOException {

        final Map<Integer, WidePreTreatmentData> preTreatmentDataPerEntry = Maps.newHashMap();
        int count = 1;

        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines.subList(1, lines.size())) {
            final String[] parts = line.split(FIELD_SEPARATOR, PRE_TREATMENT_DATA_COUNT);
            if (parts.length == PRE_TREATMENT_DATA_COUNT) {
                WidePreTreatmentData widePreTreatmentData = ImmutableWidePreTreatmentData.of(parts[PRE_TREATMENT_DATA_PATIENT_ID].replace("-", ""),
                        parts[PRE_TREATMENT_DATA_PREVIOUS_THERAPY],
                        parts[PRE_TREATMENT_DATA_DRUG1],
                        parts[PRE_TREATMENT_DATA_DRUG2],
                        parts[PRE_TREATMENT_DATA_DRUG3],
                        parts[PRE_TREATMENT_DATA_DRUG4],
                        parts[PRE_TREATMENT_DATA_DATE_LAST_SYSTEMIC_THERAPY]);

                preTreatmentDataPerEntry.put(count, widePreTreatmentData);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE pre treatment data csv: " + line);
            }
            count +=1;
        }
        LOGGER.info(preTreatmentDataPerEntry);
        return preTreatmentDataPerEntry;
    }

    @NotNull
    public static Map<Integer, WideBiopsyData> buildBiopsyData(@NotNull final String pathToCsv) throws IOException {

        final Map<Integer, WideBiopsyData> biopsyDataPerEntry = Maps.newHashMap();
        int count = 1;

        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines.subList(1, lines.size())) {
            final String[] parts = line.split(FIELD_SEPARATOR, BIOPSY_DATA_COUNT);
            if (parts.length == BIOPSY_DATA_COUNT) {
                WideBiopsyData wideBiopsyData = ImmutableWideBiopsyData.of(parts[BIOPSY_DATA_PATIENT_ID],
                        parts[BIOPSY_DATA_WIDE_ID].replace("-", ""),
                        parts[BIOPSY_DATA_DATA_AVAILABLE],
                        parts[BIOPSY_DATA_TISSUE_ID],
                        parts[BIOPSY_DATA_DATE],
                        parts[BIOPSY_DATA_WGS_SUCCESFUL]);

                biopsyDataPerEntry.put(count, wideBiopsyData);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in WIDE biopsy date csv: " + line);
            }
            count +=1;
        }
        LOGGER.info(biopsyDataPerEntry);
        return biopsyDataPerEntry;
    }
}
