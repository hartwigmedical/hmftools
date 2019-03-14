package com.hartwig.hmftools.common.hospital;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HospitalModelFactory {

    private static final Logger LOGGER = LogManager.getLogger(HospitalModelFactory.class);

    private static final int HOSPITAL_ID_COLUMN = 0;
    private static final int HOSPITAL_COLUMN = 1;
    private static final int CPCT_PI_COLUMN = 3;
    private static final int CPCT_RECIPIENTS_COLUMN = 4;
    private static final int DRUP_PI_COLUMN = 6;
    private static final int DRUP_RECIPIENTS_COLUMN = 7;
    private static final int ADDRESS_NAME_COLUMN = 10;
    private static final int ADDRESS_ZIP_COLUMN = 11;
    private static final int ADDRESS_CITY_COLUMN = 12;

    private static final int FIELD_COUNT = 13;

    private static final int SAMPLE_ID_COLUMN_MANUAL = 0;
    private static final int HOSPITAL_COLUMN_MANUAL = 1;

    private static final int FIELD_COUNT_MANUAL = 2;

    private static final String FIELD_SEPARATOR = ",";

    private HospitalModelFactory() {
    }

    @NotNull
    public static HospitalModel readFromCSV(@NotNull final String pathToCsv, @NotNull final String pathToCsvManual) throws IOException {
        final Map<String, HospitalData> hospitalPerId = Maps.newHashMap();
        final Map<String, HospitalData> hospitalPerHospital = Maps.newHashMap();

        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                HospitalData hospital = ImmutableHospitalData.of(parts[CPCT_RECIPIENTS_COLUMN],
                        parts[DRUP_RECIPIENTS_COLUMN],
                        parts[ADDRESS_NAME_COLUMN],
                        parts[ADDRESS_ZIP_COLUMN],
                        parts[ADDRESS_CITY_COLUMN],
                        parts[CPCT_PI_COLUMN],
                        parts[DRUP_PI_COLUMN]);

                hospitalPerId.put(parts[HOSPITAL_ID_COLUMN], hospital);
                hospitalPerHospital.put(parts[HOSPITAL_COLUMN], hospital);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: " + line);
            }
        }
        return ImmutableHospitalModel.of(hospitalPerId, hospitalPerHospital, readFromCSVManual(pathToCsvManual));
    }

    @NotNull
    private static Map<String, HospitalSampleMapping>  readFromCSVManual(@NotNull final String pathToCsv) throws IOException {
        final Map<String, HospitalSampleMapping> hospitalPerIdManual = Maps.newHashMap();

        final List<String> lines = FileReader.build().readLines(new File(pathToCsv).toPath());
        for (final String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT_MANUAL);
            if (parts.length == FIELD_COUNT_MANUAL) {
                HospitalSampleMapping hospitalManual = ImmutableHospitalSampleMapping.of(parts[HOSPITAL_COLUMN_MANUAL]);

                hospitalPerIdManual.put(parts[SAMPLE_ID_COLUMN_MANUAL], hospitalManual);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in hospital csv: " + line);
            }
        }
        return hospitalPerIdManual;
    }
}
