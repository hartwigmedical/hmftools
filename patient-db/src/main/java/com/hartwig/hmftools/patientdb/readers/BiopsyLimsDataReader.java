package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsData;
import com.hartwig.hmftools.common.lims.LimsModel;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class BiopsyLimsDataReader {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyLimsDataReader.class);

    private static final DateTimeFormatter LIMS_DATE_FORMATTER = DateTimeFormatter.ofPattern("dd-MM-yyyy");
    private static final DateTimeFormatter OLD_LIMS_DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final LimsModel limsModel;
    @NotNull
    private final LimsModel limsOldModel;
    @NotNull
    private final LimsModel limsUmcuModel;

    BiopsyLimsDataReader(@NotNull final String limsCsv, @NotNull final String limsOldCsv,
            @NotNull final String limsUmcuCsv) throws IOException, EmptyFileException {
        LOGGER.info("Reading lims file: " + limsCsv);
        this.limsModel = Lims.buildModelFromCsv(limsCsv, LIMS_DATE_FORMATTER);
        LOGGER.info("Reading lims file: " + limsOldCsv);
        this.limsOldModel = Lims.buildModelFromCsv(limsOldCsv, OLD_LIMS_DATE_FORMATTER);
        LOGGER.info("Reading lims file: " + limsUmcuCsv);
        this.limsUmcuModel = Lims.buildModelFromCsv(limsUmcuCsv, OLD_LIMS_DATE_FORMATTER);
    }

    @NotNull
    List<BiopsyLimsData> read(@NotNull final List<String> sampleIds) {
        final List<BiopsyLimsData> limsBiopsies = Lists.newArrayList();
        sampleIds.forEach(sampleId -> {
            final LimsData hmfLimsData = limsModel.findDataPerSample(sampleId);
            final LimsData oldLimsData = limsOldModel.findDataPerSample(sampleId);
            final LimsData umcuLimsData = limsUmcuModel.findDataPerSample(sampleId);
            if (Utils.anyNotNull(hmfLimsData, oldLimsData, umcuLimsData)) {
                final LocalDate arrivalDate = determineArrivalDate(hmfLimsData, oldLimsData, umcuLimsData);
                final LocalDate samplingDate = determineSamplingDate(hmfLimsData, oldLimsData, umcuLimsData);
                limsBiopsies.add(new BiopsyLimsData(sampleId, arrivalDate, samplingDate));
            } else {
                LOGGER.warn("Missing LIMS data for sample: " + sampleId);
            }
        });
        return limsBiopsies;
    }

    @NotNull
    private static LocalDate determineArrivalDate(@Nullable final LimsData hmfLimsData,
            @Nullable final LimsData oldLimsData, @Nullable final LimsData umcuLimsData) {
        if (hmfLimsData != null) {
            return hmfLimsData.arrivalDate();
        } else {
            assert oldLimsData != null;
            final LocalDate hmfArrivalDate = oldLimsData.arrivalDate();
            final LocalDate umcuArrivalDate = umcuLimsData == null ? null : umcuLimsData.arrivalDate();
            return umcuArrivalDate == null ?
                    hmfArrivalDate :
                    umcuArrivalDate.isBefore(hmfArrivalDate) ? umcuArrivalDate : hmfArrivalDate;
        }
    }

    @Nullable
    private static LocalDate determineSamplingDate(@Nullable final LimsData hmfLimsData,
            @Nullable final LimsData oldLimsData, @Nullable final LimsData umcuLimsData) {
        final LocalDate limsSamplingDate = hmfLimsData == null ? null : hmfLimsData.samplingDate();
        if (limsSamplingDate != null)
            return limsSamplingDate;
        else {
            final LocalDate hmfSamplingDate = oldLimsData == null ? null : oldLimsData.samplingDate();
            final LocalDate umcuSamplingDate = umcuLimsData == null ? null : umcuLimsData.samplingDate();
            return umcuSamplingDate == null ?
                    hmfSamplingDate :
                    (hmfSamplingDate == null ?
                            umcuSamplingDate :
                            umcuSamplingDate.isBefore(hmfSamplingDate) ? umcuSamplingDate : hmfSamplingDate);
        }
    }
}
