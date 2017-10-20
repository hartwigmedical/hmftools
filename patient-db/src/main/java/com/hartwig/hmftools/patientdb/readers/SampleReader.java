package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.lims.LimsCsv;
import com.hartwig.hmftools.common.lims.LimsData;
import com.hartwig.hmftools.common.lims.LimsJsonModel;
import com.hartwig.hmftools.common.lims.LimsModel;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class SampleReader {
    private static final Logger LOGGER = LogManager.getLogger(SampleReader.class);

    private static final DateTimeFormatter LIMS_DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final LimsJsonModel limsModel;
    @NotNull
    private final LimsModel limsOldModel;
    @NotNull
    private final LimsModel limsUmcuModel;

    SampleReader(@NotNull final String limsJson, @NotNull final String limsOldCsv, @NotNull final String limsUmcuCsv)
            throws IOException, EmptyFileException {
        LOGGER.info("Reading lims file: " + limsJson);
        this.limsModel = LimsJsonModel.readModelFromFile(limsJson);
        LOGGER.info("Reading lims file: " + limsOldCsv);
        this.limsOldModel = LimsCsv.buildModelFromCsv(limsOldCsv, LIMS_DATE_FORMATTER);
        LOGGER.info("Reading lims file: " + limsUmcuCsv);
        this.limsUmcuModel = LimsCsv.buildModelFromCsv(limsUmcuCsv, LIMS_DATE_FORMATTER);
    }

    @NotNull
    List<SampleData> read(@NotNull final List<String> sampleIds) {
        final List<SampleData> limsBiopsies = Lists.newArrayList();
        sampleIds.forEach(sampleId -> {
            final LimsData oldLimsData = limsOldModel.findDataPerSample(sampleId);
            final LimsData umcuLimsData = limsUmcuModel.findDataPerSample(sampleId);
            if (limsModel.dataPerSample().containsKey(sampleId) || Utils.anyNotNull(oldLimsData, umcuLimsData)) {
                final LocalDate arrivalDate = determineArrivalDate(limsModel.arrivalDateForSample(sampleId), oldLimsData, umcuLimsData);
                final LocalDate samplingDate = determineSamplingDate(limsModel.samplingDateForSample(sampleId), oldLimsData, umcuLimsData);
                limsBiopsies.add(new SampleData(sampleId, arrivalDate, samplingDate));
            } else {
                LOGGER.warn("Missing LIMS data for sample: " + sampleId);
            }
        });
        return limsBiopsies;
    }

    @NotNull
    private static LocalDate determineArrivalDate(@Nullable final LocalDate hmfLimsArrivalDate, @Nullable final LimsData oldLimsData,
            @Nullable final LimsData umcuLimsData) {
        if (hmfLimsArrivalDate != null) {
            return hmfLimsArrivalDate;
        } else {
            assert oldLimsData != null;
            final LocalDate hmfArrivalDate = oldLimsData.arrivalDate();
            final LocalDate umcuArrivalDate = umcuLimsData == null ? null : umcuLimsData.arrivalDate();
            return umcuArrivalDate == null ? hmfArrivalDate : umcuArrivalDate.isBefore(hmfArrivalDate) ? umcuArrivalDate : hmfArrivalDate;
        }
    }

    @Nullable
    private static LocalDate determineSamplingDate(@Nullable final LocalDate hmfLimsSamplingDate, @Nullable final LimsData oldLimsData,
            @Nullable final LimsData umcuLimsData) {
        if (hmfLimsSamplingDate != null) {
            return hmfLimsSamplingDate;
        } else {
            final LocalDate hmfSamplingDate = oldLimsData == null ? null : oldLimsData.samplingDate();
            final LocalDate umcuSamplingDate = umcuLimsData == null ? null : umcuLimsData.samplingDate();
            return umcuSamplingDate == null
                    ? hmfSamplingDate
                    : (hmfSamplingDate == null
                            ? umcuSamplingDate
                            : umcuSamplingDate.isBefore(hmfSamplingDate) ? umcuSamplingDate : hmfSamplingDate);
        }
    }
}
