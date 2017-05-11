package com.hartwig.hmftools.patientdb.readers;

import java.io.IOException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsModel;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsyLimsDataReader {
    private static final Logger LOGGER = LogManager.getLogger(BiopsyLimsDataReader.class);

    private static final DateTimeFormatter newLimsDateFormatter = DateTimeFormatter.ofPattern("dd-MM-yyyy");
    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final LimsModel limsModel;
    @NotNull
    private final LimsModel limsOldModel;
    @NotNull
    private final LimsModel limsUmcuModel;

    BiopsyLimsDataReader(@NotNull final String limsCsv, @NotNull final String limsOldCsv,
            @NotNull final String limsUmcuCsv) throws IOException, EmptyFileException {
        LOGGER.info("reading lims file: " + limsCsv);
        this.limsModel = Lims.buildOldModelFromCsv(limsCsv);
        LOGGER.info("reading lims file: " + limsOldCsv);
        this.limsOldModel = Lims.buildOldModelFromCsv(limsOldCsv);
        LOGGER.info("reading lims file: " + limsUmcuCsv);
        this.limsUmcuModel = Lims.buildOldModelFromCsv(limsUmcuCsv);
    }

    @NotNull
    public List<BiopsyLimsData> read(@NotNull final List<String> sampleIds) {
        return sampleIds.stream().map(sampleId -> new BiopsyLimsData(sampleId, determineArrivalDate(sampleId),
                determineSamplingDate(sampleId))).collect(Collectors.toList());
    }

    @NotNull
    private LocalDate determineArrivalDate(@NotNull final String sampleId) {
        final LocalDate limsArrivalDate = Utils.getDate(limsModel.findArrivalDateForSample(sampleId),
                newLimsDateFormatter);
        if (limsArrivalDate != null)
            return limsArrivalDate;
        else {
            final LocalDate umcuArrivalDate = Utils.getDate(limsUmcuModel.findArrivalDateForSample(sampleId),
                    dateFormatter);
            final LocalDate hmfArrivalDate = LocalDate.parse(limsOldModel.findArrivalDateForSample(sampleId),
                    dateFormatter);
            return umcuArrivalDate == null ?
                    hmfArrivalDate :
                    umcuArrivalDate.isBefore(hmfArrivalDate) ? umcuArrivalDate : hmfArrivalDate;
        }
    }

    @Nullable
    private LocalDate determineSamplingDate(@NotNull final String sampleId) {
        final LocalDate limsSamplingDate = Utils.getDate(limsModel.findSamplingDateForSample(sampleId),
                newLimsDateFormatter);
        if (limsSamplingDate != null)
            return limsSamplingDate;
        else {
            final LocalDate umcuSamplingDate = Utils.getDate(limsUmcuModel.findSamplingDateForSample(sampleId),
                    dateFormatter);
            final LocalDate hmfSamplingDate = Utils.getDate(limsOldModel.findSamplingDateForSample(sampleId),
                    dateFormatter);
            return umcuSamplingDate == null ?
                    hmfSamplingDate :
                    (hmfSamplingDate == null ?
                            umcuSamplingDate :
                            umcuSamplingDate.isBefore(hmfSamplingDate) ? umcuSamplingDate : hmfSamplingDate);
        }
    }
}
