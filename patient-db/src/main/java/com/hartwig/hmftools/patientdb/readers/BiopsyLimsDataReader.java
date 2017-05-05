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

import org.jetbrains.annotations.NotNull;

public class BiopsyLimsDataReader {

    private static final DateTimeFormatter newLimsDateFormatter = DateTimeFormatter.ofPattern("dd-MM-yyyy");
    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final LimsModel limsModel;
    @NotNull
    private final LimsModel limsOldModel;

    BiopsyLimsDataReader(@NotNull final String limsCsv, @NotNull final String limsOldCsv)
            throws IOException, EmptyFileException {
        this.limsModel = Lims.buildOldModelFromCsv(limsCsv);
        this.limsOldModel = Lims.buildOldModelFromCsv(limsOldCsv);
    }

    @NotNull
    public List<BiopsyLimsData> read(@NotNull final List<String> sampleIds) {
        return sampleIds.stream().map(
                sampleId -> new BiopsyLimsData(sampleId, determineArrivalDate(sampleId))).collect(Collectors.toList());
    }

    @NotNull
    private LocalDate determineArrivalDate(@NotNull final String sampleId) {
        final LocalDate limsArrivalDate = Utils.getDate(limsModel.findArrivalDateForSample(sampleId),
                newLimsDateFormatter);
        if (limsArrivalDate != null)
            return limsArrivalDate;
        else
            return LocalDate.parse(limsOldModel.findArrivalDateForSample(sampleId), dateFormatter);
    }
}
