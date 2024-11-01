package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CUPPA;

import java.io.IOException;
import java.time.LocalDateTime;

import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.CuppaRecord;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep7;

public class CuppaDAO
{

    @NotNull
    private final DSLContext context;

    CuppaDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    private static Double parseDouble(Double value)
    {
        if(Double.isNaN(value))
        {
            return null;
        }

        if(value.equals(Double.POSITIVE_INFINITY))
        {
            return Double.MAX_VALUE;
        }

        if(value.equals(Double.NEGATIVE_INFINITY))
        {
            return Double.MIN_VALUE;
        }

        return value;
    }

    private static String parseClfName(Categories.ClfName category)
    {
        if(category.equals(Categories.ClfName.NONE))
        {
            return null;
        }
        return category.toString();
    }

    void deleteCuppaForSample(@NotNull String sample)
    {
        context.delete(CUPPA).where(CUPPA.SAMPLEID.eq(sample)).execute();
    }

    void writeCuppa2(
            @NotNull final String sample,
            @NotNull CuppaPredictions cuppaPredictions,
            @NotNull final int topNProbs
    ) throws IOException
    {
        deleteCuppaForSample(sample);

        @NotNull
        InsertValuesStep7<CuppaRecord, LocalDateTime, String, String, String, Double, Integer, Byte> inserter = context.insertInto(CUPPA,
                CUPPA.MODIFIED,
                CUPPA.SAMPLEID,
                CUPPA.CLFNAME,
                CUPPA.CANCERTYPE,
                CUPPA.PROB,
                CUPPA.RANK,
                CUPPA.ISOLDCUPPAOUTPUT
        );

        LocalDateTime timestamp = LocalDateTime.now();

        cuppaPredictions = cuppaPredictions
                .subsetByDataType(Categories.DataType.PROB)
                .getTopPredictions(topNProbs)
                .sortByRank();

        for(CuppaPredictionEntry cuppaPredictionEntry : cuppaPredictions.PredictionEntries)
        {
            inserter.values(
                    timestamp,
                    sample,
                    parseClfName(cuppaPredictionEntry.ClfName),
                    cuppaPredictionEntry.CancerType,
                    parseDouble(cuppaPredictionEntry.DataValue),
                    cuppaPredictionEntry.Rank,
                    (byte) 0
            );
        }

        inserter.execute();
    }

    void writeCuppa(@NotNull String sample, @NotNull String cancerType, double likelihood)
    {
        deleteCuppaForSample(sample);
        LocalDateTime timestamp = LocalDateTime.now();

        context.insertInto(CUPPA,
                CUPPA.MODIFIED,
                CUPPA.SAMPLEID,
                CUPPA.CLFNAME,
                CUPPA.CANCERTYPE,
                CUPPA.PROB,
                CUPPA.RANK,
                CUPPA.ISOLDCUPPAOUTPUT
        ).values(
                timestamp,
                sample,
                null,
                cancerType,
                likelihood,
                1,
                (byte) 1
        ).execute();
    }
}