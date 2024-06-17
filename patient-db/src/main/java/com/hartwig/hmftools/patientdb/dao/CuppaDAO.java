package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CUPPA;

import java.time.LocalDateTime;

import com.hartwig.hmftools.common.cuppa.ClassifierName;
import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.cuppa.DataType;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.CuppaRecord;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep7;

public class CuppaDAO
{
    private final DSLContext context;

    CuppaDAO(final DSLContext context)
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

    private static String parseClfName(ClassifierName category)
    {
        if(category.equals(ClassifierName.NONE))
        {
            return null;
        }
        return category.toString();
    }

    void deleteCuppaForSample(final String sample)
    {
        context.delete(CUPPA).where(CUPPA.SAMPLEID.eq(sample)).execute();
    }

    void writeCuppa2(final String sample, final CuppaPredictions cuppaPredictions, final int topNProbs)
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

        CuppaPredictions cuppaPredictionsSorted = cuppaPredictions
                .subsetByDataType(DataType.PROB)
                .getTopPredictions(topNProbs)
                .sortByRank();

        for(CuppaPredictionEntry cuppaPredictionEntry : cuppaPredictionsSorted.PredictionEntries)
        {
            inserter.values(
                    timestamp,
                    cuppaPredictionEntry.SampleId,
                    parseClfName(cuppaPredictionEntry.ClassifierName),
                    cuppaPredictionEntry.CancerType,
                    parseDouble(cuppaPredictionEntry.DataValue),
                    cuppaPredictionEntry.Rank,
                    (byte) 0
            );
        }

        inserter.execute();
    }

    void writeCuppa(final String sample, final String cancerType, double likelihood)
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