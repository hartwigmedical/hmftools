package com.hartwig.hmftools.patientdb.dao;

import java.time.LocalDateTime;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CUPPA;

import com.hartwig.hmftools.common.cuppa2.Categories;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictionEntry;
import com.hartwig.hmftools.common.cuppa2.CuppaPredictions;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.CuppaRecord;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;

public class CuppaDAO {

    @NotNull private final DSLContext context;

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

    void deleteCuppaForSample(@NotNull String sample) {
        context.delete(CUPPA).where(CUPPA.SAMPLEID.eq(sample)).execute();
    }

    void writeCuppa2(
            @NotNull final String sample,
            @NotNull final CuppaPredictions cuppaPredictions
    ){
        deleteCuppaForSample(sample);

        @NotNull
        InsertValuesStep8<CuppaRecord, LocalDateTime, String, String, String, String, Double, Integer, Integer> inserter = context.insertInto(CUPPA,
                CUPPA.MODIFIED,
                CUPPA.SAMPLEID,
                CUPPA.DATATYPE,
                CUPPA.CLFNAME,
                CUPPA.CANCERTYPE,
                CUPPA.DATAVALUE,
                CUPPA.RANK,
                CUPPA.RANKGROUP
        );

        LocalDateTime timestamp = LocalDateTime.now();

        for(CuppaPredictionEntry cuppaPredictionEntry : cuppaPredictions.PredictionEntries)
        {
            inserter.values(
                    timestamp,
                    cuppaPredictionEntry.SampleId,
                    cuppaPredictionEntry.DataType.toString(),
                    parseClfName(cuppaPredictionEntry.ClfName),
                    cuppaPredictionEntry.CancerType,
                    parseDouble(cuppaPredictionEntry.DataValue),
                    cuppaPredictionEntry.Rank,
                    cuppaPredictionEntry.RankGroup
            );
        }

        inserter.execute();
    }

    void writeCuppa(@NotNull String sample, @NotNull String cancerType, double likelihood) {
        deleteCuppaForSample(sample);
        LocalDateTime timestamp = LocalDateTime.now();

        context.insertInto(CUPPA,
                        CUPPA.MODIFIED,
                        CUPPA.SAMPLEID,
                        CUPPA.CANCERTYPE,
                        CUPPA.DATAVALUE,
                        CUPPA.DATATYPE
                ).values(
                        timestamp,
                        sample,
                        cancerType,
                        likelihood,
                        Categories.DataType.PROB.toString()
                ).execute();
    }
}
