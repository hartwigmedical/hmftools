package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.COPYNUMBERCHROMOSOMEARM;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.purple.cnchromosome.CnPerChromosomeArmData;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep5;

public class CnPerChromosomeArmDAO {

    @NotNull
    private final DSLContext context;

    CnPerChromosomeArmDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeChromosomeCopyNumber(@NotNull String sample, @NotNull List<CnPerChromosomeArmData> cnPerChromosomeArmData) {
        deleteCopyNumberChromosomeArmForSample(sample);
        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<CnPerChromosomeArmData> cnPerChromosomeArm : Iterables.partition(cnPerChromosomeArmData, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep5 inserter = context.insertInto(COPYNUMBERCHROMOSOMEARM,
                    COPYNUMBERCHROMOSOMEARM.MODIFIED,
                    COPYNUMBERCHROMOSOMEARM.SAMPLEID,
                    COPYNUMBERCHROMOSOMEARM.CHROMOSOME,
                    COPYNUMBERCHROMOSOMEARM.CHROMOSOMEARM,
                    COPYNUMBERCHROMOSOMEARM.COPYNUMBER);
            cnPerChromosomeArm.forEach(x -> addCopynumberChromosomeRecord(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addCopynumberChromosomeRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep5 inserter,
            @NotNull String sample, @NotNull CnPerChromosomeArmData cnPerChromosomeArmData) {
        inserter.values(timestamp,
                sample,
                cnPerChromosomeArmData.chromosome(),
                cnPerChromosomeArmData.chromosomeArm(),
                cnPerChromosomeArmData.copyNumber());
    }

    void deleteCopyNumberChromosomeArmForSample(@NotNull String sample) {
         context.delete(COPYNUMBERCHROMOSOMEARM).where(COPYNUMBERCHROMOSOMEARM.SAMPLEID.eq(sample)).execute();
    }
}
