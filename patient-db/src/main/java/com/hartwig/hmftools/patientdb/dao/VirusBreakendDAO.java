package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRUSBREAKEND;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.virus.VirusBreakend;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

class VirusBreakendDAO
{
    private final DSLContext context;

    VirusBreakendDAO(final DSLContext context)
    {
        this.context = context;
    }

    void writeVirusBreakend(final String sample, final List<VirusBreakend> virusBreakends)
    {
        deleteVirusBreakendForSample(sample);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for(List<VirusBreakend> virusBreakend : Iterables.partition(virusBreakends, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStepN inserter = context.insertInto(VIRUSBREAKEND,
                    VIRUSBREAKEND.MODIFIED,
                    VIRUSBREAKEND.SAMPLEID,
                    VIRUSBREAKEND.TAXIDGENUS,
                    VIRUSBREAKEND.NAMEGENUS,
                    VIRUSBREAKEND.READSGENUSTREE,
                    VIRUSBREAKEND.TAXIDSPECIES,
                    VIRUSBREAKEND.NAMESPECIES,
                    VIRUSBREAKEND.READSSPECIESTREE,
                    VIRUSBREAKEND.TAXIDASSIGNED,
                    VIRUSBREAKEND.NAMEASSIGNED,
                    VIRUSBREAKEND.READSASSIGNEDTREE,
                    VIRUSBREAKEND.READSASSIGNEDDIRECT,
                    VIRUSBREAKEND.REFERENCE,
                    VIRUSBREAKEND.REFERENCETAXID,
                    VIRUSBREAKEND.REFERENCEKMERCOUNT,
                    VIRUSBREAKEND.ALTERNATEKMERCOUNT,
                    VIRUSBREAKEND.RNAME,
                    VIRUSBREAKEND.STARTPOS,
                    VIRUSBREAKEND.ENDPOS,
                    VIRUSBREAKEND.NUMREADS,
                    VIRUSBREAKEND.COVBASES,
                    VIRUSBREAKEND.COVERAGE,
                    VIRUSBREAKEND.MEANDEPTH,
                    VIRUSBREAKEND.MEANBASEQ,
                    VIRUSBREAKEND.MEANMAPQ,
                    VIRUSBREAKEND.INTEGRATIONS,
                    VIRUSBREAKEND.QCSTATUS);
            virusBreakend.forEach(x -> addVirusBreakend(timestamp, inserter, sample, x));
            inserter.execute();
        }
    }

    private static void addVirusBreakend(
            final Timestamp timestamp, final InsertValuesStepN inserter, final String sample, final VirusBreakend virusBreakend)
    {
        inserter.values(timestamp,
                sample,
                virusBreakend.taxidGenus(),
                virusBreakend.nameAssigned(),
                virusBreakend.readsGenusTree(),
                virusBreakend.taxidSpecies(),
                virusBreakend.nameSpecies(),
                virusBreakend.readsSpeciesTree(),
                virusBreakend.taxidAssigned(),
                virusBreakend.nameAssigned(),
                virusBreakend.readsAssignedTree(),
                virusBreakend.readsAssignedDirect(),
                virusBreakend.reference(),
                virusBreakend.referenceTaxid(),
                virusBreakend.referenceKmerCount(),
                virusBreakend.alternateKmerCount(),
                virusBreakend.RName(),
                virusBreakend.startPos(),
                virusBreakend.endPos(),
                virusBreakend.numReads(),
                virusBreakend.covBases(),
                virusBreakend.coverage(),
                virusBreakend.meanDepth(),
                virusBreakend.meanBaseQ(),
                virusBreakend.meanMapQ(),
                virusBreakend.integrations(),
                virusBreakend.qcStatus());
    }

    void deleteVirusBreakendForSample(final String sample)
    {
        context.delete(VIRUSBREAKEND).where(VIRUSBREAKEND.SAMPLEID.eq(sample)).execute();
    }
}
