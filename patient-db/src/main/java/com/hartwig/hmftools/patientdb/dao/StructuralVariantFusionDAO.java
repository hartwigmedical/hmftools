package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVFUSION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRALINSERTION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Viralinsertion;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep17;
import org.jooq.InsertValuesStep19;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.types.UInteger;

public class StructuralVariantFusionDAO {

    @NotNull
    private final DSLContext context;

    public StructuralVariantFusionDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void deleteAnnotationsForSample(@NotNull String sampleId) {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    @NotNull
    private InsertValuesStep19 createBreakendInserter() {
        return context.insertInto(SVBREAKEND,
                SVBREAKEND.MODIFIED,
                SVBREAKEND.SAMPLEID,
                SVBREAKEND.SVID,
                SVBREAKEND.STARTBREAKEND,
                SVBREAKEND.GENE,
                SVBREAKEND.TRANSCRIPTID,
                SVBREAKEND.CANONICALTRANSCRIPT,
                SVBREAKEND.GENEORIENTATION,
                SVBREAKEND.DISRUPTIVE,
                SVBREAKEND.REPORTEDDISRUPTION,
                SVBREAKEND.UNDISRUPTEDCOPYNUMBER,
                SVBREAKEND.REGIONTYPE,
                SVBREAKEND.CODINGCONTEXT,
                SVBREAKEND.BIOTYPE,
                SVBREAKEND.EXONICBASEPHASE,
                SVBREAKEND.NEXTSPLICEEXONRANK,
                SVBREAKEND.NEXTSPLICEEXONPHASE,
                SVBREAKEND.NEXTSPLICEDISTANCE,
                SVBREAKEND.TOTALEXONCOUNT);
    }

    public void writeBreakendsAndFusions(@NotNull String sampleId, @NotNull List<LinxBreakend> breakends,
            @NotNull List<LinxFusion> fusions) {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        // a map of breakend DB Ids to transcripts for the fusion DB record foreign key to the breakend table
        Map<Integer, Integer> breakendIdToDbIdMap = Maps.newHashMap();

        InsertValuesStep19 inserter = createBreakendInserter();
        List<LinxBreakend> insertedBreakends = Lists.newArrayList();

        for (int i = 0; i < breakends.size(); ++i) {
            LinxBreakend breakend = breakends.get(i);

            inserter.values(timestamp,
                    sampleId,
                    breakend.svId(),
                    breakend.isStart(),
                    breakend.gene(),
                    breakend.transcriptId(),
                    breakend.canonical(),
                    breakend.geneOrientation(),
                    breakend.disruptive(),
                    breakend.reportedDisruption(),
                    DatabaseUtil.decimal(breakend.undisruptedCopyNumber()),
                    breakend.regionType(),
                    breakend.codingContext(),
                    breakend.biotype(),
                    breakend.exonicBasePhase(),
                    breakend.nextSpliceExonRank(),
                    breakend.nextSpliceExonPhase(),
                    breakend.nextSpliceDistance(),
                    breakend.totalExonCount());

            insertedBreakends.add(breakend);

            // batch-insert transcripts since there can be many more than the batch size per sample
            if (insertedBreakends.size() >= DB_BATCH_INSERT_SIZE || i == breakends.size() - 1) {
                List<UInteger> ids = inserter.returning(SVBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if (ids.size() != insertedBreakends.size()) {
                    throw new RuntimeException("Not all transcripts were inserted successfully");
                }

                for (int j = 0; j < ids.size(); j++) {
                    breakendIdToDbIdMap.put(insertedBreakends.get(j).id(), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                insertedBreakends.clear();
            }
        }

        for (List<LinxFusion> batch : Iterables.partition(fusions, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep17 fusionInserter = context.insertInto(SVFUSION,
                    SVFUSION.MODIFIED,
                    SVFUSION.SAMPLEID,
                    SVFUSION.FIVEPRIMEBREAKENDID,
                    SVFUSION.THREEPRIMEBREAKENDID,
                    SVFUSION.NAME,
                    SVFUSION.REPORTED,
                    SVFUSION.REPORTEDTYPE,
                    SVFUSION.PHASED,
                    SVFUSION.CHAINLENGTH,
                    SVFUSION.CHAINLINKS,
                    SVFUSION.CHAINTERMINATED,
                    SVFUSION.DOMAINSKEPT,
                    SVFUSION.DOMAINSLOST,
                    SVFUSION.SKIPPEDEXONSUP,
                    SVFUSION.SKIPPEDEXONSDOWN,
                    SVFUSION.FUSEDEXONUP,
                    SVFUSION.FUSEDEXONDOWN);

            for (LinxFusion fusion : batch) {
                Integer fivePrimeId = breakendIdToDbIdMap.get(fusion.fivePrimeBreakendId());
                Integer threePrimeId = breakendIdToDbIdMap.get(fusion.threePrimeBreakendId());

                if (fivePrimeId == null || threePrimeId == null) {
                    return;
                }

                fusionInserter.values(timestamp,
                        sampleId,
                        fivePrimeId,
                        threePrimeId,
                        fusion.name(),
                        fusion.reported(),
                        fusion.reportedType(),
                        fusion.phased(),
                        fusion.chainLength(),
                        fusion.chainLinks(),
                        fusion.chainTerminated(),
                        DatabaseUtil.checkStringLength(fusion.domainsKept(), SVFUSION.DOMAINSKEPT),
                        DatabaseUtil.checkStringLength(fusion.domainsLost(), SVFUSION.DOMAINSLOST),
                        fusion.skippedExonsUp(),
                        fusion.skippedExonsDown(),
                        fusion.fusedExonUp(),
                        fusion.fusedExonDown());
            }

            fusionInserter.execute();
        }
    }

    @NotNull
    List<LinxFusion> readFusions(@NotNull String sample)
    {
        List<LinxFusion> fusionList = Lists.newArrayList();

        Result<Record> result = context.select().from(SVFUSION).where(SVFUSION.SAMPLEID.eq(sample)).fetch();

        for (Record record : result)
        {
            LinxFusion fusion = ImmutableLinxFusion.builder()
                    .fivePrimeBreakendId(record.getValue(SVFUSION.FIVEPRIMEBREAKENDID).intValue())
                    .threePrimeBreakendId(record.getValue(SVFUSION.THREEPRIMEBREAKENDID).intValue())
                    .name(record.getValue(SVFUSION.NAME))
                    .reported(record.getValue(SVFUSION.REPORTED) == 1)
                    .reportedType(record.getValue(SVFUSION.REPORTEDTYPE))
                    .phased(record.getValue(SVFUSION.PHASED) == 1)
                    .chainLength(record.getValue(SVFUSION.CHAINLENGTH))
                    .chainLinks(record.getValue(SVFUSION.CHAINLINKS))
                    .chainTerminated(record.getValue(SVFUSION.CHAINTERMINATED) == 1)
                    .domainsKept(record.getValue(SVFUSION.DOMAINSKEPT))
                    .domainsLost(record.getValue(SVFUSION.DOMAINSLOST))
                    .skippedExonsUp(record.getValue(SVFUSION.SKIPPEDEXONSUP))
                    .skippedExonsDown(record.getValue(SVFUSION.SKIPPEDEXONSDOWN))
                    .fusedExonUp(record.getValue(SVFUSION.FUSEDEXONUP))
                    .fusedExonDown(record.getValue(SVFUSION.FUSEDEXONDOWN))
                    .build();

            fusionList.add(fusion);
        }

        return fusionList;
    }

    @NotNull
    List<LinxViralInsertion> readViralInsertions(@NotNull String sample)
    {
        List<LinxViralInsertion> viralinsertions = Lists.newArrayList();

        Result<Record> result = context.select().from(VIRALINSERTION).where(VIRALINSERTION.SAMPLEID.eq(sample)).fetch();

        for (Record record : result)
        {
            LinxViralInsertion viralInsert = new LinxViralInsertion(
                    sample, record.getValue(VIRALINSERTION.SVID),
                    record.getValue(VIRALINSERTION.VIRUSID), record.getValue(VIRALINSERTION.VIRUSNAME));

            viralinsertions.add(viralInsert);
        }

        return viralinsertions;
    }

}

