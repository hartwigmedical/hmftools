package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVBREAKEND;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVFUSION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.FusionPhasedType;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep18;
import org.jooq.InsertValuesStep19;
import org.jooq.InsertValuesStep21;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.types.UInteger;

public class StructuralVariantFusionDAO
{
    @NotNull
    private final DSLContext context;

    public StructuralVariantFusionDAO(final DSLContext context)
    {
        this.context = context;
    }

    public void deleteAnnotationsForSample(final String sampleId)
    {
        context.delete(SVFUSION).where(SVFUSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sampleId)).execute();
    }

    @NotNull
    private InsertValuesStep21 createBreakendInserter()
    {
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
                SVBREAKEND.EXONUP,
                SVBREAKEND.EXONDOWN,
                SVBREAKEND.EXONICBASEPHASE,
                SVBREAKEND.NEXTSPLICEEXONRANK,
                SVBREAKEND.NEXTSPLICEEXONPHASE,
                SVBREAKEND.NEXTSPLICEDISTANCE,
                SVBREAKEND.TOTALEXONCOUNT
                );
    }

    public void writeBreakendsAndFusions(final String sampleId, final List<LinxBreakend> breakends, final List<LinxFusion> fusions)
    {
        deleteAnnotationsForSample(sampleId);

        Timestamp timestamp = new Timestamp(new Date().getTime());

        // a map of breakend DB Ids to transcripts for the fusion DB record foreign key to the breakend table
        Map<Integer, Integer> breakendIdToDbIdMap = Maps.newHashMap();

        InsertValuesStep21 inserter = createBreakendInserter();
        List<LinxBreakend> insertedBreakends = Lists.newArrayList();

        for(int i = 0; i < breakends.size(); ++i)
        {
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
                    breakend.codingType(),
                    breakend.biotype(),
                    breakend.exonUp(),
                    breakend.exonDown(),
                    breakend.exonicBasePhase(),
                    breakend.nextSpliceExonRank(),
                    breakend.nextSpliceExonPhase(),
                    breakend.nextSpliceDistance(),
                    breakend.totalExonCount()
                    );

            insertedBreakends.add(breakend);

            // batch-insert transcripts since there can be many more than the batch size per sample
            if(insertedBreakends.size() >= DB_BATCH_INSERT_SIZE || i == breakends.size() - 1)
            {
                List<UInteger> ids = inserter.returning(SVBREAKEND.ID).fetch().getValues(0, UInteger.class);

                if(ids.size() != insertedBreakends.size())
                {
                    throw new RuntimeException("Not all transcripts were inserted successfully");
                }

                for(int j = 0; j < ids.size(); j++)
                {
                    breakendIdToDbIdMap.put(insertedBreakends.get(j).id(), ids.get(j).intValue());
                }

                inserter = createBreakendInserter();
                insertedBreakends.clear();
            }
        }

        for(List<LinxFusion> batch : Iterables.partition(fusions, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep19 fusionInserter = context.insertInto(SVFUSION,
                    SVFUSION.MODIFIED,
                    SVFUSION.SAMPLEID,
                    SVFUSION.FIVEPRIMEBREAKENDID,
                    SVFUSION.THREEPRIMEBREAKENDID,
                    SVFUSION.NAME,
                    SVFUSION.REPORTED,
                    SVFUSION.REPORTEDTYPE,
                    SVFUSION.REPORTEDREASON,
                    SVFUSION.PHASED,
                    SVFUSION.LIKELIHOOD,
                    SVFUSION.CHAINLENGTH,
                    SVFUSION.CHAINLINKS,
                    SVFUSION.CHAINTERMINATED,
                    SVFUSION.DOMAINSKEPT,
                    SVFUSION.DOMAINSLOST,
                    SVFUSION.SKIPPEDEXONSUP,
                    SVFUSION.SKIPPEDEXONSDOWN,
                    SVFUSION.FUSEDEXONUP,
                    SVFUSION.FUSEDEXONDOWN);

            for(LinxFusion fusion : batch)
            {
                Integer fivePrimeId = breakendIdToDbIdMap.get(fusion.fivePrimeBreakendId());
                Integer threePrimeId = breakendIdToDbIdMap.get(fusion.threePrimeBreakendId());

                if(fivePrimeId == null || threePrimeId == null)
                {
                    return;
                }

                fusionInserter.values(timestamp,
                        sampleId,
                        fivePrimeId,
                        threePrimeId,
                        fusion.name(),
                        fusion.reported(),
                        fusion.reportedType(),
                        DatabaseUtil.checkStringLength(fusion.reportableReasons(), SVFUSION.REPORTEDREASON),
                        fusion.phased().toString(),
                        fusion.likelihood().toString(),
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
    List<LinxFusion> readFusions(final String sample)
    {
        List<LinxFusion> fusionList = Lists.newArrayList();

        Result<Record> result = context.select().from(SVFUSION).where(SVFUSION.SAMPLEID.eq(sample)).fetch();

        for(Record record : result)
        {
            LinxFusion fusion = ImmutableLinxFusion.builder()
                    .fivePrimeBreakendId(record.getValue(SVFUSION.FIVEPRIMEBREAKENDID).intValue())
                    .threePrimeBreakendId(record.getValue(SVFUSION.THREEPRIMEBREAKENDID).intValue())
                    .name(record.getValue(SVFUSION.NAME))
                    .reported(record.getValue(SVFUSION.REPORTED) == 1)
                    .reportedType(record.getValue(SVFUSION.REPORTEDTYPE))
                    .likelihood(FusionLikelihoodType.valueOf(record.getValue(SVFUSION.LIKELIHOOD)))
                    .phased(FusionPhasedType.valueOf(record.getValue(SVFUSION.PHASED)))
                    .chainLength(record.getValue(SVFUSION.CHAINLENGTH))
                    .chainLinks(record.getValue(SVFUSION.CHAINLINKS))
                    .chainTerminated(record.getValue(SVFUSION.CHAINTERMINATED) == 1)
                    .domainsKept(record.getValue(SVFUSION.DOMAINSKEPT))
                    .domainsLost(record.getValue(SVFUSION.DOMAINSLOST))
                    .skippedExonsUp(record.getValue(SVFUSION.SKIPPEDEXONSUP))
                    .skippedExonsDown(record.getValue(SVFUSION.SKIPPEDEXONSDOWN))
                    .fusedExonUp(record.getValue(SVFUSION.FUSEDEXONUP))
                    .fusedExonDown(record.getValue(SVFUSION.FUSEDEXONDOWN))
                    .geneStart("")
                    .geneContextStart("")
                    .geneTranscriptStart("")
                    .geneEnd("")
                    .geneContextEnd("")
                    .geneTranscriptEnd("")
                    .junctionCopyNumber(0.0)
                    .build();

            fusionList.add(fusion);
        }

        return fusionList;
    }

    @NotNull
    List<LinxBreakend> readBreakends(final String sample)
    {
        List<LinxBreakend> breakendList = Lists.newArrayList();

        Result<Record> result = context.select().from(SVBREAKEND).where(SVBREAKEND.SAMPLEID.eq(sample)).fetch();

        for(Record record : result)
        {
            LinxBreakend breakend = ImmutableLinxBreakend.builder()
                    .id(record.getValue(SVBREAKEND.ID).intValue())
                    .svId(record.getValue(SVBREAKEND.SVID))
                    .isStart(record.getValue(SVBREAKEND.STARTBREAKEND) == 1)
                    .gene(record.getValue(SVBREAKEND.GENE))
                    .transcriptId(record.getValue(SVBREAKEND.TRANSCRIPTID))
                    .canonical(record.getValue(SVBREAKEND.CANONICALTRANSCRIPT) == 1)
                    .geneOrientation(record.getValue(SVBREAKEND.GENEORIENTATION))
                    .disruptive(record.getValue(SVBREAKEND.DISRUPTIVE) == 1)
                    .reportedDisruption(record.getValue(SVBREAKEND.REPORTEDDISRUPTION) == 1)
                    .undisruptedCopyNumber(record.getValue(SVBREAKEND.UNDISRUPTEDCOPYNUMBER))
                    .regionType(TranscriptRegionType.valueOf(record.getValue(SVBREAKEND.REGIONTYPE)))
                    .codingType(TranscriptCodingType.valueOf(record.getValue(SVBREAKEND.CODINGCONTEXT)))
                    .biotype(record.getValue(SVBREAKEND.BIOTYPE))
                    .exonicBasePhase(record.getValue(SVBREAKEND.EXONICBASEPHASE))
                    .nextSpliceExonRank(record.getValue(SVBREAKEND.NEXTSPLICEEXONRANK) == null ?
                            0 : record.getValue(SVBREAKEND.NEXTSPLICEEXONRANK).intValue())
                    .nextSpliceExonPhase(record.getValue(SVBREAKEND.NEXTSPLICEEXONPHASE))
                    .nextSpliceDistance(record.getValue(SVBREAKEND.NEXTSPLICEDISTANCE))
                    .totalExonCount(record.getValue(SVBREAKEND.TOTALEXONCOUNT))
                    .type(StructuralVariantType.BND)
                    .chromosome("")
                    .orientation(0)
                    .strand(0)
                    .chrBand("")
                    .exonUp(-1)
                    .exonDown(-1)
                    .junctionCopyNumber(0)
                    .build();

            breakendList.add(breakend);
        }

        return breakendList;
    }
}

