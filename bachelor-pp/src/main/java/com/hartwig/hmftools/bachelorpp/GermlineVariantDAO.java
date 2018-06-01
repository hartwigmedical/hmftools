package com.hartwig.hmftools.bachelorpp;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStepN;

public class GermlineVariantDAO {

    @NotNull
    private final DSLContext context;

    private static final String PASS = "PASS";

    public GermlineVariantDAO (@NotNull final DSLContext context) {
        this.context = context;
    }

    public void write(final String sampleId, final List<BachelorGermlineVariant> bachRecords) {

        // first remove any existing records for this patient
        context.delete(GERMLINEVARIANT).where(GERMLINEVARIANT.SAMPLEID.eq(sampleId)).execute();

        final Timestamp timestamp = new Timestamp(new Date().getTime());

        InsertValuesStepN inserter = context.insertInto(GERMLINEVARIANT,
                GERMLINEVARIANT.SAMPLEID,
                GERMLINEVARIANT.CHROMOSOME,
                GERMLINEVARIANT.POSITION,
                GERMLINEVARIANT.FILTER,
                GERMLINEVARIANT.TYPE,
                GERMLINEVARIANT.REF,
                GERMLINEVARIANT.ALT,
                GERMLINEVARIANT.GENE,
                GERMLINEVARIANT.GENESEFFECTED,
                GERMLINEVARIANT.COSMICID,
                GERMLINEVARIANT.DBSNPID,
                GERMLINEVARIANT.WORSTEFFECT,
                GERMLINEVARIANT.WORSTCODINGEFFECT,
                GERMLINEVARIANT.WORSTEFFECTTRANSCRIPT,
                GERMLINEVARIANT.CANONICALEFFECT,
                GERMLINEVARIANT.CANONICALCODINGEFFECT,
                GERMLINEVARIANT.ALLELEREADCOUNT,
                GERMLINEVARIANT.TOTALREADCOUNT,
                GERMLINEVARIANT.ADJUSTEDCOPYNUMBER,
                GERMLINEVARIANT.ADJUSTEDVAF,
                GERMLINEVARIANT.HIGHCONFIDENCE,
                GERMLINEVARIANT.TRINUCLEOTIDECONTEXT,
                GERMLINEVARIANT.MICROHOMOLOGY,
                GERMLINEVARIANT.REPEATSEQUENCE,
                GERMLINEVARIANT.REPEATCOUNT,
                GERMLINEVARIANT.CLONALITY,
                GERMLINEVARIANT.BIALLELIC,
                GERMLINEVARIANT.HOTSPOT,
                GERMLINEVARIANT.MAPPABILITY,
                GERMLINEVARIANT.GERMLINESTATUS,
                GERMLINEVARIANT.MINORALLELEPLOIDY,
                GERMLINEVARIANT.PROGRAM,
                GERMLINEVARIANT.SOURCE,
                GERMLINEVARIANT.MODIFIED);

        for(final BachelorGermlineVariant bachRecord : bachRecords)
        {
            final EnrichedSomaticVariant region = bachRecord.getEnrichedVariant();
            inserter.values(
                    sampleId,
                    bachRecord.chromosome(),
                    bachRecord.position(),
                    PASS,
                    region.type(),
                    region.ref(),
                    region.alt(),
                    region.gene(),
                    region.genesEffected(),
                    region.cosmicID() == null ? "" : region.cosmicID(),
                    region.dbsnpID() == null ? "" : region.dbsnpID(),
                    region.worstEffect(),
                    region.worstCodingEffect(),
                    region.worstEffectTranscript(),
                    region.canonicalEffect(),
                    region.canonicalCodingEffect(),
                    bachRecord.getAltCount(),
                    bachRecord.getAltCount() + bachRecord.getRefCount(),
                    DatabaseUtil.decimal(region.adjustedCopyNumber()),
                    DatabaseUtil.decimal(region.adjustedVAF()),
                    region.highConfidenceRegion(),
                    region.trinucleotideContext(),
                    region.microhomology(),
                    region.repeatSequence(),
                    region.repeatCount(),
                    region.clonality(),
                    region.biallelic(),
                    region.hotspot(),
                    DatabaseUtil.decimal(region.mappability()),
                    region.germlineStatus(),
                    DatabaseUtil.decimal(region.minorAllelePloidy()),
                    bachRecord.program(),
                    bachRecord.source(),
                    timestamp);
        }

        inserter.execute();
    }
}
