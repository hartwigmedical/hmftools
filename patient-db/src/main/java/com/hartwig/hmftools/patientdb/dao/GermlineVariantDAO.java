package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTGERMLINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVCLUSTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Germlinevariant.GERMLINEVARIANT;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep8;
import org.jooq.InsertValuesStepN;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineVariantDAO
{
    @NotNull
    private final DSLContext context;

    public GermlineVariantDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    @NotNull
    public BufferedWriter<VariantContext> writer(String tumorSample, String referenceSample, String rnaSample)
    {
        BufferedWriterConsumer<VariantContext> consumer = new BufferedWriterConsumer<VariantContext>()
        {
            @Override
            public void initialise()
            {
                deleteGermlineVariantsForSample(tumorSample);
            }

            @Override
            public void accept(final Timestamp timestamp, final List<VariantContext> entries)
            {
                writeAll(timestamp, tumorSample, referenceSample, rnaSample, entries);
            }
        };

        return new BufferedWriter<>(consumer);
    }

    private void writeAll(@NotNull final Timestamp timestamp, String tumorSample, String referenceSample, String rnaSample,
            @NotNull List<VariantContext> variants)
    {
        final InsertValuesStepN inserter = createInserter();
        variants.forEach(variant -> addRecord(timestamp, inserter, tumorSample, referenceSample, rnaSample, variant));
        inserter.execute();
    }

    void deleteGermlineVariantsForSample(@NotNull String sampleId)
    {
        context.delete(GERMLINEVARIANT).where(GERMLINEVARIANT.SAMPLEID.eq(sampleId)).execute();
        context.delete(STRUCTURALVARIANTGERMLINE).where(STRUCTURALVARIANTGERMLINE.SAMPLEID.eq(sampleId)).execute();
    }

    @NotNull
    private InsertValuesStepN createInserter()
    {
        return context.insertInto(GERMLINEVARIANT,
                GERMLINEVARIANT.MODIFIED,
                GERMLINEVARIANT.SAMPLEID,
                GERMLINEVARIANT.CHROMOSOME,
                GERMLINEVARIANT.POSITION,
                GERMLINEVARIANT.FILTER,
                GERMLINEVARIANT.TYPE,
                GERMLINEVARIANT.REF,
                GERMLINEVARIANT.ALT,
                GERMLINEVARIANT.QUAL,
                GERMLINEVARIANT.TIER,
                GERMLINEVARIANT.GERMLINEGENOTYPE,
                GERMLINEVARIANT.GERMLINEALLELEREADCOUNT,
                GERMLINEVARIANT.GERMLINETOTALREADCOUNT,
                GERMLINEVARIANT.RNAALLELEREADCOUNT,
                GERMLINEVARIANT.RNATOTALREADCOUNT,
                GERMLINEVARIANT.TUMORALLELEREADCOUNT,
                GERMLINEVARIANT.TUMORTOTALREADCOUNT,
                GERMLINEVARIANT.LOCALPHASESET,
                GERMLINEVARIANT.ADJUSTEDVAF,
                GERMLINEVARIANT.VARIANTCOPYNUMBER,
                GERMLINEVARIANT.COPYNUMBER,
                GERMLINEVARIANT.BIALLELIC,
                GERMLINEVARIANT.MINORALLELECOPYNUMBER,
                GERMLINEVARIANT.CLINVARINFO,
                GERMLINEVARIANT.PATHOGENICITY,
                GERMLINEVARIANT.PATHOGENIC,
                GERMLINEVARIANT.GENE,
                GERMLINEVARIANT.GENESAFFECTED,
                GERMLINEVARIANT.CANONICALEFFECT,
                GERMLINEVARIANT.CANONICALCODINGEFFECT,
                GERMLINEVARIANT.CANONICALHGVSCODINGIMPACT,
                GERMLINEVARIANT.CANONICALHGVSPROTEINIMPACT,
                GERMLINEVARIANT.SPLICEREGION,
                GERMLINEVARIANT.OTHERTRANSCRIPTEFFECTS,
                GERMLINEVARIANT.WORSTCODINGEFFECT,
                GERMLINEVARIANT.MICROHOMOLOGY,
                GERMLINEVARIANT.REPEATSEQUENCE,
                GERMLINEVARIANT.REPEATCOUNT,
                GERMLINEVARIANT.TRINUCLEOTIDECONTEXT,
                GERMLINEVARIANT.HOTSPOT,
                GERMLINEVARIANT.MAPPABILITY,
                GERMLINEVARIANT.REPORTED);
    }

    private static void addRecord(
            Timestamp timestamp, InsertValuesStepN inserter, String tumorSample, String referenceSample,
            String rnaSample, VariantContext variantContext)
    {
        final VariantContextDecorator decorator = new VariantContextDecorator(variantContext);
        final AllelicDepth tumorDepth = decorator.allelicDepth(tumorSample);
        final AllelicDepth referenceDepth = decorator.allelicDepth(referenceSample);
        final AllelicDepth rnaDepth = decorator.allelicDepth(rnaSample);
        final VariantImpact variantImpact = decorator.variantImpact();
        final PathogenicSummary pathogenicSummary = decorator.pathogenicSummary();

        inserter.values(timestamp,
                tumorSample,
                decorator.chromosome(),
                decorator.position(),
                decorator.filter(),
                decorator.type(),
                decorator.ref(),
                decorator.alt(),
                decorator.qual(),
                decorator.tier(),
                decorator.genotypeStatus(referenceSample).simplifiedDisplay(),
                referenceDepth.alleleReadCount(),
                referenceDepth.totalReadCount(),
                rnaDepth.alleleReadCount(),
                rnaDepth.totalReadCount(),
                tumorDepth.alleleReadCount(),
                tumorDepth.totalReadCount(),
                decorator.localPhaseSet(),
                decorator.adjustedVaf(),
                decorator.variantCopyNumber(),
                decorator.adjustedCopyNumber(),
                decorator.biallelic(),
                decorator.minorAlleleCopyNumber(),
                pathogenicSummary.clinvarInfo(),
                pathogenicSummary.pathogenicity().toString(),
                decorator.isPathogenic(),
                variantImpact.gene(),
                variantImpact.GenesAffected,
                variantImpact.CanonicalEffect,
                variantImpact.CanonicalCodingEffect != CodingEffect.UNDEFINED ? variantImpact.CanonicalCodingEffect : Strings.EMPTY,
                variantImpact.CanonicalHgvsCoding,
                variantImpact.CanonicalHgvsProtein,
                variantImpact.CanonicalSpliceRegion,
                variantImpact.OtherReportableEffects,
                variantImpact.WorstCodingEffect != CodingEffect.UNDEFINED ? variantImpact.WorstCodingEffect : Strings.EMPTY,
                decorator.microhomology(),
                decorator.repeatSequence(),
                decorator.repeatCount(),
                decorator.trinucleotideContext(),
                decorator.hotspot().toString(),
                decorator.mappability(),
                decorator.reported()
        );
    }

    public void writeGermlineSVs(final String sample, final List<LinxGermlineSv> germlineSVs)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(STRUCTURALVARIANTGERMLINE).where(STRUCTURALVARIANTGERMLINE.SAMPLEID.eq(sample)).execute();

        InsertValuesStepN inserter = context.insertInto(STRUCTURALVARIANTGERMLINE,
                STRUCTURALVARIANTGERMLINE.SAMPLEID,
                STRUCTURALVARIANTGERMLINE.MODIFIED,
                STRUCTURALVARIANTGERMLINE.CHROMOSOMESTART,
                STRUCTURALVARIANTGERMLINE.CHROMOSOMEEND,
                STRUCTURALVARIANTGERMLINE.POSITIONSTART,
                STRUCTURALVARIANTGERMLINE.POSITIONEND,
                STRUCTURALVARIANTGERMLINE.ORIENTATIONSTART,
                STRUCTURALVARIANTGERMLINE.ORIENTATIONEND,
                STRUCTURALVARIANTGERMLINE.GENE,
                STRUCTURALVARIANTGERMLINE.TYPE,
                STRUCTURALVARIANTGERMLINE.FILTER,
                STRUCTURALVARIANTGERMLINE.EVENT,
                STRUCTURALVARIANTGERMLINE.QUALSCORE,
                STRUCTURALVARIANTGERMLINE.GERMLINEFRAGMENTS,
                STRUCTURALVARIANTGERMLINE.GERMLINEREFERENCEFRAGMENTSSTART,
                STRUCTURALVARIANTGERMLINE.GERMLINEREFERENCEFRAGMENTSEND,
                STRUCTURALVARIANTGERMLINE.TUMORFRAGMENTS,
                STRUCTURALVARIANTGERMLINE.TUMORREFERENCEFRAGMENTSSTART,
                STRUCTURALVARIANTGERMLINE.TUMORREFERENCEFRAGMENTSEND,
                STRUCTURALVARIANTGERMLINE.INSERTSEQUENCE,
                STRUCTURALVARIANTGERMLINE.CLUSTERID,
                STRUCTURALVARIANTGERMLINE.CLUSTERCOUNT,
                STRUCTURALVARIANTGERMLINE.RESOLVEDTYPE,
                STRUCTURALVARIANTGERMLINE.LINKEDBYSTART,
                STRUCTURALVARIANTGERMLINE.LINKEDBYEND,
                STRUCTURALVARIANTGERMLINE.COHORTFREQUENCY,
                STRUCTURALVARIANTGERMLINE.REPORTED);

        for(LinxGermlineSv germlineSV : germlineSVs)
        {
            addRecord(timestamp, inserter, sample, germlineSV);
        }

        inserter.execute();
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStepN inserter, final String sample, final LinxGermlineSv germlineSV)
    {
        inserter.values(sample,
                timestamp,
                germlineSV.ChromosomeStart,
                germlineSV.Type != SGL ? germlineSV.ChromosomeEnd : null,
                germlineSV.PositionStart,
                germlineSV.Type != SGL ? germlineSV.PositionEnd : null,
                germlineSV.OrientStart,
                germlineSV.Type != SGL ? germlineSV.OrientEnd : null,
                germlineSV.GeneName,
                germlineSV.Type,
                germlineSV.Filter,
                germlineSV.EventId,
                DatabaseUtil.decimal(germlineSV.QualScore),
                germlineSV.GermlineFragments,
                germlineSV.GermlineReferenceFragmentsStart,
                germlineSV.GermlineReferenceFragmentsEnd,
                germlineSV.TumorFragments,
                germlineSV.TumorReferenceFragmentsStart,
                germlineSV.TumorReferenceFragmentsEnd,
                DatabaseUtil.checkStringLength(germlineSV.InsertSequence, STRUCTURALVARIANTGERMLINE.INSERTSEQUENCE),
                germlineSV.ClusterId,
                germlineSV.ClusterCount,
                germlineSV.ResolvedType,
                germlineSV.LinkedByStart,
                germlineSV.LinkedByEnd,
                germlineSV.CohortFrequency,
                germlineSV.Reported);
    }

}
