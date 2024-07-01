package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.genotype.GenotypeStatus.UNKNOWN;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.checkStringLength;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GERMLINEVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.STRUCTURALVARIANTGERMLINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVBREAKENDGERMLINE;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.pathogenic.PathogenicSummary;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableGermlineVariantImpl;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.apache.logging.log4j.util.Strings;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;
import org.jooq.InsertValuesStepN;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.TableField;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineVariantDAO
{
    private final DSLContext context;

    public GermlineVariantDAO(final DSLContext context)
    {
        this.context = context;
    }

    public BufferedWriter<VariantContext> writer(String tumorSample, String referenceSample, String rnaSample)
    {
        BufferedWriterConsumer<VariantContext> consumer = new BufferedWriterConsumer<>()
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

    private void writeAll(final Timestamp timestamp, String tumorSample, String referenceSample, String rnaSample,
            final List<VariantContext> variants)
    {
        final InsertValuesStepN inserter = createInserter();
        variants.forEach(variant -> addRecord(timestamp, inserter, tumorSample, referenceSample, rnaSample, variant));
        inserter.execute();
    }

    public void deleteGermlineVariantsForSample(final String sampleId)
    {
        context.delete(GERMLINEVARIANT).where(GERMLINEVARIANT.SAMPLEID.eq(sampleId)).execute();
    }

    public void deleteGermlineStructuralVariantsForSample(final String sampleId)
    {
        context.delete(STRUCTURALVARIANTGERMLINE).where(STRUCTURALVARIANTGERMLINE.SAMPLEID.eq(sampleId)).execute();
    }

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
                referenceDepth.AlleleReadCount,
                referenceDepth.TotalReadCount,
                rnaDepth.AlleleReadCount,
                rnaDepth.TotalReadCount,
                tumorDepth.AlleleReadCount,
                tumorDepth.TotalReadCount,
                checkStringLength(decorator.localPhaseSetsToString(), GERMLINEVARIANT.LOCALPHASESET),
                decorator.adjustedVaf(),
                decorator.variantCopyNumber(),
                decorator.adjustedCopyNumber(),
                decorator.biallelic(),
                decorator.minorAlleleCopyNumber(),
                pathogenicSummary.ClinvarInfo,
                pathogenicSummary.Status.toString(),
                decorator.isPathogenic(),
                variantImpact.GeneName,
                variantImpact.GenesAffected,
                variantImpact.CanonicalEffect,
                variantImpact.CanonicalCodingEffect != CodingEffect.UNDEFINED ? variantImpact.CanonicalCodingEffect : Strings.EMPTY,
                checkTrimHgsvString(variantImpact.CanonicalHgvsCoding, GERMLINEVARIANT.CANONICALHGVSCODINGIMPACT),
                checkTrimHgsvString(variantImpact.CanonicalHgvsProtein, GERMLINEVARIANT.CANONICALHGVSPROTEINIMPACT),
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

    protected static String checkTrimHgsvString(final String hgvsStr, final TableField<?, String> field)
    {
        int maxLength = field.getDataType().length();

        if(hgvsStr.length() < maxLength)
            return hgvsStr;

        String trimmedStr = hgvsStr.substring(0, maxLength - 4);
        return trimmedStr + "...";
    }

    public void writeGermlineSVs(final String sample, final List<LinxGermlineSv> germlineSVs)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(STRUCTURALVARIANTGERMLINE).where(STRUCTURALVARIANTGERMLINE.SAMPLEID.eq(sample)).execute();

        InsertValuesStepN inserter = context.insertInto(STRUCTURALVARIANTGERMLINE,
                STRUCTURALVARIANTGERMLINE.SAMPLEID,
                STRUCTURALVARIANTGERMLINE.MODIFIED,
                STRUCTURALVARIANTGERMLINE.SVID,
                STRUCTURALVARIANTGERMLINE.VCFID,
                STRUCTURALVARIANTGERMLINE.CHROMOSOMESTART,
                STRUCTURALVARIANTGERMLINE.CHROMOSOMEEND,
                STRUCTURALVARIANTGERMLINE.POSITIONSTART,
                STRUCTURALVARIANTGERMLINE.POSITIONEND,
                STRUCTURALVARIANTGERMLINE.ORIENTATIONSTART,
                STRUCTURALVARIANTGERMLINE.ORIENTATIONEND,
                STRUCTURALVARIANTGERMLINE.TYPE,
                STRUCTURALVARIANTGERMLINE.FILTER,
                STRUCTURALVARIANTGERMLINE.EVENT,
                STRUCTURALVARIANTGERMLINE.QUALSCORE,
                STRUCTURALVARIANTGERMLINE.HOMOLOGYSEQUENCESTART,
                STRUCTURALVARIANTGERMLINE.HOMOLOGYSEQUENCEEND,
                STRUCTURALVARIANTGERMLINE.JUNCTIONCOPYNUMBER,
                STRUCTURALVARIANTGERMLINE.ADJUSTEDAFSTART,
                STRUCTURALVARIANTGERMLINE.ADJUSTEDAFEND,
                STRUCTURALVARIANTGERMLINE.ADJUSTEDCOPYNUMBERSTART,
                STRUCTURALVARIANTGERMLINE.ADJUSTEDCOPYNUMBEREND,
                STRUCTURALVARIANTGERMLINE.ADJUSTEDCOPYNUMBERCHANGESTART,
                STRUCTURALVARIANTGERMLINE.ADJUSTEDCOPYNUMBERCHANGEEND,
                STRUCTURALVARIANTGERMLINE.GERMLINEFRAGMENTS,
                STRUCTURALVARIANTGERMLINE.GERMLINEREFERENCEFRAGMENTSSTART,
                STRUCTURALVARIANTGERMLINE.GERMLINEREFERENCEFRAGMENTSEND,
                STRUCTURALVARIANTGERMLINE.TUMORFRAGMENTS,
                STRUCTURALVARIANTGERMLINE.TUMORREFERENCEFRAGMENTSSTART,
                STRUCTURALVARIANTGERMLINE.TUMORREFERENCEFRAGMENTSEND,
                STRUCTURALVARIANTGERMLINE.INSERTSEQUENCE,
                STRUCTURALVARIANTGERMLINE.INSERTSEQUENCEALIGNMENTS,
                STRUCTURALVARIANTGERMLINE.INSERTSEQUENCEREPEATCLASS,
                STRUCTURALVARIANTGERMLINE.INSERTSEQUENCEREPEATTYPE,
                STRUCTURALVARIANTGERMLINE.CLUSTERID,
                STRUCTURALVARIANTGERMLINE.CLUSTERCOUNT,
                STRUCTURALVARIANTGERMLINE.RESOLVEDTYPE,
                STRUCTURALVARIANTGERMLINE.LINKEDBYSTART,
                STRUCTURALVARIANTGERMLINE.LINKEDBYEND,
                STRUCTURALVARIANTGERMLINE.COHORTFREQUENCY);

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
                germlineSV.SvId,
                germlineSV.VcfId,
                germlineSV.ChromosomeStart,
                germlineSV.Type != SGL ? germlineSV.ChromosomeEnd : null,
                germlineSV.PositionStart,
                germlineSV.Type != SGL ? germlineSV.PositionEnd : null,
                germlineSV.OrientStart,
                germlineSV.Type != SGL ? germlineSV.OrientEnd : null,
                germlineSV.Type,
                checkStringLength(germlineSV.Filter, STRUCTURALVARIANTGERMLINE.FILTER),
                germlineSV.EventId,
                DatabaseUtil.decimal(germlineSV.QualScore),
                germlineSV.HomologyStart,
                germlineSV.HomologyEnd,
                DatabaseUtil.decimal(germlineSV.JunctionCopyNumber),
                DatabaseUtil.decimal(germlineSV.AdjustedAFStart),
                DatabaseUtil.decimal(germlineSV.AdjustedAFEnd),
                DatabaseUtil.decimal(germlineSV.AdjustedCopyNumberStart),
                DatabaseUtil.decimal(germlineSV.AdjustedCopyNumberEnd),
                DatabaseUtil.decimal(germlineSV.AdjustedCopyNumberChangeStart),
                DatabaseUtil.decimal(germlineSV.AdjustedCopyNumberChangeEnd),
                germlineSV.GermlineFragments,
                germlineSV.GermlineReferenceFragmentsStart,
                germlineSV.GermlineReferenceFragmentsEnd,
                germlineSV.TumorFragments,
                germlineSV.TumorReferenceFragmentsStart,
                germlineSV.TumorReferenceFragmentsEnd,
                checkStringLength(germlineSV.InsertSequence, STRUCTURALVARIANTGERMLINE.INSERTSEQUENCE),
                checkStringLength(germlineSV.InsertSequenceAlignments, STRUCTURALVARIANTGERMLINE.INSERTSEQUENCEALIGNMENTS),
                checkStringLength(germlineSV.InsertSequenceRepeatClass, STRUCTURALVARIANTGERMLINE.INSERTSEQUENCEREPEATCLASS),
                checkStringLength(germlineSV.InsertSequenceRepeatType, STRUCTURALVARIANTGERMLINE.INSERTSEQUENCEREPEATTYPE),
                germlineSV.ClusterId,
                germlineSV.ClusterCount,
                germlineSV.ResolvedType,
                germlineSV.LinkedByStart,
                germlineSV.LinkedByEnd,
                germlineSV.CohortFrequency);
    }

    public void writeGermlineBreakends(final String sample, final List<LinxBreakend> germlineBreakends)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVBREAKENDGERMLINE).where(SVBREAKENDGERMLINE.SAMPLEID.eq(sample)).execute();

        InsertValuesStep19 inserter = context.insertInto(SVBREAKENDGERMLINE,
                SVBREAKENDGERMLINE.SAMPLEID,
                SVBREAKENDGERMLINE.MODIFIED,
                SVBREAKENDGERMLINE.SVID,
                SVBREAKENDGERMLINE.STARTBREAKEND,
                SVBREAKENDGERMLINE.GENE,
                SVBREAKENDGERMLINE.TRANSCRIPTID,
                SVBREAKENDGERMLINE.CANONICALTRANSCRIPT,
                SVBREAKENDGERMLINE.GENEORIENTATION,
                SVBREAKENDGERMLINE.DISRUPTIVE,
                SVBREAKENDGERMLINE.REPORTEDDISRUPTION,
                SVBREAKENDGERMLINE.UNDISRUPTEDCOPYNUMBER,
                SVBREAKENDGERMLINE.REGIONTYPE,
                SVBREAKENDGERMLINE.CODINGTYPE,
                SVBREAKENDGERMLINE.BIOTYPE,
                SVBREAKENDGERMLINE.EXONICBASEPHASE,
                SVBREAKENDGERMLINE.NEXTSPLICEEXONRANK,
                SVBREAKENDGERMLINE.NEXTSPLICEEXONPHASE,
                SVBREAKENDGERMLINE.NEXTSPLICEDISTANCE,
                SVBREAKENDGERMLINE.TOTALEXONCOUNT);

        for(LinxBreakend germlineBreakend : germlineBreakends)
        {
            addRecord(timestamp, inserter, sample, germlineBreakend);
        }

        inserter.execute();
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep19 inserter, final String sample, final LinxBreakend germlineBreakend)
    {
        inserter.values(sample,
                timestamp,
                germlineBreakend.svId(),
                germlineBreakend.isStart(),
                germlineBreakend.gene(),
                germlineBreakend.transcriptId(),
                germlineBreakend.canonical(),
                germlineBreakend.geneOrientation(),
                germlineBreakend.disruptive(),
                germlineBreakend.reportedDisruption(),
                DatabaseUtil.decimal(germlineBreakend.undisruptedCopyNumber()),
                germlineBreakend.regionType(),
                germlineBreakend.codingType(),
                germlineBreakend.biotype(),
                germlineBreakend.exonicBasePhase(),
                germlineBreakend.nextSpliceExonRank(),
                germlineBreakend.nextSpliceExonPhase(),
                germlineBreakend.nextSpliceDistance(),
                germlineBreakend.totalExonCount());
    }

    public List<GermlineVariant> read(final String sample)
    {
        List<GermlineVariant> variants = Lists.newArrayList();

        Result<Record> result = context.select().from(GERMLINEVARIANT).where(GERMLINEVARIANT.SAMPLEID.eq(sample)).fetch();

        for(Record record : result)
        {
            variants.add(buildFromRecord(record));
        }

        return variants;
    }

    public static GermlineVariant buildFromRecord(final Record record)
    {
        Integer rnaAlleleReadCount = record.getValue(GERMLINEVARIANT.RNAALLELEREADCOUNT);
        Integer rnaTotalCount = record.getValue(GERMLINEVARIANT.RNATOTALREADCOUNT);

        AllelicDepth rnaAllelicDepth = rnaAlleleReadCount != null && rnaTotalCount != null ?
                new AllelicDepth(rnaTotalCount, rnaAlleleReadCount) : null;

        return ImmutableGermlineVariantImpl.builder()
                .chromosome(record.getValue(GERMLINEVARIANT.CHROMOSOME))
                .position(record.getValue(GERMLINEVARIANT.POSITION))
                .filter(record.getValue(GERMLINEVARIANT.FILTER))
                .type(VariantType.valueOf(record.getValue(GERMLINEVARIANT.TYPE)))
                .ref(record.getValue(GERMLINEVARIANT.REF))
                .alt(record.getValue(GERMLINEVARIANT.ALT))
                .gene(record.getValue(GERMLINEVARIANT.GENE))
                .genesAffected(record.getValue(GERMLINEVARIANT.GENESAFFECTED))
                .worstCodingEffect(record.getValue(GERMLINEVARIANT.WORSTCODINGEFFECT).isEmpty()
                        ? CodingEffect.UNDEFINED
                        : CodingEffect.valueOf(record.getValue(GERMLINEVARIANT.WORSTCODINGEFFECT)))
                .canonicalTranscript("")
                .canonicalEffect(record.getValue(GERMLINEVARIANT.CANONICALEFFECT))
                .canonicalCodingEffect(record.getValue(GERMLINEVARIANT.CANONICALCODINGEFFECT).isEmpty()
                        ? CodingEffect.UNDEFINED
                        : CodingEffect.valueOf(record.getValue(GERMLINEVARIANT.CANONICALCODINGEFFECT)))
                .canonicalHgvsCodingImpact(record.getValue(GERMLINEVARIANT.CANONICALHGVSCODINGIMPACT))
                .canonicalHgvsProteinImpact(record.getValue(GERMLINEVARIANT.CANONICALHGVSPROTEINIMPACT))
                .spliceRegion(DatabaseUtil.byteToBoolean(record.getValue(GERMLINEVARIANT.SPLICEREGION)))
                .otherReportedEffects(DatabaseUtil.valueNotNull(record.getValue(GERMLINEVARIANT.OTHERTRANSCRIPTEFFECTS)))
                .allelicDepth(new AllelicDepth(
                        record.getValue(GERMLINEVARIANT.GERMLINETOTALREADCOUNT), record.getValue(GERMLINEVARIANT.GERMLINEALLELEREADCOUNT)))
                .adjustedCopyNumber(record.getValue(GERMLINEVARIANT.COPYNUMBER))
                .adjustedVAF(record.getValue(GERMLINEVARIANT.ADJUSTEDVAF))
                .variantCopyNumber(record.getValue(GERMLINEVARIANT.VARIANTCOPYNUMBER))
                .biallelic(DatabaseUtil.byteToBoolean(record.getValue(GERMLINEVARIANT.BIALLELIC)))
                .reported(DatabaseUtil.byteToBoolean(record.getValue(GERMLINEVARIANT.REPORTED)))
                .trinucleotideContext(record.getValue(GERMLINEVARIANT.TRINUCLEOTIDECONTEXT))
                .microhomology(record.getValue(GERMLINEVARIANT.MICROHOMOLOGY))
                .repeatSequence(record.getValue(GERMLINEVARIANT.REPEATSEQUENCE))
                .repeatCount(record.getValue(GERMLINEVARIANT.REPEATCOUNT))
                .hotspot(Hotspot.valueOf(record.getValue(GERMLINEVARIANT.HOTSPOT)))
                .mappability(record.getValue(GERMLINEVARIANT.MAPPABILITY))
                .minorAlleleCopyNumber(record.getValue(GERMLINEVARIANT.MINORALLELECOPYNUMBER))
                .tier(VariantTier.fromString(record.get(GERMLINEVARIANT.TIER)))
                .rnaDepth(rnaAllelicDepth)
                .qual(record.get(GERMLINEVARIANT.QUAL))
                .genotypeStatus(UNKNOWN)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .clinvarInfo(record.get(GERMLINEVARIANT.CLINVARINFO))
                .pathogenicity(record.get(GERMLINEVARIANT.PATHOGENICITY))
                .pathogenic(record.get(GERMLINEVARIANT.PATHOGENIC).intValue() == 1)
                .build();
    }

}
