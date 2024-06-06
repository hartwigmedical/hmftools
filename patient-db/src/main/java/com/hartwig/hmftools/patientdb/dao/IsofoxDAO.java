package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.rna.RnaQcFilter.qcFiltersToString;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.checkStringLength;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.GENEEXPRESSION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.NOVELSPLICEJUNCTION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.RNAFUSION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.RNASTATISTICS;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaStatistics;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep10;
import org.jooq.InsertValuesStep15;
import org.jooq.InsertValuesStep20;

public class IsofoxDAO
{
    private final DSLContext context;

    private static final int DB_BATCH_INSERT_SIZE = 10000;

    public IsofoxDAO(final DSLContext context) {
        this.context = context;
    }

    void deleteSampleData(final String sampleId)
    {
        context.delete(RNASTATISTICS).where(RNASTATISTICS.SAMPLEID.eq(sampleId)).execute();
        context.delete(GENEEXPRESSION).where(GENEEXPRESSION.SAMPLEID.eq(sampleId)).execute();
        context.delete(NOVELSPLICEJUNCTION).where(NOVELSPLICEJUNCTION.SAMPLEID.eq(sampleId)).execute();
        context.delete(RNAFUSION).where(RNAFUSION.SAMPLEID.eq(sampleId)).execute();
    }

    public void writeRnaStatistics(final String sampleId, final RnaStatistics statistics)
    {
        context.delete(RNASTATISTICS).where(RNASTATISTICS.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        InsertValuesStep15 inserter = context.insertInto(RNASTATISTICS,
                RNASTATISTICS.MODIFIED,
                RNASTATISTICS.SAMPLEID,
                RNASTATISTICS.QCSTATUS,
                RNASTATISTICS.READLENGTH,
                RNASTATISTICS.TOTALFRAGMENTS,
                RNASTATISTICS.DUPLICATES,
                RNASTATISTICS.SPLICEDPERCENT,
                RNASTATISTICS.UNSPLICEDPERCENT,
                RNASTATISTICS.ALTERNATESPLICEPERCENT,
                RNASTATISTICS.CHIMERICPERCENT,
                RNASTATISTICS.FRAGMENTLENGTHPCT05,
                RNASTATISTICS.FRAGMENTLENGTHPCT50,
                RNASTATISTICS.FRAGMENTLENGTHPCT95,
                RNASTATISTICS.ENRICHEDGENEPERCENT,
                RNASTATISTICS.MEDIANGCRATIO);

        inserter.values(
                timestamp,
                sampleId,
                DatabaseUtil.checkStringLength(qcFiltersToString(statistics.qcStatus()), RNASTATISTICS.QCSTATUS),
                statistics.readLength(),
                statistics.totalFragments(),
                statistics.duplicateFragments(),
                DatabaseUtil.decimal(statistics.splicedFragmentPerc()),
                DatabaseUtil.decimal(statistics.unsplicedFragmentPerc()),
                DatabaseUtil.decimal(statistics.altFragmentPerc()),
                DatabaseUtil.decimal(statistics.chimericFragmentPerc()),
                DatabaseUtil.decimal(statistics.fragmentLength5thPercent()),
                DatabaseUtil.decimal(statistics.fragmentLength50thPercent()),
                DatabaseUtil.decimal(statistics.fragmentLength95thPercent()),
                DatabaseUtil.decimal(statistics.enrichedGenePercent()),
                DatabaseUtil.decimal(statistics.medianGCRatio()));
        inserter.execute();
    }

    public void writeGeneExpressions(final String sampleId, final List<GeneExpression> geneExpressions)
    {
        context.delete(GENEEXPRESSION).where(GENEEXPRESSION.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<GeneExpression> batch : Iterables.partition(geneExpressions, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep10 inserter = context.insertInto(GENEEXPRESSION,
                    GENEEXPRESSION.MODIFIED,
                    GENEEXPRESSION.SAMPLEID,
                    GENEEXPRESSION.GENE,
                    GENEEXPRESSION.TPM,
                    GENEEXPRESSION.SPLICEDFRAGMENTS,
                    GENEEXPRESSION.UNSPLICEDFRAGMENTS,
                    GENEEXPRESSION.MEDIANTPMCANCER,
                    GENEEXPRESSION.PERCENTILECANCER,
                    GENEEXPRESSION.MEDIANTPMCOHORT,
                    GENEEXPRESSION.PERCENTILECOHORT);

            batch.forEach(entry -> addRecord(timestamp, inserter, sampleId, entry));
            inserter.execute();
        }
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep10 inserter, final String sampleId, final GeneExpression geneExpression)
    {
        inserter.values(
                timestamp,
                sampleId,
                checkStringLength(geneExpression.geneName(), GENEEXPRESSION.GENE),
                geneExpression.tpm(),
                geneExpression.splicedFragments(),
                geneExpression.unsplicedFragments(),
                geneExpression.medianTpmCancer(),
                DatabaseUtil.decimal(geneExpression.percentileCancer()),
                geneExpression.medianTpmCohort(),
                DatabaseUtil.decimal(geneExpression.percentileCohort()));
    }

    public void writeNovelSpliceJunctions(final String sampleId, final List<NovelSpliceJunction> novelJunctions)
    {
        context.delete(NOVELSPLICEJUNCTION).where(NOVELSPLICEJUNCTION.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<NovelSpliceJunction> batch : Iterables.partition(novelJunctions, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep15 inserter = context.insertInto(NOVELSPLICEJUNCTION,
                    NOVELSPLICEJUNCTION.MODIFIED,
                    NOVELSPLICEJUNCTION.SAMPLEID,
                    NOVELSPLICEJUNCTION.GENE,
                    NOVELSPLICEJUNCTION.CHROMOSOME,
                    NOVELSPLICEJUNCTION.JUNCTIONSTART,
                    NOVELSPLICEJUNCTION.JUNCTIONEND,
                    NOVELSPLICEJUNCTION.TYPE,
                    NOVELSPLICEJUNCTION.FRAGMENTCOUNT,
                    NOVELSPLICEJUNCTION.DEPTHSTART,
                    NOVELSPLICEJUNCTION.DEPTHEND,
                    NOVELSPLICEJUNCTION.REGIONSTART,
                    NOVELSPLICEJUNCTION.REGIONEND,
                    NOVELSPLICEJUNCTION.BASESSTART,
                    NOVELSPLICEJUNCTION.BASESEND,
                    NOVELSPLICEJUNCTION.COHORTFREQUENCY);

            batch.forEach(entry -> addRecord(timestamp, inserter, sampleId, entry));
            inserter.execute();
        }
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep15 inserter, final String sampleId, final NovelSpliceJunction novelJunction)
    {
        inserter.values(
                timestamp,
                sampleId,
                checkStringLength(novelJunction.geneName(), NOVELSPLICEJUNCTION.GENE),
                novelJunction.chromosome(),
                novelJunction.junctionStart(),
                novelJunction.junctionEnd(),
                novelJunction.type(),
                novelJunction.fragmentCount(),
                novelJunction.depthStart(),
                novelJunction.depthEnd(),
                novelJunction.regionStart(),
                novelJunction.regionEnd(),
                novelJunction.basesStart(),
                novelJunction.basesEnd(),
                novelJunction.cohortFrequency());
    }

    public void writeRnaFusions(final String sampleId, final List<RnaFusion> rnaFusions)
    {
        context.delete(RNAFUSION).where(RNAFUSION.SAMPLEID.eq(sampleId)).execute();

        Timestamp timestamp = new Timestamp(new Date().getTime());

        for (List<RnaFusion> batch : Iterables.partition(rnaFusions, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep20 inserter = context.insertInto(RNAFUSION,
                    RNAFUSION.MODIFIED,
                    RNAFUSION.SAMPLEID,
                    RNAFUSION.NAME,
                    RNAFUSION.CHROMOSOMEUP,
                    RNAFUSION.CHROMOSOMEDOWN,
                    RNAFUSION.POSITIONUP,
                    RNAFUSION.POSITIONDOWN,
                    RNAFUSION.ORIENTATIONUP,
                    RNAFUSION.ORIENTATIONDOWN,
                    RNAFUSION.JUNCTIONTYPEUP,
                    RNAFUSION.JUNCTIONTYPEDOWN,
                    RNAFUSION.SVTYPE,
                    RNAFUSION.SPLITFRAGMENTS,
                    RNAFUSION.REALIGNEDFRAGS,
                    RNAFUSION.DISCORDANTFRAGS,
                    RNAFUSION.DEPTHUP,
                    RNAFUSION.DEPTHDOWN,
                    RNAFUSION.MAXANCHORLENGTHUP,
                    RNAFUSION.MAXANCHORLENGTHDOWN,
                    RNAFUSION.COHORTFREQUENCY);

            batch.forEach(entry -> addRecord(timestamp, inserter, sampleId, entry));
            inserter.execute();
        }
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep20 inserter, final String sampleId, final RnaFusion rnaFusion)
    {
        inserter.values(
                timestamp,
                sampleId,
                checkStringLength(rnaFusion.name(), RNAFUSION.NAME),
                rnaFusion.chromosomeUp(),
                rnaFusion.chromosomeDown(),
                rnaFusion.positionUp(),
                rnaFusion.positionDown(),
                rnaFusion.orientationUp(),
                rnaFusion.orientationDown(),
                rnaFusion.junctionTypeUp(),
                rnaFusion.junctionTypeDown(),
                rnaFusion.svType(),
                rnaFusion.splitFragments(),
                rnaFusion.realignedFrags(),
                rnaFusion.discordantFrags(),
                rnaFusion.depthUp(),
                rnaFusion.depthDown(),
                rnaFusion.maxAnchorLengthUp(),
                rnaFusion.maxAnchorLengthDown(),
                rnaFusion.cohortFrequency());
    }
}

