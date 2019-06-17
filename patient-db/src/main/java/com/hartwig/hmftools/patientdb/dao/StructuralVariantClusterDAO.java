package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVANNOTATION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVCLUSTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVLINK;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRALINSERTION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep17;
import org.jooq.InsertValuesStep22;
import org.jooq.InsertValuesStep5;
import org.jooq.InsertValuesStep9;

class StructuralVariantClusterDAO
{

    @NotNull
    private final DSLContext context;

    StructuralVariantClusterDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    public void writeClusters(final String sample, final List<LinxCluster> clusters)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).execute();

        for (List<LinxCluster> batch : Iterables.partition(clusters, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep9 inserter = context.insertInto(SVCLUSTER,
                    SVCLUSTER.SAMPLEID,
                    SVCLUSTER.MODIFIED,
                    SVCLUSTER.CLUSTERID,
                    SVCLUSTER.RESOLVEDTYPE,
                    SVCLUSTER.SYNTHETIC,
                    SVCLUSTER.SUBCLONAL,
                    SVCLUSTER.SUBTYPE,
                    SVCLUSTER.CLUSTERCOUNT,
                    SVCLUSTER.CLUSTERDESC);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep9 inserter, final String sample, final LinxCluster cluster)
    {
        //noinspection unchecked
        inserter.values(sample,
                timestamp,
                cluster.clusterId(),
                cluster.resolvedType(),
                cluster.synthetic(),
                cluster.subClonal(),
                cluster.subType(),
                cluster.clusterCount(),
                DatabaseUtil.checkStringLength(cluster.clusterDesc(), SVCLUSTER.CLUSTERDESC));
    }

    public void writeSvData(final String sample, final List<LinxSvData> svData)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).execute();

        for (List<LinxSvData> batch : Iterables.partition(svData, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep22 inserter = context.insertInto(SVANNOTATION,
                    SVANNOTATION.SAMPLEID,
                    SVANNOTATION.MODIFIED,
                    SVANNOTATION.SVID,
                    SVANNOTATION.CLUSTERID,
                    SVANNOTATION.CLUSTERREASON,
                    SVANNOTATION.FRAGILESITESTART,
                    SVANNOTATION.FRAGILESITEEND,
                    SVANNOTATION.ISFOLDBACK,
                    SVANNOTATION.LINETYPESTART,
                    SVANNOTATION.LINETYPEEND,
                    SVANNOTATION.PLOIDYMIN,
                    SVANNOTATION.PLOIDYMAX,
                    SVANNOTATION.GENESTART,
                    SVANNOTATION.GENEEND,
                    SVANNOTATION.REPLICATIONTIMINGSTART,
                    SVANNOTATION.REPLICATIONTIMINGEND,
                    SVANNOTATION.LOCALTOPOLOGYIDSTART,
                    SVANNOTATION.LOCALTOPOLOGYIDEND,
                    SVANNOTATION.LOCALTOPOLOGYSTART,
                    SVANNOTATION.LOCALTOPOLOGYEND,
                    SVANNOTATION.LOCALTICOUNTSTART,
                    SVANNOTATION.LOCALTICOUNTEND);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep22 inserter, final String sample, final LinxSvData svData)
    {
        //noinspection unchecked
        inserter.values(sample,
                timestamp,
                svData.svId(),
                svData.clusterId(),
                svData.clusterReason(),
                svData.fragileSiteStart(),
                svData.fragileSiteEnd(),
                svData.isFoldback(),
                svData.lineTypeStart(),
                svData.lineTypeEnd(),
                DatabaseUtil.decimal(svData.ploidyMin()),
                DatabaseUtil.decimal(svData.ploidyMax()),
                DatabaseUtil.checkStringLength(svData.geneStart(), SVANNOTATION.GENESTART),
                DatabaseUtil.checkStringLength(svData.geneEnd(), SVANNOTATION.GENEEND),
                DatabaseUtil.decimal(svData.replicationTimingStart()),
                DatabaseUtil.decimal(svData.replicationTimingEnd()),
                svData.localTopologyIdStart(),
                svData.localTopologyIdEnd(),
                svData.localTopologyStart(),
                svData.localTopologyEnd(),
                svData.localTICountStart(),
                svData.localTICountEnd());
    }

    public void writeLinks(final String sample, final List<LinxLink> links)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVLINK).where(SVLINK.SAMPLEID.eq(sample)).execute();

        for (List<LinxLink> batch : Iterables.partition(links, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep17 inserter = context.insertInto(SVLINK,
                    SVLINK.SAMPLEID,
                    SVLINK.MODIFIED,
                    SVLINK.CLUSTERID,
                    SVLINK.CHAINID,
                    SVLINK.CHAININDEX,
                    SVLINK.CHAINLINKCOUNT,
                    SVLINK.LOWERSVID,
                    SVLINK.UPPERSVID,
                    SVLINK.LOWERBREAKENDISSTART,
                    SVLINK.UPPERBREAKENDISSTART,
                    SVLINK.CHROMOSOME,
                    SVLINK.ARM,
                    SVLINK.ASSEMBLED,
                    SVLINK.TRAVERSEDSVCOUNT,
                    SVLINK.LINKLENGTH,
                    SVLINK.PLOIDY,
                    SVLINK.PSEUDOGENEINFO);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep17 inserter, final String sample, final LinxLink link)
    {
        //noinspection unchecked
        inserter.values(sample,
                timestamp,
                link.clusterId(),
                link.chainId(),
                link.chainIndex(),
                link.chainCount(),
                link.lowerSvId(),
                link.upperSvId(),
                link.lowerBreakendIsStart(),
                link.upperBreakendIsStart(),
                link.chromosome(),
                link.arm(),
                link.assembled(),
                link.traversedSVCount(),
                link.length(),
                DatabaseUtil.decimal(link.ploidy()),
                link.pseudogeneInfo());
    }

    public void writeViralInserts(final String sample, final List<LinxViralInsertFile> inserts)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(VIRALINSERTION).where(VIRALINSERTION.SAMPLEID.eq(sample)).execute();

        InsertValuesStep5 inserter = context.insertInto(VIRALINSERTION,
                VIRALINSERTION.SAMPLEID,
                VIRALINSERTION.MODIFIED,
                VIRALINSERTION.SVID,
                VIRALINSERTION.VIRUSID,
                VIRALINSERTION.VIRUSNAME);

        for (LinxViralInsertFile record : inserts)
        {
            addRecord(timestamp, inserter, sample, record);
        }

        inserter.execute();
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep5 inserter, final String sample, final LinxViralInsertFile insert)
    {
        //noinspection unchecked
        inserter.values(sample, timestamp, insert.SvId, insert.VirusId, insert.VirusName);
    }

    public void deleteClusterDataForSample(@NotNull String sample)
    {
        context.delete(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).execute();
        context.delete(SVLINK).where(SVLINK.SAMPLEID.eq(sample)).execute();
        context.delete(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).execute();
        context.delete(VIRALINSERTION).where(VIRALINSERTION.SAMPLEID.eq(sample)).execute();
    }
}
