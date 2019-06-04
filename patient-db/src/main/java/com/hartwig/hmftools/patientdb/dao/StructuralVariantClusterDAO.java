package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLUSTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVLINK;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVLINXDATA;
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
import org.jooq.InsertValuesStep16;
import org.jooq.InsertValuesStep22;
import org.jooq.InsertValuesStep5;
import org.jooq.InsertValuesStep9;

public class StructuralVariantClusterDAO
{
    @NotNull
    private final DSLContext context;

    public StructuralVariantClusterDAO(@NotNull final DSLContext context)
    {
        this.context = context;
    }

    public void writeClusters(final String sample, final List<LinxCluster> clusters)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(CLUSTER).where(CLUSTER.SAMPLEID.eq(sample)).execute();

        InsertValuesStep9 inserter = context.insertInto(CLUSTER,
                CLUSTER.SAMPLEID,
                CLUSTER.MODIFIED,
                CLUSTER.CLUSTERID,
                CLUSTER.RESOLVEDTYPE,
                CLUSTER.SYNTHETIC,
                CLUSTER.SUBCLONAL,
                CLUSTER.SUBTYPE,
                CLUSTER.CLUSTERCOUNT,
                CLUSTER.CLUSTERDESC);

        for (List<LinxCluster> batch : Iterables.partition(clusters, DB_BATCH_INSERT_SIZE))
        {
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep9 inserter, final String sample, final LinxCluster cluster)
    {
        inserter.values(sample,
                timestamp,
                cluster.clusterId(),
                cluster.resolvedType(),
                cluster.synthetic(),
                cluster.subClonal(),
                cluster.subType(),
                cluster.clusterCount(),
                DatabaseUtil.checkStringLength(cluster.clusterDesc(), CLUSTER.CLUSTERDESC));
    }

    public void writeSvData(final String sample, final List<LinxSvData> svData)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVLINXDATA).where(SVLINXDATA.SAMPLEID.eq(sample)).execute();

        InsertValuesStep22 inserter = context.insertInto(SVLINXDATA,
                SVLINXDATA.SAMPLEID,
                SVLINXDATA.MODIFIED,
                SVLINXDATA.SVID,
                SVLINXDATA.CLUSTERID,
                SVLINXDATA.CLUSTERREASON,
                SVLINXDATA.FRAGILESITESTART,
                SVLINXDATA.FRAGILESITEEND,
                SVLINXDATA.ISFOLDBACK,
                SVLINXDATA.LINETYPESTART,
                SVLINXDATA.LINETYPEEND,
                SVLINXDATA.PLOIDYMIN,
                SVLINXDATA.PLOIDYMAX,
                SVLINXDATA.GENESTART,
                SVLINXDATA.GENEEND,
                SVLINXDATA.REPLICATIONTIMINGSTART,
                SVLINXDATA.REPLICATIONTIMINGEND,
                SVLINXDATA.LOCALTOPOLOGYIDSTART,
                SVLINXDATA.LOCALTOPOLOGYIDEND,
                SVLINXDATA.LOCALTOPOLOGYSTART,
                SVLINXDATA.LOCALTOPOLOGYEND,
                SVLINXDATA.LOCALTICOUNTSTART,
                SVLINXDATA.LOCALTICOUNTEND);


        for (List<LinxSvData> batch : Iterables.partition(svData, DB_BATCH_INSERT_SIZE))
        {
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep22 inserter, final String sample, final LinxSvData svData)
    {
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
                svData.ploidyMin(),
                svData.ploidyMax(),
                DatabaseUtil.checkStringLength(svData.geneStart(), SVLINXDATA.GENESTART),
                DatabaseUtil.checkStringLength(svData.geneEnd(), SVLINXDATA.GENEEND),
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

        InsertValuesStep16 inserter = context.insertInto(SVLINK,
                SVLINK.SAMPLEID,
                SVLINK.MODIFIED,
                SVLINK.CLUSTERID,
                SVLINK.CHAINID,
                SVLINK.CHAININDEX,
                SVLINK.CHAINLINKCOUNT,
                SVLINK.LOWERBREAKENDID,
                SVLINK.UPPERBREAKENDID,
                SVLINK.LOWERBREAKENDISSTART,
                SVLINK.UPPERBREAKENDISSTART,
                SVLINK.ARM,
                SVLINK.ASSEMBLED,
                SVLINK.TRAVERSEDSVCOUNT,
                SVLINK.LINKLENGTH,
                SVLINK.PLOIDY,
                SVLINK.PSEUDOGENEINFO);

        for (List<LinxLink> batch : Iterables.partition(links, DB_BATCH_INSERT_SIZE))
        {
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep16 inserter, final String sample, final LinxLink link)
    {
        inserter.values(sample,
                timestamp,
                link.clusterId(),
                link.chainId(),
                link.chainIndex(),
                link.chainCount(),
                link.lowerBreakendId(),
                link.upperBreakendId(),
                link.lowerBreakendIsStart(),
                link.upperBreakendIsStart(),
                link.arm(),
                link.assembled(),
                link.traversedSVCount(),
                link.length(),
                link.ploidy(),
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

        for (List<LinxViralInsertFile> batch : Iterables.partition(inserts, DB_BATCH_INSERT_SIZE))
        {
            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(Timestamp timestamp, InsertValuesStep5 inserter, final String sample, final LinxViralInsertFile insert)
    {
        inserter.values(sample,
                timestamp,
                insert.SvId,
                insert.VirusId,
                insert.VirusName);
    }

    public void deleteClusterDataForSample(@NotNull String sample)
    {
        context.delete(CLUSTER).where(CLUSTER.SAMPLEID.eq(sample)).execute();
        context.delete(SVLINK).where(SVLINK.SAMPLEID.eq(sample)).execute();
        context.delete(SVLINXDATA).where(SVLINXDATA.SAMPLEID.eq(sample)).execute();
        context.delete(VIRALINSERTION).where(VIRALINSERTION.SAMPLEID.eq(sample)).execute();
    }



}
