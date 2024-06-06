package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVANNOTATION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVCLUSTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVDRIVER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVLINK;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep19;
import org.jooq.InsertValuesStep20;
import org.jooq.InsertValuesStep5;
import org.jooq.InsertValuesStep8;
import org.jooq.Record;
import org.jooq.Result;

class StructuralVariantClusterDAO
{
    private final DSLContext context;

    StructuralVariantClusterDAO(final DSLContext context) {
        this.context = context;
    }

    public void writeClusters(final String sample, final List<LinxCluster> clusters)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).execute();

        for (List<LinxCluster> batch : Iterables.partition(clusters, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep8 inserter = context.insertInto(SVCLUSTER,
                    SVCLUSTER.SAMPLEID,
                    SVCLUSTER.MODIFIED,
                    SVCLUSTER.CLUSTERID,
                    SVCLUSTER.CATEGORY,
                    SVCLUSTER.SYNTHETIC,
                    SVCLUSTER.RESOLVEDTYPE,
                    SVCLUSTER.CLUSTERCOUNT,
                    SVCLUSTER.CLUSTERDESC);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep8 inserter, final String sample, final LinxCluster cluster)
    {
        inserter.values(sample,
                timestamp,
                cluster.clusterId(),
                cluster.category(),
                cluster.synthetic(),
                cluster.resolvedType(),
                cluster.clusterCount(),
                DatabaseUtil.checkStringLength(cluster.clusterDesc(), SVCLUSTER.CLUSTERDESC));
    }

    public void writeSvData(final String sample, final List<LinxSvAnnotation> svData)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).execute();

        for (List<LinxSvAnnotation> batch : Iterables.partition(svData, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep20 inserter = context.insertInto(SVANNOTATION,
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
                    SVANNOTATION.JUNCTIONCOPYNUMBERMIN,
                    SVANNOTATION.JUNCTIONCOPYNUMBERMAX,
                    SVANNOTATION.GENESTART,
                    SVANNOTATION.GENEEND,
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

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep20 inserter, final String sample, final LinxSvAnnotation svData)
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
                DatabaseUtil.decimal(svData.junctionCopyNumberMin()),
                DatabaseUtil.decimal(svData.junctionCopyNumberMax()),
                DatabaseUtil.checkStringLength(svData.geneStart(), SVANNOTATION.GENESTART),
                DatabaseUtil.checkStringLength(svData.geneEnd(), SVANNOTATION.GENEEND),
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
            InsertValuesStep19 inserter = context.insertInto(SVLINK,
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
                    SVLINK.JUNCTIONCOPYNUMBER,
                    SVLINK.JUNCTIONCOPYNUMBERUNCERTAINTY,
                    SVLINK.PSEUDOGENEINFO,
                    SVLINK.ECDNA);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep19 inserter, final String sample, final LinxLink link)
    {
        inserter.values(sample,
                timestamp,
                link.clusterId(),
                link.chainId(),
                DatabaseUtil.checkStringLength(link.chainIndex(), SVLINK.CHAININDEX),
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
                link.junctionCopyNumber(),
                link.junctionCopyNumberUncertainty(),
                link.pseudogeneInfo(),
                link.ecDna());
    }

    public void writeDrivers(final String sample, final List<LinxDriver> drivers)
    {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVDRIVER).where(SVDRIVER.SAMPLEID.eq(sample)).execute();

        for (List<LinxDriver> batch : Iterables.partition(drivers, DB_BATCH_INSERT_SIZE))
        {
            InsertValuesStep5 inserter = context.insertInto(SVDRIVER,
                    SVDRIVER.SAMPLEID,
                    SVDRIVER.MODIFIED,
                    SVDRIVER.CLUSTERID,
                    SVDRIVER.GENE,
                    SVDRIVER.EVENTTYPE);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(
            final Timestamp timestamp, final InsertValuesStep5 inserter, final String sample, final LinxDriver driver)
    {
        inserter.values(sample,
                timestamp,
                driver.clusterId() == -1 ? null : driver.clusterId(),
                driver.gene(),
                DatabaseUtil.checkStringLength(driver.eventType().toString(), SVDRIVER.EVENTTYPE));
    }

    @NotNull
    public List<LinxSvAnnotation> readAnnotations(final String sample)
    {
        List<LinxSvAnnotation> svAnnotations = Lists.newArrayList();

        Result<Record> result = context.select().from(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).fetch();

        for (Record record : result)
        {
            LinxSvAnnotation svData = ImmutableLinxSvAnnotation.builder()
                    .svId(record.getValue(SVANNOTATION.SVID))
                    .vcfId("")
                    .clusterId(record.getValue(SVANNOTATION.CLUSTERID))
                    .clusterReason(record.getValue(SVANNOTATION.CLUSTERREASON))
                    .fragileSiteStart(record.getValue(SVANNOTATION.FRAGILESITESTART) == 1)
                    .fragileSiteEnd(record.getValue(SVANNOTATION.FRAGILESITEEND) == 1)
                    .isFoldback(record.getValue(SVANNOTATION.ISFOLDBACK) == 1)
                    .lineTypeStart(record.getValue(SVANNOTATION.LINETYPESTART))
                    .lineTypeEnd(record.getValue(SVANNOTATION.LINETYPEEND))
                    .junctionCopyNumberMin(record.getValue(SVANNOTATION.JUNCTIONCOPYNUMBERMIN))
                    .junctionCopyNumberMax(record.getValue(SVANNOTATION.JUNCTIONCOPYNUMBERMAX))
                    .geneStart(record.getValue(SVANNOTATION.GENESTART))
                    .geneEnd(record.getValue(SVANNOTATION.GENEEND))
                    .localTopologyIdStart(record.getValue(SVANNOTATION.LOCALTOPOLOGYIDSTART))
                    .localTopologyIdEnd(record.getValue(SVANNOTATION.LOCALTOPOLOGYIDEND))
                    .localTopologyStart(record.getValue(SVANNOTATION.LOCALTOPOLOGYSTART))
                    .localTopologyEnd(record.getValue(SVANNOTATION.LOCALTOPOLOGYEND))
                    .localTICountStart(record.getValue(SVANNOTATION.LOCALTICOUNTSTART))
                    .localTICountEnd(record.getValue(SVANNOTATION.LOCALTICOUNTEND))
                    .build();

            svAnnotations.add(svData);
        }

        return svAnnotations;
    }

    @NotNull
    public List<LinxCluster> readClusters(final String sample)
    {
        List<LinxCluster> clusterList = Lists.newArrayList();

        Result<Record> result = context.select().from(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).fetch();

        for (Record record : result)
        {
            LinxCluster cluster = ImmutableLinxCluster.builder()
                    .clusterId(record.getValue(SVCLUSTER.CLUSTERID))
                    .category(record.getValue(SVCLUSTER.CATEGORY))
                    .synthetic(record.getValue(SVCLUSTER.SYNTHETIC) == 1)
                    .resolvedType(record.getValue(SVCLUSTER.RESOLVEDTYPE))
                    .clusterCount(record.getValue(SVCLUSTER.CLUSTERCOUNT))
                    .clusterDesc(record.getValue(SVCLUSTER.CLUSTERDESC))
                    .build();

            clusterList.add(cluster);
        }

        return clusterList;
    }

    @NotNull
    public List<LinxDriver> readSvDrivers(final String sample)
    {
        List<LinxDriver> driverList = Lists.newArrayList();

        Result<Record> result = context.select().from(SVDRIVER).where(SVDRIVER.SAMPLEID.eq(sample)).fetch();

        for (Record record : result)
        {
            LinxDriver driver = ImmutableLinxDriver.builder()
                    .clusterId(DatabaseUtil.valueNotNull(record.getValue(SVDRIVER.CLUSTERID)))
                    .eventType(DriverEventType.valueOf(record.getValue(SVDRIVER.EVENTTYPE)))
                    .gene(record.getValue(SVDRIVER.GENE))
                    .build();

            driverList.add(driver);
        }

        return driverList;
    }

    public void deleteClusterDataForSample(final String sample)
    {
        context.delete(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).execute();
        context.delete(SVLINK).where(SVLINK.SAMPLEID.eq(sample)).execute();
        context.delete(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).execute();
        context.delete(SVDRIVER).where(SVDRIVER.SAMPLEID.eq(sample)).execute();
    }
}
