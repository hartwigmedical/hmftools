package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.Config.DB_BATCH_INSERT_SIZE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVANNOTATION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVCLUSTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVDRIVER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVLINK;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRALINSERTION;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep18;
import org.jooq.InsertValuesStep22;
import org.jooq.InsertValuesStep5;
import org.jooq.InsertValuesStep9;

class StructuralVariantClusterDAO {

    @NotNull
    private final DSLContext context;

    StructuralVariantClusterDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeClusters(@NotNull String sample, @NotNull List<LinxCluster> clusters) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).execute();

        for (List<LinxCluster> batch : Iterables.partition(clusters, DB_BATCH_INSERT_SIZE)) {
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

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep9 inserter, @NotNull String sample,
            @NotNull LinxCluster cluster) {
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

    void writeSvData(@NotNull String sample, @NotNull List<LinxSvData> svData) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).execute();

        for (List<LinxSvData> batch : Iterables.partition(svData, DB_BATCH_INSERT_SIZE)) {
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

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep22 inserter, @NotNull String sample,
            @NotNull LinxSvData svData) {
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

    void writeLinks(@NotNull String sample, @NotNull List<LinxLink> links) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVLINK).where(SVLINK.SAMPLEID.eq(sample)).execute();

        for (List<LinxLink> batch : Iterables.partition(links, DB_BATCH_INSERT_SIZE)) {
            InsertValuesStep18 inserter = context.insertInto(SVLINK,
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
                    SVLINK.PLOIDYUNCERTAINTY,
                    SVLINK.PSEUDOGENEINFO);

            batch.forEach(entry -> addRecord(timestamp, inserter, sample, entry));
            inserter.execute();
        }
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep18 inserter, @NotNull String sample,
            @NotNull LinxLink link) {
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
                link.ploidy(),
                link.ploidyUncertainty(),
                link.pseudogeneInfo());
    }

    void writeViralInserts(@NotNull String sample, @NotNull List<LinxViralInsertFile> inserts) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(VIRALINSERTION).where(VIRALINSERTION.SAMPLEID.eq(sample)).execute();

        InsertValuesStep5 inserter = context.insertInto(VIRALINSERTION,
                VIRALINSERTION.SAMPLEID,
                VIRALINSERTION.MODIFIED,
                VIRALINSERTION.SVID,
                VIRALINSERTION.VIRUSID,
                VIRALINSERTION.VIRUSNAME);

        for (LinxViralInsertFile record : inserts) {
            addRecord(timestamp, inserter, sample, record);
        }

        inserter.execute();
    }

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep5 inserter, @NotNull String sample,
            @NotNull LinxViralInsertFile insert) {
        inserter.values(sample, timestamp, insert.SvId, insert.VirusId, insert.VirusName);
    }

    void writeDrivers(@NotNull String sample, @NotNull List<LinxDriver> drivers) {
        Timestamp timestamp = new Timestamp(new Date().getTime());

        context.delete(SVDRIVER).where(SVDRIVER.SAMPLEID.eq(sample)).execute();

        for (List<LinxDriver> batch : Iterables.partition(drivers, DB_BATCH_INSERT_SIZE)) {
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

    private static void addRecord(@NotNull Timestamp timestamp, @NotNull InsertValuesStep5 inserter, @NotNull String sample,
            @NotNull LinxDriver driver) {
        inserter.values(sample,
                timestamp,
                driver.clusterId() == -1 ? null : driver.clusterId(),
                driver.gene(),
                DatabaseUtil.checkStringLength(driver.eventType(), SVDRIVER.EVENTTYPE));
    }

    void deleteClusterDataForSample(@NotNull String sample) {
        context.delete(SVCLUSTER).where(SVCLUSTER.SAMPLEID.eq(sample)).execute();
        context.delete(SVLINK).where(SVLINK.SAMPLEID.eq(sample)).execute();
        context.delete(SVANNOTATION).where(SVANNOTATION.SAMPLEID.eq(sample)).execute();
        context.delete(SVDRIVER).where(SVDRIVER.SAMPLEID.eq(sample)).execute();
        context.delete(VIRALINSERTION).where(VIRALINSERTION.SAMPLEID.eq(sample)).execute();
    }
}
