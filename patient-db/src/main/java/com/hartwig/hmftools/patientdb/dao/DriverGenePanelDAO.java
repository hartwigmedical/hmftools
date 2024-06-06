package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERGENEPANEL;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep16;

public class DriverGenePanelDAO {

    private final DSLContext context;

    public DriverGenePanelDAO(final DSLContext context) {
        this.context = context;
    }

    void writeDriverGenes(final List<DriverGene> driverGenes)
    {
        context.truncate(DRIVERGENEPANEL).execute();
        Timestamp timestamp = new Timestamp(new Date().getTime());
        InsertValuesStep16 inserter = context.insertInto(DRIVERGENEPANEL,
                DRIVERGENEPANEL.MODIFIED,
                DRIVERGENEPANEL.GENE,
                DRIVERGENEPANEL.REPORTMISSENSE,
                DRIVERGENEPANEL.REPORTNONSENSE,
                DRIVERGENEPANEL.REPORTSPLICE,
                DRIVERGENEPANEL.REPORTDELETION,
                DRIVERGENEPANEL.REPORTDISRUPTION,
                DRIVERGENEPANEL.REPORTAMPLIFICATION,
                DRIVERGENEPANEL.REPORTSOMATICHOTSPOT,
                DRIVERGENEPANEL.REPORTGERMLINEVARIANT,
                DRIVERGENEPANEL.REPORTGERMLINEHOTSPOT,
                DRIVERGENEPANEL.LIKELIHOODTYPE,
                DRIVERGENEPANEL.REPORTGERMLINEDISRUPTION,
                DRIVERGENEPANEL.REPORTGERMLINEDELETION,
                DRIVERGENEPANEL.ADDITIONALREPORTEDTRANSCRIPTS,
                DRIVERGENEPANEL.REPORTPGX);

        for(DriverGene driverGene : driverGenes)
        {
            StringJoiner altTrans = new StringJoiner(";");
            driverGene.additionalReportedTranscripts().forEach(x -> altTrans.add(x));

            inserter.values(timestamp,
                    driverGene.gene(),
                    driverGene.reportMissenseAndInframe(),
                    driverGene.reportNonsenseAndFrameshift(),
                    driverGene.reportSplice(),
                    driverGene.reportDeletion(),
                    driverGene.reportDisruption(),
                    driverGene.reportAmplification(),
                    driverGene.reportSomaticHotspot(),
                    driverGene.reportGermlineVariant().toString(),
                    driverGene.reportGermlineHotspot().toString(),
                    driverGene.likelihoodType().toString(),
                    driverGene.reportGermlineDisruption().toString(),
                    driverGene.reportGermlineDeletion().toString(),
                    DatabaseUtil.checkStringLength(altTrans.toString(), DRIVERGENEPANEL.ADDITIONALREPORTEDTRANSCRIPTS),
                    driverGene.reportPGX());
        }

        inserter.execute();
    }
}
