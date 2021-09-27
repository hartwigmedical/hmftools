package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERGENEPANEL;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.InsertValuesStep12;

public class DriverGenePanelDAO {

    @NotNull
    private final DSLContext context;

    public DriverGenePanelDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeDriverGenes(@NotNull final List<DriverGene> driverGenes) {
        context.truncate(DRIVERGENEPANEL).execute();
        Timestamp timestamp = new Timestamp(new Date().getTime());
        InsertValuesStep12 inserter = context.insertInto(DRIVERGENEPANEL,
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
                DRIVERGENEPANEL.LIKELIHOODTYPE);

        for (DriverGene driverGene : driverGenes) {
            inserter.values(timestamp,
                    driverGene.gene(),
                    driverGene.reportMissenseAndInframe(),
                    driverGene.reportNonsenseAndFrameshift(),
                    driverGene.reportSplice(),
                    driverGene.reportDeletion(),
                    driverGene.reportDisruption(),
                    driverGene.reportAmplification(),
                    driverGene.reportSomaticHotspot(),
                    driverGene.reportGermlineVariant(),
                    driverGene.reportGermlineHotspot(),
                    driverGene.likelihoodType().toString());
        }

        inserter.execute();
    }
}
