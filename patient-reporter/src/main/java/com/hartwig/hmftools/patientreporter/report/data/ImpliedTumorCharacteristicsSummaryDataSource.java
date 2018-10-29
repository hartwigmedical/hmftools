package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.report.components.ChordSection;
import com.hartwig.hmftools.patientreporter.report.components.MicrosatelliteSection;
import com.hartwig.hmftools.patientreporter.report.components.MutationalBurdenSection;
import com.hartwig.hmftools.patientreporter.report.components.MutationalLoadSection;
import com.hartwig.hmftools.patientreporter.report.pages.EvidencePage;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class ImpliedTumorCharacteristicsSummaryDataSource {

    public static final FieldBuilder<?> TUMOR_PURITY = field("Tumor Purity", String.class);
    public static final FieldBuilder<?> AVERAGE_TUMOR_PLOIDY = field("Average tumor ploidy", String.class);
    public static final FieldBuilder<?> MUTATIONAL_LOAD = field("Mutational Load", String.class);
    public static final FieldBuilder<?> MUTATIONAL_BURDEN = field("Mutational Burden", String.class);
    public static final FieldBuilder<?> MICROSATELLITE_INSTALBILITY = field("Microsatellite (in)stability", String.class);
    public static final FieldBuilder<?> HR_DEFICIENCY = field("HR_Deficiency", String.class);

    private ImpliedTumorCharacteristicsSummaryDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] ImpliedTumorCharacteristicsSummaryFields() {
        return new FieldBuilder<?>[] { TUMOR_PURITY, AVERAGE_TUMOR_PLOIDY, MUTATIONAL_LOAD, MUTATIONAL_BURDEN, MICROSATELLITE_INSTALBILITY, HR_DEFICIENCY };
    }

    @NotNull
    public static JRDataSource fromImpliedTumorCharacteristicsSummary(@NotNull AnalysedPatientReport report) {
        final DRDataSource impliedTumorCharacteristicsSummaryDataSource = new DRDataSource(TUMOR_PURITY.getName(), AVERAGE_TUMOR_PLOIDY.getName(),
                MUTATIONAL_LOAD.getName(),
                MUTATIONAL_BURDEN.getName(),
                MICROSATELLITE_INSTALBILITY.getName(),
                HR_DEFICIENCY.getName());
        impliedTumorCharacteristicsSummaryDataSource.add(EvidencePage.impliedPurityString(report),
                "tumor ploidy",
                MutationalLoadSection.interpret(report.tumorMutationalLoad(), report.fitStatus()),
                MutationalBurdenSection.interpret(report.tumorMutationalBurden(), report.fitStatus()),
                MicrosatelliteSection.interpret(report.microsatelliteIndelsPerMb(), report.fitStatus()),
                ChordSection.interpret(report.chordAnalysis().hrdValue(), report.fitStatus()));

        return impliedTumorCharacteristicsSummaryDataSource;
    }
}