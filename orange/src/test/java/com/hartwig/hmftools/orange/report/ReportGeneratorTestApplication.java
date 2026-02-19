package com.hartwig.hmftools.orange.report;

import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.CHORD_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.CUPPA_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.DRIVER_GENE_PANEL_TSV;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.ISOFOX_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.LILAC_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.LINX_GERMLINE_DATA_DIRECTORY;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.LINX_PLOT_DIRECTORY;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.LINX_SOMATIC_DATA_DIRECTORY;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.MELANOMA_DOID;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.PEACH_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.PIPELINE_VERSION_FILE;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.PURPLE_DATA_DIRECTORY;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.PURPLE_PLOT_DIRECTORY;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.REFERENCE_SAMPLE_BAM_METRICS_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.REFERENCE_SAMPLE_ID;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.REFERENCE_SAMPLE_REDUX_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.SIGNATURES_ETIOLOGY_TSV;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.SIGS_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.TUMOR_SAMPLE_BAM_METRICS_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.TUMOR_SAMPLE_ID;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.TUMOR_SAMPLE_REDUX_DIR;
import static com.hartwig.hmftools.orange.TestOrangeConfigFactory.VIRUS_DIR;

import java.io.File;
import java.nio.file.Files;
import java.time.LocalDate;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.OrangeAlgo;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.Nullable;

public class ReportGeneratorTestApplication
{
    private static final Logger LOGGER = LogManager.getLogger(ReportGeneratorTestApplication.class);

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";

    private static final boolean USE_MOCK_DATA_FOR_REPORT = false;
    private static final boolean REMOVE_UNREPORTED_VARIANTS = false;

    // Switch LIMIT_JSON_OUTPUT to true if you want to generate the real.orange.json test resource in orange-datamodel!
    private static final boolean LIMIT_JSON_OUTPUT = false;
    private static final Set<PurpleQCStatus> OVERRIDE_QC_STATUS = null;
    private static final boolean TUMOR_ONLY = false;

    public static void main(String[] args) throws Exception
    {
        Configurator.setRootLevel(Level.DEBUG);

        OrangeConfig config = buildConfig();

        ReportWriter writer = ReportWriterFactory.createToDiskWriter(config);

        if(!new File(REPORT_BASE_DIR).isDirectory())
        {
            LOGGER.warn("{} is not a directory. Can't write to disk", REPORT_BASE_DIR);
        }
        else
        {
            LOGGER.info("Deleting plot dir");
            deleteDir(new File(REPORT_BASE_DIR + File.separator + "plot"));
            writer.write(buildReport(config));
        }
    }

    private static OrangeRecord buildReport(final OrangeConfig config) throws Exception
    {
        if(USE_MOCK_DATA_FOR_REPORT)
        {
            LOGGER.info("Using mock data for report");
            return TestOrangeReportFactory.createProperTestReport();
        }

        OrangeRecord report = OrangeAlgo.fromConfig(config).run(config);
        if (OVERRIDE_QC_STATUS != null)
        {
            LOGGER.info("Overriding QC status to {}", OVERRIDE_QC_STATUS);
            report = overwritePurpleQCStatus(report, OVERRIDE_QC_STATUS);
        }

        OrangeRecord filtered;
        if(REMOVE_UNREPORTED_VARIANTS)
        {
            filtered = removeUnreported(report);
        }
        else
        {
            filtered = report;
        }

        OrangeRecord finalReport = ImmutableOrangeRecord.builder().from(filtered).sampleId("Test").build();

        return finalReport;
    }

    private static OrangeRecord overwritePurpleQCStatus(final OrangeRecord report, final Set<PurpleQCStatus> newStatus)
    {
        return ImmutableOrangeRecord.builder().from(report)
                .purple(ImmutablePurpleRecord.builder().from(report.purple())
                        .fit(ImmutablePurpleFit.builder()
                                .from(report.purple().fit())
                                .qc(ImmutablePurpleQC.builder()
                                        .from(report.purple().fit().qc())
                                        .status(newStatus)
                                        .build())
                                .build())
                        .build())
                .build();
    }

    private static OrangeConfig buildConfig()
    {
        return new OrangeConfig(
                ExperimentType.WHOLE_GENOME, TUMOR_SAMPLE_ID, TUMOR_ONLY ? null : REFERENCE_SAMPLE_ID, null,
                RefGenomeVersion.V37, Collections.emptySet(), LocalDate.now(),
                null, MELANOMA_DOID, SIGNATURES_ETIOLOGY_TSV, DRIVER_GENE_PANEL_TSV,
                PIPELINE_VERSION_FILE, PURPLE_DATA_DIRECTORY, PURPLE_PLOT_DIRECTORY, LINX_SOMATIC_DATA_DIRECTORY,
                LINX_GERMLINE_DATA_DIRECTORY, LINX_PLOT_DIRECTORY, TUMOR_SAMPLE_BAM_METRICS_DIR, REFERENCE_SAMPLE_BAM_METRICS_DIR,
                TUMOR_SAMPLE_REDUX_DIR, REFERENCE_SAMPLE_REDUX_DIR, LILAC_DIR, CHORD_DIR, CUPPA_DIR, PEACH_DIR, SIGS_DIR, VIRUS_DIR,
                ISOFOX_DIR, false, true);
    }

    private static OrangeRecord removeUnreported(final OrangeRecord report)
    {
        ImmutableOrangeRecord.Builder builder = ImmutableOrangeRecord.builder()
                .from(report)
                .purple(ImmutablePurpleRecord.builder()
                        .from(report.purple())
                        .somaticCopyNumbers(Lists.newArrayList())
                        .somaticGeneCopyNumbers(retainReportableCopyNumbers(report.purple().somaticGeneCopyNumbers(),
                                report.purple().somaticDrivers()))
                        .build())
                .linx(ImmutableLinxRecord.builder()
                        .from(report.linx())
                        .build());

        if(report.isofox() != null)
        {
            builder.isofox(ImmutableIsofoxRecord.builder()
                    .from(report.isofox())
                    .allGeneExpressions(Lists.newArrayList())
                    .allFusions(Lists.newArrayList())
                    .allNovelSpliceJunctions(Lists.newArrayList())
                    .build());
        }

        return builder.build();
    }

    private static List<PurpleGeneCopyNumber> retainReportableCopyNumbers(
            final List<PurpleGeneCopyNumber> geneCopyNumbers, final List<PurpleDriver> drivers)
    {
        List<String> copyNumberDriverGenes = Lists.newArrayList();
        for(PurpleDriver driver : drivers)
        {
            if(driver.type() == PurpleDriverType.AMP || driver.type() == PurpleDriverType.PARTIAL_AMP
                    || driver.type() == PurpleDriverType.DEL
                    || driver.type() == PurpleDriverType.GERMLINE_DELETION)
            {
                copyNumberDriverGenes.add(driver.gene());
            }
        }

        List<PurpleGeneCopyNumber> reportable = Lists.newArrayList();
        for(PurpleGeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(copyNumberDriverGenes.contains(geneCopyNumber.gene()))
            {
                reportable.add(geneCopyNumber);
            }
        }
        return reportable;
    }

    @Nullable
    private static List<LinxSvAnnotation> retainReportableStructuralVariants(@Nullable List<LinxSvAnnotation> structuralVariants,
            @Nullable List<LinxBreakend> reportableBreakends)
    {
        if(structuralVariants == null || reportableBreakends == null)
        {
            return null;
        }

        List<LinxSvAnnotation> reportable = Lists.newArrayList();
        for(LinxSvAnnotation structuralVariant : structuralVariants)
        {
            if(isReportableSv(structuralVariant, reportableBreakends))
            {
                reportable.add(structuralVariant);
            }
        }
        return reportable;
    }

    private static boolean isReportableSv(final LinxSvAnnotation structuralVariant, final List<LinxBreakend> reportableBreakends)
    {
        for(LinxBreakend breakend : reportableBreakends)
        {
            if(breakend.svId() == structuralVariant.svId())
            {
                return true;
            }
        }
        return false;
    }

    private static void deleteDir(final File file)
    {
        File[] contents = file.listFiles();
        if(contents != null)
        {
            for(File content : contents)
            {
                if(!Files.isSymbolicLink(content.toPath()))
                {
                    deleteDir(content);
                }
            }
        }
        file.delete();
    }
}
