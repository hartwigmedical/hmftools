package com.hartwig.hmftools.common.virusinterpreter;

import static com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus.NO_ABNORMALITIES;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.viridian.ClonalityStatus;
import com.hartwig.hmftools.common.viridian.ReportQcStatus;
import com.hartwig.hmftools.common.viridian.VirusReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class VirusInterpreterDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(VirusInterpreterDataLoader.class);

    public static VirusInterpreterData load(final String basePath, final String sampleId) throws IOException
    {
        return VirusInterpreterDataLoader.load(AnnotatedVirusFile.generateFileName(basePath, sampleId));
    }

    public static VirusInterpreterData load(final String annotatedVirusTsv) throws IOException
    {
        List<AnnotatedVirus> viruses = AnnotatedVirusFile.read(annotatedVirusTsv);

        List<AnnotatedVirus> reportable = Lists.newArrayList();
        for(AnnotatedVirus virus : viruses)
        {
            if(virus.reported())
            {
                reportable.add(virus);
            }
        }

        LOGGER.debug(
                " loaded {} annotated viruses, reportable({}) from {}",
                viruses.size(), reportable.size(), annotatedVirusTsv);

        return ImmutableVirusInterpreterData.builder().allViruses(viruses).reportableViruses(reportable).build();
    }

    public static List<VirusReport> convertToViridianReports(final VirusInterpreterData virusInterpreterData)
    {
        List<VirusReport> virusReports = Lists.newArrayList();

        for(AnnotatedVirus annotatedVirus : virusInterpreterData.allViruses())
        {
            VirusReport virusReport = new VirusReport(
                    annotatedVirus.name(),
                    annotatedVirus.name(),
                    annotatedVirus.qcStatus() == NO_ABNORMALITIES ? ReportQcStatus.NO_ABNORMALITIES : ReportQcStatus.AMBIGUOUS_VIRUS,
                    annotatedVirus.interpretation(),
                    annotatedVirus.integrations(),
                    annotatedVirus.percentageCovered() * 0.01,
                    annotatedVirus.meanCoverage(),
                    null,
                    ClonalityStatus.NOT_PROFILED,
                    annotatedVirus.reported());

            virusReports.add(virusReport);
        }

        return virusReports;
    }
}
