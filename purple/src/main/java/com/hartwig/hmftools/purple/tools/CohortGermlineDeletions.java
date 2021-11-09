package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.drivers.GermlineDeletionFrequency.COHORT_DEL_FREQ_FILE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.purple.drivers.GermlineDeletionFrequency;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CohortGermlineDeletions
{
    private final GermlineDeletionFrequency mCohortFrequencies;
    private final String mCohortDeletionsFile;
    private final List<DriverGene> mDriverGenes;
    private final String mOutputDir;
    private final BufferedWriter mWriter;

    public static final String COHORT_DEL_FILE = "cohort_germline_del_file";

    public CohortGermlineDeletions(final CommandLine cmd)
    {
        mCohortDeletionsFile = cmd.getOptionValue(COHORT_DEL_FILE);
        mOutputDir = parseOutputDir(cmd);

        mCohortFrequencies = new GermlineDeletionFrequency(cmd.getOptionValue(COHORT_DEL_FREQ_FILE));

        mDriverGenes = loadDriverGenes(cmd);
        mWriter = initialiseWriter(mOutputDir);
    }

    public void run()
    {
        if(mCohortDeletionsFile == null)
        {
            PPL_LOGGER.error("missing germline deletion file");
            return;
        }

        processCohortFile(mCohortDeletionsFile);
    }

    private void processCohortFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            int sampleIdIndex = fieldsIndexMap.get("SampleId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int regionStartIndex = fieldsIndexMap.get("RegionStart");
            int regionEndIndex = fieldsIndexMap.get("RegionEnd");
            int germlineStatusIndex = fieldsIndexMap.get("GermlineStatus");
            int regionBafCountIndex = fieldsIndexMap.get("RegionBafCount");
            int regionDWCIndex = fieldsIndexMap.get("RegionDWC");
            int regionObsTumorRatioIndex = fieldsIndexMap.get("RegionObsTumorRatio");
            int regionObsNormalRatioIndex = fieldsIndexMap.get("RegionObsNormalRatio");
            int regionRefNormalisedCNIndex = fieldsIndexMap.get("RegionRefNormalisedCN");
            int purpleCNIndex = fieldsIndexMap.get("PurpleCN");
            int purpleDWCIndex = fieldsIndexMap.get("PurpleDWC");
            int gcContentIndex = fieldsIndexMap.get("GcContent");
            int majorAlleleCNIndex = fieldsIndexMap.get("MajorAlleleCN");
            int minorAlleleCNIndex = fieldsIndexMap.get("MinorAlleleCN");
            int geneNameIndex = fieldsIndexMap.get("GeneName");

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(",", -1);

                GermlineDeletion germlineDeletion = new GermlineDeletion(
                        values[sampleIdIndex], values[chromosomeIndex], Integer.parseInt(values[regionStartIndex]),
                        Integer.parseInt(values[regionEndIndex]), GermlineStatus.valueOf(values[germlineStatusIndex]),
                        Integer.parseInt(values[regionBafCountIndex]), Integer.parseInt(values[regionDWCIndex]),
                        Double.parseDouble(values[regionObsTumorRatioIndex]), Double.parseDouble(values[regionObsNormalRatioIndex]),
                        Double.parseDouble(values[regionRefNormalisedCNIndex]),
                        Double.parseDouble(values[purpleCNIndex]), Integer.parseInt(values[purpleDWCIndex]),
                        Double.parseDouble(values[gcContentIndex]), Double.parseDouble(values[majorAlleleCNIndex]),
                        Double.parseDouble(values[minorAlleleCNIndex]), values[geneNameIndex]);
            }
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to read cohort germline deletions file: {}", e.toString());
        }
    }

    private void processGermlineDeletion(final GermlineDeletion germlineDeletion)
    {

    }

    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String fileName = outputDir + "purple_cohort_germline_deletions.csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId,GeneName,Reportable,Chromosome,RegionStart,RegionEnd");
            writer.write(",RegionObsTumorRatio,RegionObsNormalRatio,RegionRefNormalisedCN");
            writer.write(",PurpleCN,PurpleDWC,GcContent,MajorAlleleCN,MinorAlleleCN");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    private List<DriverGene> loadDriverGenes(final CommandLine cmd)
    {
        if(DriverGenePanelConfig.isConfigured(cmd))
        {
            try
            {
                return DriverGenePanelConfig.driverGenes(cmd);
            }
            catch (IOException e)
            {
                PPL_LOGGER.error("invalid driver gene panel file({})", cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));
            }
        }

        return Lists.newArrayList();
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        options.addOption(COHORT_DEL_FREQ_FILE, true, "Germline cohort deletions frequency file");
        options.addOption(COHORT_DEL_FILE, true, "Germline cohort deletions file");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        CohortGermlineDeletions cohortGermlineDeletions = new CohortGermlineDeletions(cmd);
        cohortGermlineDeletions.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
