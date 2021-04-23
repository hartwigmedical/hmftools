package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RecurrentVariantFinder
{
    private final CohortConfig mConfig;
    private final String mSomaticVariantsFile;
    private BufferedWriter mWriter;

    protected static final String SOMATIC_VARIANT_FILE = "somatic_variants_file";

    public RecurrentVariantFinder(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mSomaticVariantsFile = cmd.getOptionValue(SOMATIC_VARIANT_FILE);
        mWriter = null;
        initialiseWriter();
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(SOMATIC_VARIANT_FILE, true, "Write out recurrent somatic variants for subsequent non-DB loading");
    }

    public void processCohortSomaticVariants()
    {
        findRecurrentVariants();
        closeBufferedWriter(mWriter);
    }

    private static final int LOG_NEXT = 100000;

    private void findRecurrentVariants()
    {
        if(!Files.exists(Paths.get(mSomaticVariantsFile)))
        {
            ISF_LOGGER.error("invalid somatic variant file({})", mSomaticVariantsFile);
            return;
        }

        // variants are order by chromosome and position, so can be analysed as a unit at a specific location
        // but they then need to be looked at for their specific base changes

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mSomaticVariantsFile));

            String line = fileReader.readLine();

            if (line == null)
                return;

            // mysql: [Warning] Using a password on the command line interface can be insecure.
            if(line.contains("mysql"))
            {
                line = fileReader.readLine();
            }

            final String delimiter = line.contains(DELIMITER) ? DELIMITER : "\t";
            final String[] fieldNames = line.split(delimiter, -1);
            if(fieldNames.length < 11)
            {
                ISF_LOGGER.error("invalid header: {}", line);
                return;
            }

            int sampleIdIndex = 0;
            int geneIndex = 1;
            int chrIndex = 2;
            int posIndex = 3;
            int typeIndex = 4;
            int refIndex = 5;
            int altIndex = 6;
            int effectIndex = 7;
            int proteinIndex = 8;
            int contextIndex = 9;
            int lpsIndex = 10;

            //sampleId	gene	chromosome	position	type	ref	alt	canonicalEffect	canonicalHgvsCodingImpact	trinucleotideContext	localPhaseSet
            if(!(fieldNames[sampleIdIndex].equalsIgnoreCase("sampleId")
            && fieldNames[chrIndex].equalsIgnoreCase("chromosome")
            && fieldNames[posIndex].equalsIgnoreCase("position")
            && fieldNames[refIndex].equalsIgnoreCase("ref") && fieldNames[altIndex].equalsIgnoreCase("alt")))
            {
                ISF_LOGGER.error("invalid header: {}", line);
                return;
            }

            List<SpliceVariant> chrPosVariants = Lists.newArrayList(); // for the current chromosome and position
            int processed = 0;
            int recurrent = 0;
            int nextLog = LOG_NEXT;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(delimiter);
                final String sampleId = items[sampleIdIndex];

                // must be recurrent in the specified cohort
                if(!mConfig.SampleData.SampleIds.contains(sampleId))
                    continue;

                final String chromosome = items[chrIndex];
                final int position = Integer.parseInt(items[posIndex]);
                final String ref = items[refIndex];
                final String alt = items[altIndex];

                if(chrPosVariants.isEmpty() || !(chrPosVariants.get(0).matchesLocation(chromosome, position)))
                {
                    for(final SpliceVariant variant : chrPosVariants)
                    {
                        if(variant.SampleIds.size() > 1)
                        {
                            writeRecurrentVariant(variant);
                            ++recurrent;
                        }
                    }

                    chrPosVariants.clear();
                }

                // add sample to any matching variant, otherwise register a new one at this location
                SpliceVariant currentVariant = chrPosVariants.stream()
                        .filter(x -> x.matches(chromosome, position, ref, alt)).findFirst().orElse(null);

                if(currentVariant == null)
                {
                    String localPhaseSet = items[lpsIndex];
                    currentVariant = new SpliceVariant(
                            items[geneIndex], chromosome, position, VariantType.valueOf(items[typeIndex]), ref, alt,
                            items[effectIndex], items[proteinIndex], items[contextIndex],
                            localPhaseSet.equals("NULL") ? -1 : Integer.parseInt(localPhaseSet));

                    chrPosVariants.add(currentVariant);
                }

                currentVariant.addSampleId(sampleId);

                ++processed;

                if(processed >= nextLog)
                {
                    nextLog += LOG_NEXT;
                    ISF_LOGGER.info("loaded {} somatic variants, recurrent({})", processed, recurrent);
                }
            }

            for(final SpliceVariant variant : chrPosVariants)
            {
                if(variant.SampleIds.size() > 1)
                {
                    writeRecurrentVariant(variant);
                    ++recurrent;
                }
            }

            ISF_LOGGER.info("found {} recurrent variants", recurrent);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load somatic variant data file: {}", e.toString());
            return;
        }
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("recurrent_variants.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write(String.format("SampleId,%s", SpliceVariant.header()));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write recurrent variants: {}", e.toString());
        }
    }

    private void writeRecurrentVariant(final SpliceVariant variant)
    {
        if(variant.SampleIds.size() < 2)
            return;

        ISF_LOGGER.debug("found recurrent somatic variant({}) samples({})", variant.key(), variant.SampleIds.size());

        try
        {
            for(String sampleId : variant.SampleIds)
            {
                mWriter.write(String.format("%s,%s", sampleId, variant.toCsv()));
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write recurrent variants: {}", e.toString());
        }
    }


}
