package com.hartwig.hmftools.lilac.cohort;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.ITEM_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SampleTruthSet
{
    private final Map<String, List<HlaAllele>> mSampleAlleles;

    private static final String TRUTH_SET_FILE = "truth_set_file";

    public SampleTruthSet(final CommandLine cmd)
    {
        mSampleAlleles = Maps.newHashMap();

        if(cmd != null)
            loadData(cmd.getOptionValue(TRUTH_SET_FILE));
    }

    public Map<String,List<HlaAllele>> getSampleAlleleSet() { return mSampleAlleles; }
    public List<HlaAllele> getSampleAlleles(final String sampleId) { return mSampleAlleles.get(sampleId); }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(TRUTH_SET_FILE, true, "File with sample truth-set alleles");
    }

    private void loadData(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            for(final String line : fileData)
            {
                if(line.startsWith("SampleId"))
                    continue;

                String[] items = line.split(",");

                if(items.length != 7)
                {
                    LL_LOGGER.error("invalid truth-set record: {}", line);
                    return;
                }

                String sampleId = items[0];

                List<HlaAllele> alleles = Lists.newArrayList();
                for(int index = 1; index < 7; ++index)
                {
                    alleles.add(HlaAllele.fromString(items[index]));
                }

                mSampleAlleles.put(sampleId, alleles);
            }

            LL_LOGGER.info("loaded {} sample truth-set allele from file({})", mSampleAlleles.size(), filename);
        }
        catch (IOException e)
        {
            LL_LOGGER.error("failed to read sample truth-set file({}): {}", filename, e.toString());
            return;
        }

    }
}
