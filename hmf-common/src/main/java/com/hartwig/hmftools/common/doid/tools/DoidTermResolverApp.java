package com.hartwig.hmftools.common.doid.tools;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidNode;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DoidTermResolverApp
{
    private static final Logger LOGGER = LogManager.getLogger(DoidTermResolverApp.class);

    private static final String DOID_JSON = "doid_json";
    private static final String DOID_TO_RESOLVE = "doid_to_resolve";

    public static void main(String[] args) throws ParseException, IOException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        run(cmd.getOptionValue(DOID_JSON), cmd.getOptionValue(DOID_TO_RESOLVE));
    }

    private static void run(@NotNull String doidJsonPath, @NotNull String doid) throws IOException
    {
        DoidEntry entry = DiseaseOntology.readDoidOwlEntryFromDoidJson(doidJsonPath);
        DoidNode node = findDoidNode(entry.nodes(), doid);
        if(node != null)
        {
            LOGGER.info("Term for doid {} is '{}'", doid, node.doidTerm());
        }
        else
        {
            LOGGER.info("Could not resolve doid node with id '{}'", doid);
        }
    }

    @Nullable
    private static DoidNode findDoidNode(@NotNull List<DoidNode> nodes, @NotNull String doid)
    {
        for(DoidNode node : nodes)
        {
            if(node.doid().equals(doid))
            {
                return node;
            }
        }

        return null;
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();

        options.addOption(DOID_JSON, true, "Path towards to the json file of the doid ID of primary tumors.");
        options.addOption(DOID_TO_RESOLVE, true, "The DOID for which to resolve the term");

        return options;
    }
}
