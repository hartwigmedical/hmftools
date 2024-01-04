package com.hartwig.hmftools.esvee;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Ignore;
import org.junit.Test;

@SuppressWarnings("NewClassNamingConvention") // Named specifically to avoid common build tools finding it
@Ignore
public class IntegrationTestCreator
{
    public static final String COLO_SOMATIC_BAM = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003T.sv_prep_sorted.bam";
    public static final String COLO_GERMLINE_BAM = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003R.sv_prep_sorted.bam";
    public static final String COLO_JUNCTION_FILE = "/Users/james/code/data/COLO829/sv-prep2/COLO829v003T.sv_prep.junctions.tsv";

    public void createIntegrationTest(final String inputSomaticBAM, final String inputGermlineBAM, final String inputJunctionsFile, final String outputBaseName, final List<Pair<Integer, Integer>> locations)
    {
        final String outputBAM = outputBaseName + ".bam";
        final String outputJunctions = outputBaseName + ".tsv";


    }

    @Test
    public void createIntegrationTest()
    {
        createIntegrationTest(COLO_SOMATIC_BAM, COLO_GERMLINE_BAM, COLO_JUNCTION_FILE, "colo_3_6_6_3", List.of(

        ));
    }
}
