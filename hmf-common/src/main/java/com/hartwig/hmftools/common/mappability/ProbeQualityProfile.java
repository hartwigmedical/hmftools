package com.hartwig.hmftools.common.mappability;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Reads the data written by GeneUtils ProbeQualityProfiler.
// See that class for more context.
public class ProbeQualityProfile
{
    public static final String PROBE_QUALITY_FILE_CONFIG = "probe_quality_profile.tsv";
    public static final String PROBE_QUALITY_FILE_DESC = "Genome regions to unmap reads";

    // Must match the config used for generating the file
    private static final int BASE_WINDOW_LENGTH = 40;
    private static final int BASE_WINDOW_SPACING = 20;

    private static final String QUALITY_SCORE_FIELD = "QualityScore";

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(PROBE_QUALITY_FILE_CONFIG, false, PROBE_QUALITY_FILE_DESC);
    }

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityProfile.class);

    public static Map<String, List<ProbeQualityWindow>> loadProbeQualityProfile(final String filePath)
    {
        LOGGER.debug("Loading probe quality profile file: {}", filePath);
        long startTimeMs = System.currentTimeMillis();
        Map<String, List<ProbeQualityWindow>> result = new HashMap<>();
        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            for(DelimFileReader.Row row : reader)
            {
                String chromosome = row.get(FLD_CHROMOSOME);
                int start = row.getInt(FLD_POSITION_START);
                double qualityScore = row.getDouble(QUALITY_SCORE_FIELD);
                int end = start + BASE_WINDOW_LENGTH;
                ProbeQualityWindow window = new ProbeQualityWindow(start, end, (float) qualityScore);
                result.computeIfAbsent(chromosome, k -> new ArrayList<>()).add(window);
            }
        }
        LOGGER.info("Loading complete, mins({})", runTimeMinsStr(startTimeMs));
        return result;
    }
}
