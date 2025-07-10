package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TARGET;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENERAL_GC_TOLERANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_QUALITY_REJECT;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering a list of arbitrary regions provided by the user.
public class CustomRegions
{
    private static final ProbeSourceType PROBE_SOURCE = ProbeSourceType.CUSTOM;

    private static final ProbeEvalCriteria PROBE_EVAL_CRITERIA =
            new ProbeEvalCriteria(PROBE_QUALITY_REJECT, GENERAL_GC_TARGET, GENERAL_GC_TOLERANCE);

    private static final String FLD_EXTRA_INFO = "ExtraInfo";

    private static final Logger LOGGER = LogManager.getLogger(CustomRegions.class);

    public static ProbeGenerationResult generateProbes(final String customRegionFile, final ProbeGenerator probeGenerator)
    {
        LOGGER.info("Generating custom region probes");

        List<CustomRegion> regions = loadCustomRegionsFile(customRegionFile);

        ProbeGenerationResult result = regions.stream()
                .map(region -> generateProbes(region, probeGenerator))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);

        LOGGER.info("Done generating custom region probes");
        return result;
    }

    private record CustomRegion(
            ChrBaseRegion region,
            // Arbitrary descriptor for the user.
            String extraInfo
    )
    {
    }

    private static List<CustomRegion> loadCustomRegionsFile(final String filePath)
    {
        try(DelimFileReader reader = new DelimFileReader(filePath))
        {
            int chromosomeIdx = reader.getColumnIndex(FLD_CHROMOSOME);
            int posStartIdx = reader.getColumnIndex(FLD_POSITION_START);
            int posEndIdx = reader.getColumnIndex(FLD_POSITION_END);
            int extraInfoIdx = reader.getColumnIndex(FLD_EXTRA_INFO);

            List<CustomRegion> regions = reader.stream().map(row ->
            {
                String chromosome = row.get(chromosomeIdx);
                int start = row.getInt(posStartIdx);
                int end = row.getInt(posEndIdx);
                String extraInfo = row.get(extraInfoIdx);
                ChrBaseRegion baseRegion = new ChrBaseRegion(chromosome, start, end);
                return new CustomRegion(baseRegion, extraInfo);
            }).toList();

            LOGGER.info("Loaded {} custom regions from {}", regions.size(), filePath);
            return regions;
        }
    }

    private static ProbeGenerationResult generateProbes(final CustomRegion region, final ProbeGenerator probeGenerator)
    {
        ProbeSourceInfo source = new ProbeSourceInfo(PROBE_SOURCE, region.extraInfo());
        return probeGenerator.coverWholeRegion(region.region(), source, PROBE_EVAL_CRITERIA);
    }
}
