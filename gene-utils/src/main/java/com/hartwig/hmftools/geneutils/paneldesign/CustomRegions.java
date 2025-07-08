package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering a list of arbitrary regions provided by the user.
public class CustomRegions
{
    private static final String FLD_EXTRA_INFO = "ExtraInfo";

    private static final Logger LOGGER = LogManager.getLogger(CustomRegions.class);

    public static List<ProbeCandidate> createProbeCandidates(final String customRegionFile)
    {
        List<CustomRegion> regions = loadCustomRegionsFile(customRegionFile);
        List<ProbeCandidate> probesCandidates = regions.stream()
                .flatMap(region -> createProbeCandidates(region).stream()).toList();
        return probesCandidates;
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

            LOGGER.info("Loaded {} custom regions", regions.size());
            return regions;
        }
    }

    private static List<ProbeCandidate> createProbeCandidates(final CustomRegion region)
    {
        ProbeSourceInfo source = new ProbeSourceInfo(ProbeSource.CUSTOM, region.extraInfo());
        return RegionProbeTiling.fillRegionWithProbes(region.region(), source);
    }
}
