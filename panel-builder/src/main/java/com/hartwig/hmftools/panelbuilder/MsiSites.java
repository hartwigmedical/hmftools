package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.MSI_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.MSI_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.MSI_QUALITY_MIN;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// TODO: doc
public class MsiSites
{
    private static final TargetMetadata.Type TARGET_TYPE = TargetMetadata.Type.MSI;

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            MSI_QUALITY_MIN, MSI_GC_TARGET, MSI_GC_TOLERANCE);
    private static final ProbeSelector.Strategy PROBE_SELECT = new ProbeSelector.Strategy.MaxQuality();

    private static final Logger LOGGER = LogManager.getLogger(MsiSites.class);

    public static void generateProbes(final String msiSitesFile, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        LOGGER.info("Generating MSI probes");

        List<BasePosition> msiSites = loadMsiSitesFile(msiSitesFile);

        ProbeGenerationResult result = generateProbes(msiSites, probeGenerator, panelData);
        // Don't bother checking overlaps between MSI probes because there are few sites which are far apart.
        panelData.addResult(result);

        LOGGER.info("Done generating MSI probes");
    }

    private static List<BasePosition> loadMsiSitesFile(final String path)
    {
        LOGGER.debug("Loading MSI sites file: {}", path);

        try(DelimFileReader reader = new DelimFileReader(path))
        {
            int chrIndex = requireNonNull(reader.getColumnIndex(FLD_CHROMOSOME));
            int posIndex = requireNonNull(reader.getColumnIndex(FLD_POSITION));

            List<BasePosition> sites = reader.stream().map(row ->
            {
                String chromosome = row.get(chrIndex);
                int position = row.getInt(posIndex);
                return new BasePosition(chromosome, position);
            }).toList();

            LOGGER.info("Loaded {} MSI sites from {}", sites.size(), path);
            return sites;
        }
    }

    private static ProbeGenerationResult generateProbes(final List<BasePosition> msiSites, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        return msiSites.stream()
                .map(msiSite -> generateProbe(msiSite, probeGenerator, coverage))
                .reduce(new ProbeGenerationResult(), ProbeGenerationResult::add);
    }

    private static ProbeGenerationResult generateProbe(final BasePosition msiSite, final ProbeGenerator probeGenerator,
            final PanelCoverage coverage)
    {
        TargetMetadata metadata = createTargetMetadata(msiSite);
        ChrBaseRegion region = ChrBaseRegion.from(msiSite);
        return probeGenerator.coverRegion(region, metadata, PROBE_CRITERIA, PROBE_SELECT, coverage);
    }

    private static TargetMetadata createTargetMetadata(final BasePosition msiSite)
    {
        String extraInfo = msiSite.toString();
        return new TargetMetadata(TARGET_TYPE, extraInfo);
    }
}
